#include "TFP_SCP.h"

double get_wall_time_ilp(){
    struct timeval time;
    if(gettimeofday(&time,nullptr)){
        // HANDLE ERROR
        return 0;
    }else{
        return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_usec*0.000001); //microsegundos
    }
}

TFP_SCP::TFP_SCP(Reader *r, const char* typeSEC): rd(r){
    
    getcwd(CURRENT_DIR, 500);
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    current_day = tm.tm_mday; current_month = tm.tm_mon + 1; current_year = tm.tm_year + 1900;
    this->typeSEC = typeSEC;
    cout<< "[INFO] Graph type: " << rd->G_type << endl; 
    initILP(typeSEC);
}

TFP_SCP::~TFP_SCP(){
    env.end();
}


void TFP_SCP::initILP(const char* typeSEC){

    try{

        model = IloModel(env);
        initVariables();        

        if(typeSEC == "MTZ" || typeSEC == "mtz" || strcmp(typeSEC,"MTZ") == 0){
            cout << "[INFO] Creating Model TFP_SPC (SEC: MTZ)" << endl;
            createModel_MTZ(model, x, y, f, lambda,pi);
        }else if(typeSEC == "sign" || typeSEC == "SIGN" || strcmp(typeSEC,"SIGN") == 0){
            cout << "[INFO] Creating Model TFP_SPC (SEC: SIGN)" << endl; // constr ROSA
            createModel_SIGN(model, x, y, f, lambda,mu);
        }else{
            cout << "[ERROR] SEC not defined" << endl;
            cout << typeSEC << endl;
        }

        cplex = IloCplex(model);
        // exportILP(cplex,typeSEC);

    } catch (IloException& e) {
        cerr << "ERROR: " << e.getMessage()  << endl;
        cout << "\nError ilocplex" << endl;
        return;
    }catch (int e) {
            cerr << endl << "\nException occurred = " << e << endl;
    }


}

void TFP_SCP::initVariables(){

    // var x[u][j][s] : if worker u is in projecet j with skill s
    // x = BoolVar3Matrix(env,rd->num_vertices);
    x = BoolVar3Matrix(env,rd->num_vertices);
    for (int u=0; u<rd->num_vertices; u++){
        x[u] = BoolVarMatrix(env, rd->num_teams);
        for (int j = 0; j < rd->num_teams; j++)
            x[u][j] = IloBoolVarArray(env, rd->num_skills);
    }

    // var y[u][v][j] : if worker u and v is in the same project j
    y = BoolVar3Matrix(env,rd->num_vertices);
    for (int u=0; u<rd->num_vertices; u++){
        y[u] = BoolVarMatrix(env, rd->num_vertices);
        for (int v = u+1; v < rd->num_vertices; v++) // u<v
            y[u][v] = IloBoolVarArray(env, rd->num_teams);
    }

    // var f[u][v][p][q]: if arc pq is used in the flow associate to path(u,v)
    f = BoolVar4Matrix(env,rd->num_vertices);
    for (int u=0; u<rd->num_vertices; u++){
        f[u] = BoolVar3Matrix(env, rd->num_vertices);
        for (int v=u+1; v<rd->num_vertices; v++){ // u<v
            f[u][v] = BoolVarMatrix(env, rd->num_vertices);
            for (int p=0; p<rd->num_vertices; p++)
                f[u][v][p] = IloBoolVarArray(env, rd->num_vertices);
        }
    }


    // var lambda[u][v]
    lambda = IntVarMatrix(env,rd->num_vertices);
    for (int u = 0; u < rd->num_vertices; u++){
        lambda[u] = IloIntVarArray(env, rd->num_vertices, 0, rd->num_vertices*(rd->num_vertices-1)); // must put the bounds
    }
}


void TFP_SCP::createModel(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda){
    
    // add var x[u][j][s] in model
    for(int j=0; j<rd->num_teams; j++) // for all project j 
        for(int u=0; u< rd->num_vertices; u++) // for all worker u
            for(int s=0; s<rd->num_skills; s++){ 
                if(rd->K[u][s] > 0 && rd->R[j][s]>0){ // for all skills: s(u,j) skills of u and project j need skill u
                    char name[30];
                    sprintf(name, "x(%d,%d,%d)", u+1, j+1,s+1); // name not set all var beacause rd->K[u][s] > 0 && rd->R[j][s]>0 in definition have more vars
                    x[u][j][s].setName(name);
                    model.add(x[u][j][s]);
                }
            }

    // add var y[u][v][j]
    for(int u=0; u<rd->num_vertices; u++) // for all u 
        for(int v = u+1; v<rd->num_vertices; v++) // for all v: u<v
            for(int j = 0; j<rd->num_teams; j++){ // for all project j
                char name[30];
                sprintf(name, "y(%d,%d,%d)", u+1, v+1,j+1);
                y[u][v][j].setName(name);
                model.add(y[u][v][j]);
            }


    // add var f[u][v][p][q]
    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++) // for all u,v: u<v
            if(rd->G[u][v] == 0){ // and u,v not in E
                for(int p=0; p<rd->num_vertices; p++)
                    for(int q=0; q<rd->num_vertices; q++) // for all arc (p,q)
                        if(p!=q && rd->G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                            char name[30];
                            sprintf(name, "F(%d,%d,%d,%d)", u+1, v+1,p+1, q+1);
                            f[u][v][p][q].setName(name);
                            model.add(f[u][v][p][q]);
                        }
            }


    // add var lambda[u][v]
    for(int u = 0; u<rd->num_vertices; u++) //
        for(int v = u+1; v<rd->num_vertices; v++) // for all u,v : u<v
            if(rd->G[u][v] == 0){ // and u,v not in E
                char name[20];
                sprintf(name, "lambd(%d,%d)", u+1, v+1);
                lambda[u][v].setName(name);
                model.add(lambda[u][v]);
            }

    //object function
    objFunction(model,y,f);
    //constraints
    constr_OneTeam(model, x); // max one team with one skill
    constr_MinSkill(model, x); // min skill s per team j
    constr_LinY(model,x,y); // linearization y with x
    constr_NegEdge(model,y); // negative edge incompatibility
    constr_Flow(model,y,f); // flow const
    constr_flowSameTeam(model,y,f); // flow only if u,v is in the same team
    constr_PathComp(model,f,lambda); // every path with neg edges is pair
    //cuts
    // constr_Cuts_2(model,f); // not repeat vertice in path
}

void TFP_SCP::createModel_MTZ(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda, NumVar3Matrix mi){

    createModel(model, x, y, f, lambda);
    
//-------------------------MTZ VARIABLE---------------------------------------------
    IloEnv env = model.getEnv();
    // aloc var MTZ rd->r[u][v][i]
    pi = NumVar3Matrix(env,rd->num_vertices);
    for (int u=0; u<rd->num_vertices; u++){
        pi[u] = NumVarMatrix(env, rd->num_vertices);
        for (int v = u+1; v < rd->num_vertices; v++) // u<v
            pi[u][v] = IloNumVarArray(env, rd->num_vertices,0,rd->num_vertices, ILOFLOAT);
    }
    // add var MTZ rd->R[u][v][i]: position of vertice i in u,v-flow
    for(int u=0; u<rd->num_vertices; u++) // for all u 
        for(int v = u+1; v<rd->num_vertices; v++) // for all v: u<v
            for(int i = 0; i<rd->num_vertices; i++){ // for all vertice i
                char name[30];
                sprintf(name, "pi(%d,%d,%d)", u+1, v+1,i+1);
                pi[u][v][i].setName(name);
                model.add(pi[u][v][i]);
            }

    //constraints SEC
    constr_MTZ(model,f,pi); // break cycles using MTZ constraint

}


void TFP_SCP::createModel_SIGN(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda, BoolVar3Matrix mu){

    createModel(model, x, y, f, lambda);
    
//-------------------------SIGN VARIABLE---------------------------------------------
    IloEnv env = model.getEnv();
    // aloc var MTZ t[u][v][i]
    mu = BoolVar3Matrix(env,rd->num_vertices);
    for (int u=0; u<rd->num_vertices; u++){
        mu[u] = BoolVarMatrix(env, rd->num_vertices);
        for (int v = u+1; v < rd->num_vertices; v++) // u<v
            mu[u][v] = IloBoolVarArray(env,rd->num_vertices);
    }
    // add var SIGN t[u][v][i]: 1 if sign of the vertice p in u,v-path is positive, 0 otherwise
    for(int u=0; u<rd->num_vertices; u++) // for all u 
        for(int v = u+1; v<rd->num_vertices; v++) // for all v: u<v
            for(int i = 0; i<rd->num_vertices; i++){ // for all vertice i
                char name[30];
                sprintf(name, "mu(%d,%d,%d)", u+1, v+1,i+1);
                mu[u][v][i].setName(name);
                model.add(mu[u][v][i]);
            }

    //constraints SEC
    constr_SIGN(model,f,mu); // break cycles using SIGN constraint (rosa)

}



void TFP_SCP::objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);
    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++) // for all u,v : u<v
            for(int p=0; p<rd->num_vertices; p++)
                for(int q=0; q<rd->num_vertices; q++)
                    if(rd->G[u][v] == 0 && p!= q && rd->G[p][q] != 0){ // (u,v) not in E+ and (p,q) is a possible arc between (u,v)-path
                       objExpr += f[u][v][p][q];
                    }

    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++)
            for(int j = 0; j<rd->num_teams; j++){
                if(rd->G[u][v] >= 1){objExpr += y[u][v][j];}  // (u,v) in E+
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();
}

void TFP_SCP::constr_OneTeam(IloModel model, BoolVar3Matrix x){
    
    IloEnv env = model.getEnv();

    // each worker uses at most one skill 
    for (int u=0; u<rd->num_vertices; u++) {
        IloExpr expr(env);
        for (int j = 0; j <rd->num_teams; j++) {
            for(int s = 0; s<rd->num_skills; s++){
                if(rd->R[j][s]>0 && rd->K[u][s]>0) // team j need skill s and indidual u have skill s
                    expr += x[u][j][s];
            }
        }
        model.add(expr <= 1);   //var empty_sum to know if expr is empty or not, only add if not
        expr.end();
    }
}

void TFP_SCP::constr_MinSkill(IloModel model, BoolVar3Matrix x){
    IloEnv env = model.getEnv();

    // minimum number of individuals in project j with skill s
    for (int j=0; j <rd->num_teams; j++) {
        for (int s=0; s <rd->num_skills; s++) {
            IloExpr expr(env);
            if (rd->R[j][s] > 0) { //for all skill s in project j
                for (int u=0; u < rd->num_vertices; u++) {
                    if(rd->K[u][s] > 0) // if worker u have skill s
                        expr += x[u][j][s];
                }
                model.add(expr >= rd->R[j][s]);
                expr.end();
            }
        }
    }
}

void TFP_SCP::constr_LinY(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){

    IloEnv env = model.getEnv();

    // set var y with x (linearization)
    for (int u=0; u<rd->num_vertices; u++)
        for (int v = u+1; v < rd->num_vertices; v++) // for all u,v: u<v
            for (int j = 0; j <rd->num_teams; j++){
                IloExpr exprU(env);
                IloExpr exprV(env);
                for (int s = 0; s < rd->num_skills; s++) {
                    if(rd->K[u][s] > 0 && rd->R[j][s] > 0) // individual u have skill s and team j need skill s
                        exprU += x[u][j][s];
                    if(rd->K[v][s] > 0 && rd->R[j][s] > 0)
                        exprV += x[v][j][s];
                }
                model.add(exprU + exprV - y[u][v][j] <= 1);
                model.add(y[u][v][j]  <= exprU);
                model.add(y[u][v][j] <= exprV);
                exprU.end(); exprV.end();
            }

}

void TFP_SCP::constr_Flow(IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    
    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = u+1; v<rd->num_vertices; v++) 
            if(rd->G[u][v] == 0){                   // for all u<v and (u,v) not in E+
                for(int q = 0; q<rd->num_vertices; q++){ // for fixed vertice q
                    IloExpr sumArcs(env);
                    for(int p = 0; p<rd->num_vertices; p++) // for every arc (p,q) <=> (p,q) or (q,p) in E 
                        if(p!=q && rd->G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                            sumArcs += f[u][v][p][q];
                            sumArcs += -f[u][v][q][p];
                        }
                    if(q == u){
                        IloExpr expr(env);
                        for (int j = 0; j < rd->num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == -expr); expr.end();
                    }else if(q == v){
                        IloExpr expr(env);
                        for (int j = 0; j<rd->num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == expr); expr.end();
                    }else{
                        model.add(sumArcs == 0);
                    }
                }
            }

}

void TFP_SCP::constr_NegEdge(IloModel model, BoolVar3Matrix y){

    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = u+1; v<rd->num_vertices; v++)
            if(rd->G[u][v] == -1){
                IloExpr expr(env);
                for(int j =0; j<rd->num_teams; j++)
                    expr += y[u][v][j];
                model.add(expr == 0);
                expr.end();
            }

}

void TFP_SCP::constr_flowSameTeam(IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    // only use flow between u,v if u,v work in the same team
    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = u+1; v<rd->num_vertices; v++)
            if(rd->G[u][v] == 0){  // for all u,v in V: u<v and (u,v) not in E+
                for(int p=0; p<rd->num_vertices; p++)
                    for(int q=0; q<rd->num_vertices; q++)
                        if(p!=q && rd->G[p][q] != 0){ // (p,q) in A
                            IloExpr expr(env);
                            for (int j=0; j <rd->num_teams; j++) {
                                expr += y[u][v][j];
                            }
                            model.add(f[u][v][p][q] <= expr);
                            expr.end();
                        }
        }

}

void TFP_SCP::constr_PathComp(IloModel model, BoolVar4Matrix f, IntVarMatrix lambda){
    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = u+1; v< rd->num_vertices; v++){
            if(rd->G[u][v] == 0 ){ // (u,v) not in E+
                IloExpr expr(env);
                for(int p = 0; p<rd->num_vertices; p++)
                    for(int q = 0; q<rd->num_vertices; q++)
                        if(p!=q && rd->G[p][q] < 0){ // sum of negatives arcs(A-)
                            expr += f[u][v][p][q];
                        }
                model.add(expr - 2*lambda[u][v] == 0); // pair number of negatives edges
                expr.end();
            }
        }
}

void TFP_SCP::constr_MTZ(IloModel model, BoolVar4Matrix f, NumVar3Matrix pi){
    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices;u++)
        for(int v = u+1; v<rd->num_vertices;v++){
            for(int p=0;p<rd->num_vertices;p++)
                for(int q=0;q<rd->num_vertices;q++){
                    if(p!=q && abs(rd->G[p][q])){
                        IloExpr expr(env);
                        expr = pi[u][v][p] - pi[u][v][q] + rd->num_vertices*f[u][v][p][q];
                        model.add(expr <= (rd->num_vertices-1));
                        expr.end();
                    }
                }
        }
}

void TFP_SCP::constr_SIGN(IloModel model, BoolVar4Matrix f, BoolVar3Matrix mu){

    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices;u++)
        for(int v = u+1; v<rd->num_vertices;v++){
            for(int p=0;p<rd->num_vertices;p++)
                for(int q=0;q<rd->num_vertices;q++){
                    if(p!=q && rd->G[p][q] < 0){
                        IloExpr expr1(env); IloExpr expr2(env);
                        expr1 = f[u][v][p][q] - mu[u][v][p];
                        expr2 = 2 - f[u][v][p][q] - mu[u][v][p];
                        model.add(mu[u][v][q]  >= expr1);
                        model.add(mu[u][v][q] <= expr2);
                        expr1.end(); expr2.end();
                    }
                    if(p!=q && rd->G[p][q] > 0){
                        IloExpr expr1(env); IloExpr expr2(env);
                        expr1 = f[u][v][p][q] + mu[u][v][p] - 1;
                        expr2 = 1 - f[u][v][p][q] + mu[u][v][p];
                        model.add(mu[u][v][q] >= expr1);
                        model.add(mu[u][v][q] <= expr2);
                        expr1.end(); expr2.end();
                    } 
                }
    }
}

void TFP_SCP::printSolution(IloCplex& cplex,
               BoolVar3Matrix x,
               BoolVar3Matrix y,
               BoolVar4Matrix f,
               IntVarMatrix lambda){

    cout << "--------------------------------------------------------" << endl;

    cout << "Solution status = " << cplex.getStatus()   << endl;
    cout << "Solution value  = " << cplex.getObjValue() << endl;
    cout << "Valores das Variaveis " << endl;

    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;

    for(int j=0; j<rd->num_teams; j++)
        for(int u=0; u<rd->num_vertices; u++)
            for(int s=0; s<rd->num_skills; s++){
                if(rd->K[u][s] > 0 && rd->R[j][s]>0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        cout << "\t " << "x"<< "["<< u+1  << "," << j+1 << "," << s+1 << "]" << " -> " << " 1 \n";
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++)
            for(int j=0; j<rd->num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    cout << "\t " << "y"<< "["<< u+1  << "," << v+1 << "," << j+1 << "]" << " -> " << " 1 \n";
                }
            }

    cout << "--------------------------------------------------------" << endl;

    // for(int u = 0; u<rd->num_vertices; u++)
    //     for(int v = 0; v<rd->num_vertices; v++)
    //         if(u<v && rd->G[u][v] == 0){
    //             for(int p = 0; p<rd->num_vertices; p++)
    //                 for(int q = 0; q<rd->num_vertices; q++)
    //                     if(p!= q && rd->G[p][q] != 0){
    //                         if(cplex.getValue(f[u][v][p][q]) > 0){
    //                             cout << "\t " << "f"<< "["<< u+1  << "," << v+1 << "," << p+1 << ", " << q+1 << "]" << " -> " << " 1 \n";
    //                         }
    //                     }
    //         }           


}

void TFP_SCP::saveSolution(IloCplex& cplex,
                    BoolVar3Matrix x,
                    BoolVar3Matrix y,
                    BoolVar4Matrix f,
                    int class_type){

    char arq[1000];
    char arqv_instance[50];
    sprintf(arqv_instance, "%s_%d", instanceG,class_type);
    // sprintf(arq, "%s/results/%s.txt",CURRENT_DIR, arqv_instance);
    sprintf(arq, "%s/results/solution_%d_Vertices_%d-%d-%d_TFP_SCP_%c.ods",CURRENT_DIR, rd->num_vertices, current_year, current_month, current_day,typeSEC);

    FILE *file = fopen(arq, "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(file, "Solution value = %.1f \n", cplex.getObjValue());

    for(int j=0; j<rd->num_teams; j++)
        for(int u=0; u<rd->num_vertices; u++)
            for(int s=0; s<rd->num_skills; s++){
                if(rd->K[u][s] > 0 && rd->R[j][s]>0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        fprintf(file, "x[%d,%d,%d] = 1 \n", u+1, j+1, s+1);
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++)
            for(int j=0; j<rd->num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    fprintf(file, "y[%d,%d,%d] = 1 \n", u+1, v+1, j+1);
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = 0; v<rd->num_vertices; v++)
            if(u<v && rd->G[u][v] == 0){
                for(int p = 0; p<rd->num_vertices; p++)
                    for(int q = 0; q<rd->num_vertices; q++)
                        if(p!= q && rd->G[p][q] != 0){
                            if(cplex.getValue(f[u][v][p][q]) > 0){
                                fprintf(file, "f[%d,%d,%d,%d] = 1 \n", u+1, v+1, p+1,q+1);
                            }
                        }
            }


    fclose(file);                
}


// void TFP_SCP::saveResults(IloCplex& cplex,
//                    double timeF){
void TFP_SCP::saveResults(double timeTotal){

    char arq[1000];
    // sprintf(arq, "%s/results/%d_Vertices_2022-06-27_cycle_mtz.ods",CURRENT_DIR, rd->num_vertices);
    sprintf(arq, "%s/results/result_%d_Vertices_%d-%d-%d_TFP_SCP_%s.ods",CURRENT_DIR, rd->num_vertices, current_year, current_month, current_day,typeSEC);

    ofstream outputTable;
    outputTable.open(arq,ios:: app);
    if(outputTable.is_open()){

        outputTable << rd->instanceG << ";"; // grafo instancia
        outputTable << rd->instanceKR << ";"; // KR instancia
        outputTable << rd->num_vertices << ";";   // numero de vertices
        outputTable << rd->num_skills << ";";   // qtd de hab
        outputTable << rd->num_teams << ";";   // qtd de proj
        outputTable << cplex.getStatus() << ";"; // Status cplex
        outputTable << cplex.getObjValue() << ";"; // valor fo
        outputTable << cplex.getNnodes() << ";"; // num nos
        outputTable << cplex.getMIPRelativeGap() <<";"; // gap
        outputTable << timeTFP_SCP <<  ";"; // tempo execucao flp
        outputTable << cplex.getTime() <<  ";"; // tempo execucao cplex
        outputTable << timeTotal <<  ";"; // tempo execucao cplex
        outputTable << " \n ";


    }

    outputTable.close();




}

void TFP_SCP::printInstance(){
    rd->show();
}


void TFP_SCP::exportILP(IloCplex& cplex,const char *typeSEC){
    char ilpname[40];
    for(int i=0 ; i < 40 ; i++)
        ilpname[i] = 0;
    strcat(ilpname, "TFP_SCP-");
    strcat(ilpname, typeSEC);
    strcat(ilpname, ".lp");
    cout<< "[INFO]: Creating the the file lp: ";
    cout << ilpname << endl;
    cplex.exportModel(ilpname);
}

void TFP_SCP::solveILP(){
    double cpu0, cpu1;
    cpu0 = get_wall_time_ilp(); 
    if (!cplex.solve()){
        cpu1 = get_wall_time_ilp();
        this->timeTFP_SCP = cpu1 - cpu0;
        env.error() << "Failed to optimize LP." << endl;
        throw(-1);
    }
    // else{
    //     printSolution(cplex,x,y,f,lambda);
    // }
    cpu1 = get_wall_time_ilp();
    this->timeTFP_SCP = cpu1 - cpu0; 
}

