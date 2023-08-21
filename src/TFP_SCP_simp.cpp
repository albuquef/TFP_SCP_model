#include "TFP_SCP_simp.h"

double get_wall_time(){
    struct timeval time;
    if(gettimeofday(&time,nullptr)){
        // HANDLE ERROR
        return 0;
    }else{
        return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_usec*0.000001); //microsegundos
    }
}

TFP_SCP_SIMP::TFP_SCP_SIMP(Reader *r, const Graph & Grph,const char* method_SPP): rd(r), GRAPH(Grph){
// TFP_SCP_SIMP::TESTE(){
    
    getcwd(CURRENT_DIR, 500);
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    current_day = tm.tm_mday; current_month = tm.tm_mon + 1; current_year = tm.tm_year + 1900;
    this->method_SPP = method_SPP;

    // cout<< "[INFO] Graph type: " << rd->G_type << endl; 
    initModel(method_SPP);



}

TFP_SCP_SIMP::~TFP_SCP_SIMP(){
    env.end();
}


void TFP_SCP_SIMP::initModel(const char* method_SPP){

    try{

        model = IloModel(env);
        initVariables();        

        if(method_SPP == "DijkstraComp" || method_SPP == "DijkstraTwoLabels" || strcmp(method_SPP,"DijkstraComp") == 0){
            cout << "[INFO] Creating Weighted Graph Shortest Positive Path " << endl;
            cout << "[INFO] Method: "<< method_SPP << endl;
            cout << "[WARNING] Signed Graph must be balanced " << endl; // todo: add function to check if the graph is balanced


            double cpu0, cpu1;
            cpu0 = get_wall_time(); 
            
            GRAPH.genereateGraph_weighted("DijkstraComp");
            
            cpu1 = get_wall_time();
            this->time_GraphSPP = cpu1 - cpu0;

            createModel(model,x,y);

        }else if(method_SPP == "Matching" || method_SPP == "MinMatching" || strcmp(method_SPP,"MinMatching") == 0){
            
            cout << "[INFO] Creating Weighted Graph Shortest Positive Path " << endl;
            cout << "[INFO] Method: "<< method_SPP << endl;
            cout << "[WARNING] Signed Graph must be undirected " << endl; // todo: check if graph is not directed

            double cpu0, cpu1;
            cpu0 = get_wall_time(); 
            
            GRAPH.genereateGraph_weighted("MinMatching");
            
            cpu1 = get_wall_time();
            this->time_GraphSPP = cpu1 - cpu0;

            createModel(model,x,y);


        }else if(method_SPP == "SEC_MTZ" || method_SPP == "MTZ" || strcmp(method_SPP,"MTZ") == 0){
            cout << "[INFO] Creating Weighted Graph Shortest Positive Path " << endl;
            cout << "[INFO] ILP Method: "<< method_SPP << endl;

            // cout << "not working yet" << endl;
            
            create_GraphSPP_ILP();
            createModel(model,x,y);
            
        }else if(method_SPP == "SEC_SIGN" || method_SPP == "SIGN" || strcmp(method_SPP,"SIGN") == 0){
            cout << "[INFO] Creating Weighted Graph Shortest Positive Path " << endl;
            cout << "[INFO] ILP Method: "<< method_SPP << endl;

            cout << "not working yet" << endl;
            create_GraphSPP_ILP();
            // create_GraphSPP_MTZ(model, x, y, f, lambda,mu);

            createModel(model,x,y);
        }
        else{
            cout << "[ERROR] Method for Shortest Positive Paths not defined" << endl;
        }

        cplex = IloCplex(model);
        // exportILP(cplex,method_SPP);

    } catch (IloException& e) {
        cerr << "ERROR: " << e.getMessage()  << endl;
        cout << "\nError ilocplex" << endl;
        return;
    }catch (int e) {
            cerr << endl << "\nException occurred = " << e << endl;
    }


}


void TFP_SCP_SIMP::allocVars_SPP(IloEnv env_SPP,BoolVarMatrix f){
    // var f(u,v)[p][q]: if arc pq is used in the flow associate to path(u,v)
    for (int p=0; p<rd->num_vertices; p++)
        f[p] = IloBoolVarArray(env_SPP, rd->num_vertices);
}
void TFP_SCP_SIMP::objFunction_SPP_uv(IloModel model_SPP,BoolVarMatrix f, int u, int v){
    IloEnv env_spp = model_SPP.getEnv();
    IloExpr objExpr(env_spp);
    for(int p=0; p<rd->num_vertices; p++)
        for(int q=0; q<rd->num_vertices; q++)
            if(rd->G[u][v] < 1 && p!= q && rd->G[p][q] != 0){ // (u,v) not in E+ and (p,q) is a possible arc between (u,v)-path
                objExpr += f[p][q];
            }

    IloObjective obj = IloMinimize(env_spp, objExpr);
    model_SPP.add(obj);
    objExpr.end();
}
void TFP_SCP_SIMP::constr_Flow_uv(IloModel model_SPP,BoolVarMatrix f, int u, int v){
    IloEnv env_spp = model_SPP.getEnv();

    if(rd->G[u][v] < 1){  // for all u<v and (u,v) not in E+
        for(int q = 0; q<rd->num_vertices; q++){ // for fixed vertice q
            IloExpr sumArcs(env_spp);
            for(int p = 0; p<rd->num_vertices; p++) // for every arc (p,q) <=> (p,q) or (q,p) in E 
                if(p!=q && rd->G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                    sumArcs += f[p][q];
                    sumArcs += -f[q][p];
                }
            if(q == u){
                IloExpr expr(env_spp);
                model_SPP.add(sumArcs == -1); expr.end();
            }else if(q == v){
                IloExpr expr(env_spp);
                model_SPP.add(sumArcs == 1); expr.end();
            }else{
                model_SPP.add(sumArcs == 0);
            }
        }
    }
}
void TFP_SCP_SIMP::constr_PathComp_uv(IloModel model_SPP, BoolVarMatrix f, IloIntVar lambda, int u, int v){

    IloEnv env_spp = model_SPP.getEnv();
    
    if(rd->G[u][v] < 1 ){ // (u,v) not in E+
        IloExpr expr(env_spp);
        for(int p = 0; p<rd->num_vertices; p++)
            for(int q = 0; q<rd->num_vertices; q++)
                if(p!=q && rd->G[p][q] < 0){ // sum of negatives arcs(A-)
                    expr += f[p][q];
                }
        model_SPP.add(expr - 2*lambda == 0); // pair number of negatives edges
        expr.end();
    }
}
void TFP_SCP_SIMP::constr_BreakCycle_MTZ_uv(IloModel model_SPP, BoolVarMatrix f, IloNumVarArray pi, int u, int v){

    IloEnv env_spp = model_SPP.getEnv();
    if(rd->G[u][v] < 1 ){ // (u,v) not in E+

        for(int p=0;p<rd->num_vertices;p++)
            for(int q=0;q<rd->num_vertices;q++){
                if(p!=q && abs(rd->G[p][q])){
                    IloExpr expr(env_spp);
                    expr = pi[p] - pi[q] + rd->num_vertices*f[p][q];
                    model_SPP.add(expr <= (rd->num_vertices-1));
                    expr.end();
                }
            }
    }
    

}
void TFP_SCP_SIMP::constr_BreakCycle_SIGN_uv(IloModel model_SPP, BoolVarMatrix f, IloBoolVarArray mu, int u, int v){

    IloEnv env_spp = model_SPP.getEnv();

    for(int p=0;p<rd->num_vertices;p++)
        for(int q=0;q<rd->num_vertices;q++){
            if(p!=q && rd->G[p][q] < 0){
                IloExpr expr1(env_spp); IloExpr expr2(env_spp);
                expr1 = f[p][q] - mu[p];
                expr2 = 2 - f[p][q] - mu[p];
                model_SPP.add(mu[q]  >= expr1);
                model_SPP.add(mu[q] <= expr2);
                expr1.end(); expr2.end();
            }
            if(p!=q && rd->G[p][q] > 0){
                IloExpr expr1(env_spp); IloExpr expr2(env_spp);
                expr1 = f[p][q] + mu[p] - 1;
                expr2 = 1 - f[p][q] + mu[p];
                model_SPP.add(mu[q] >= expr1);
                model_SPP.add(mu[q] <= expr2);
                expr1.end(); expr2.end();
            } 
        }


}

void TFP_SCP_SIMP::createModel_MTZ(IloModel model_MTZ,BoolVarMatrix f,IloIntVar lambda,IloNumVarArray pi,int u, int v){
        // add var f(u,v)[p][q]
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<rd->num_vertices; p++)
            for(int q=0; q<rd->num_vertices; q++) // for all arc (p,q)
                if(p!=q && rd->G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                    char name[50];
                    sprintf(name, "F(%d,%d)%d%d", u+1, v+1,p+1, q+1);
                    f[p][q].setName(name);
                    model_MTZ.add(f[p][q]);
                }
    }

    // add var lambda(u,v)
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        char name[20];
        sprintf(name, "lambd(%d,%d)", u+1, v+1);
        lambda.setName(name);
        model_MTZ.add(lambda);
    }


    // add var r(u,v)[p]
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<rd->num_vertices; p++){
            char name[50];
            sprintf(name, "pi(%d,%d)%d", u+1, v+1,p+1);
            pi[p].setName(name);
            model_MTZ.add(pi[p]);
        }
    }

    objFunction_SPP_uv (model_MTZ,f,u,v);
    constr_Flow_uv(model_MTZ,f, u, v); // flow const 
    constr_PathComp_uv(model_MTZ,f,lambda, u,v); // every path with neg edges is pair
    constr_BreakCycle_MTZ_uv(model_MTZ,f,pi,u,v); // break cycle using MTZ
    
}
void TFP_SCP_SIMP::createModel_SIGN(IloModel model_SIGN,BoolVarMatrix f,IloIntVar lambda,IloBoolVarArray mu,int u, int v){
            // add var f(u,v)[p][q]
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<rd->num_vertices; p++)
            for(int q=0; q<rd->num_vertices; q++) // for all arc (p,q)
                if(p!=q && rd->G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                    char name[50];
                    sprintf(name, "F(%d,%d)%d%d", u+1, v+1,p+1, q+1);
                    f[p][q].setName(name);
                    model_SIGN.add(f[p][q]);
                }
    }

    // add var lambda(u,v)
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        char name[20];
        sprintf(name, "lambd(%d,%d)", u+1, v+1);
        lambda.setName(name);
        model_SIGN.add(lambda);
    }


    // add var r(u,v)[p]
    if(rd->G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<rd->num_vertices; p++){
            char name[50];
            sprintf(name, "mu(%d,%d)%d", u+1, v+1,p+1);
            mu[p].setName(name);
            model_SIGN.add(mu[p]);
        }
    }

    objFunction_SPP_uv (model_SIGN,f,u,v);
    constr_Flow_uv(model_SIGN,f, u, v); // flow const 
    constr_PathComp_uv(model_SIGN,f,lambda, u,v); // every path with neg edges is pair
    constr_BreakCycle_SIGN_uv(model_SIGN,f,mu,u,v); // break cycle using SIGN
}
void TFP_SCP_SIMP::create_GraphSPP_ILP(){


    GRAPH.Graph_SPP = (int**)(malloc(rd->num_vertices*sizeof(int*)));
    for (int u=0;u<rd->num_vertices;u++)
    {
        GRAPH.Graph_SPP[u] = (int*)malloc(rd->num_vertices*sizeof(int));
    }
    

    if (strcmp(rd->G_type,"directed") == 0){

        // int G_aux[]

        // for(int u=0; u<rd->num_vertices;u++)
        //     for(int v=0; v<rd->num_vertices; v++){
        //         if(rd->G[u][v] == 0){

        //             IloEnv env_SPP;
        //             try{
        //                 IloModel model_SPP(env_SPP);
                        
        //                 BoolVarMatrix f(env_SPP,rd->num_vertices);
        //                 IloIntVar lambda(env_SPP);
        //                 if(strcmp(method_SPP,"MTZ") == 0){
        //                     IloNumVarArray pi(env_SPP, rd->num_vertices,0,rd->num_vertices, ILOFLOAT);
        //                     allocVars_SPP(env_SPP,f);
        //                     createModel_MTZ(model_SPP,f,lambda,pi,u,v);       // can optmize var f, put all dif u,v = zero or better use 2 indices, only pq
        //                 }
        //                 if(strcmp(method_SPP,"SIGN") == 0){
        //                     IloBoolVarArray mu(env_SPP,rd->num_vertices);
        //                     allocVars_SPP(env_SPP,f);
        //                     createModel_SIGN(model_SPP,f,lambda,mu,u,v);
        //                 }
                        
        //                 IloCplex cplex_spp;
        //                 cplex_spp = IloCplex(model_SPP);
        //                 cplex_spp.setOut(env_SPP.getNullStream());
                        
        //                 if (!cplex_spp.solve()){
        //                     // infeasible -> infinite
        //                     GRAPH.Graph_SPP[u][v] = INF; 
        //                     throw(-1);
        //                 }


        //                 GRAPH.Graph_SPP[u][v] = int(cplex_spp.getObjValue());

        //                 // cout << "-------------------------------------------------- \n\n" << endl;
        //             }catch (IloException& e) {
        //                     cerr << "Concert exception caught: " << e << endl;
        //                 }
        //                 catch (...) {
        //                     // cerr << "Unknown exception caught" << endl;
        //                 }

        //             env_SPP.end();


        //         }if(rd->G[u][v] == 1){
        //             GRAPH.Graph_SPP[u][v] = 1;
        //         }if(rd->G[u][v] == -1){
        //             // incompatible -> infinite
        //             GRAPH.Graph_SPP[u][v] = INF;
        //         }

        //     }

        // // for(int u=0; u<rd->num_vertices;u++)
        // //     for(int v=u+1; v<rd->num_vertices; v++){
            





        // // }
    
    }
    else{
        for(int u=0; u<rd->num_vertices;u++)
            for(int v=u+1; v<rd->num_vertices; v++){
                if(rd->G[u][v] == 0){

                    IloEnv env_SPP;
                    try{
                        IloModel model_SPP(env_SPP);
                        
                        BoolVarMatrix f(env_SPP,rd->num_vertices);
                        IloIntVar lambda(env_SPP);
                        if(strcmp(method_SPP,"MTZ") == 0){
                            IloNumVarArray pi(env_SPP, rd->num_vertices,0,rd->num_vertices, ILOFLOAT);
                            allocVars_SPP(env_SPP,f);
                            createModel_MTZ(model_SPP,f,lambda,pi,u,v);       // can optmize var f, put all dif u,v = zero or better use 2 indices, only pq
                        }
                        if(strcmp(method_SPP,"SIGN") == 0){
                            IloBoolVarArray mu(env_SPP,rd->num_vertices);
                            allocVars_SPP(env_SPP,f);
                            createModel_SIGN(model_SPP,f,lambda,mu,u,v);
                        }
                        
                        IloCplex cplex_spp;
                        cplex_spp = IloCplex(model_SPP);
                        cplex_spp.setOut(env_SPP.getNullStream());
                        
                        if (!cplex_spp.solve()){
                            // infeasible -> infinite
                            GRAPH.Graph_SPP[u][v] = INF; 
                            GRAPH.Graph_SPP[v][u] = INF;
                            throw(-1);
                        }


                        GRAPH.Graph_SPP[u][v] = int(cplex_spp.getObjValue());
                        GRAPH.Graph_SPP[v][u] = int(cplex_spp.getObjValue());

                        // cout << "-------------------------------------------------- \n\n" << endl;
                    }catch (IloException& e) {
                            cerr << "Concert exception caught: " << e << endl;
                        }
                        catch (...) {
                            // cerr << "Unknown exception caught" << endl;
                        }

                    env_SPP.end();


                }if(rd->G[u][v] == 1){
                    GRAPH.Graph_SPP[u][v] = 1;
                    GRAPH.Graph_SPP[v][u] = 1;
                }if(rd->G[u][v] == -1){
                    // incompatible -> infinite
                    GRAPH.Graph_SPP[u][v] = INF;
                    GRAPH.Graph_SPP[v][u] = INF;
                }

            }
    }


}

void TFP_SCP_SIMP::initVariables(){

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

}

void TFP_SCP_SIMP::createModel(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){
    
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

    //object function
    objFunction(model,y);
    //constraints
    constr_OneTeam(model, x); // max one team with one skill
    constr_MinSkill(model, x); // min skill s per team j
    constr_LinY(model,x,y); // linearization y with x
    constr_Incomp(model,y); // negative edge incompatibility
}

void TFP_SCP_SIMP::objFunction (IloModel model, BoolVar3Matrix y){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);

    for(int u=0; u<rd->num_vertices; u++)
        for(int v=u+1; v<rd->num_vertices; v++)
            for(int j = 0; j<rd->num_teams; j++){
                if(rd->G[u][v] >= 0){objExpr += GRAPH.Graph_SPP[u][v]*y[u][v][j];}  // (u,v) not in E- or (u,v) in E'
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();
}

void TFP_SCP_SIMP::constr_OneTeam(IloModel model, BoolVar3Matrix x){
    
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

void TFP_SCP_SIMP::constr_MinSkill(IloModel model, BoolVar3Matrix x){
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

void TFP_SCP_SIMP::constr_LinY(IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){

    IloEnv env = model.getEnv();

    // set var y with x (linearization)
    for (int u=0; u<rd->num_vertices; u++)
        for (int v = u+1; v < rd->num_vertices; v++) // for all u,v: u<v
            // if (rd->G[u][v] >= 0  && GRAPH.Graph_SPP[u][v]>0){
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
            // }
}


// can change to constranint with variables x
void TFP_SCP_SIMP::constr_Incomp(IloModel model, BoolVar3Matrix y){

    IloEnv env = model.getEnv();

    for(int u = 0; u<rd->num_vertices; u++)
        for(int v = u+1; v<rd->num_vertices; v++)
            if(rd->G[u][v] == -1 || GRAPH.Graph_SPP[u][v]==0){
                IloExpr expr(env);
                for(int j =0; j<rd->num_teams; j++)
                    expr += y[u][v][j];
                model.add(expr == 0);
                expr.end();
            }

}


void TFP_SCP_SIMP::printSolution(IloCplex& cplex,
                            BoolVar3Matrix x,
                            BoolVar3Matrix y){
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

} 

void TFP_SCP_SIMP::saveSolution(IloCplex& cplex,
                    BoolVar3Matrix x,
                    BoolVar3Matrix y,
                    int class_type){

    char arq[1000];
    char arqv_instance[50];
    sprintf(arqv_instance, "%s_%d", instanceG,class_type);
    // sprintf(arq, "%s/results/%s.txt",CURRENT_DIR, arqv_instance);
    sprintf(arq, "%s/results/%d_Vertices_%d-%d-%d_TFP_SCP_SIMP_%c.csv",CURRENT_DIR, rd->num_vertices, current_year, current_month, current_day,method_SPP);

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

    fclose(file);                
}

// void TFP_SCP_SIMP::saveResults(IloCplex& cplex,
//                    double timeTFP, double time_WeightedGraph){
void TFP_SCP_SIMP::saveResults(double timeTotal){

    char arq[1000];
    // sprintf(arq, "%s/results/%d_Vertices_2022-06-27_cycle_mtz.ods",CURRENT_DIR, rd->num_vertices);
    sprintf(arq, "%s/results/result_%d_Vertices_%d-%d-%d_TFP_SCP_SIMP_%s.ods",CURRENT_DIR, rd->num_vertices, current_year, current_month, current_day,method_SPP);

    cout << arq << endl;

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
        outputTable << time_GraphSPP <<  ";"; // tempo execucao Grsph_SPP
        outputTable << timeTFP <<  ";"; // tempo execucao tfp (felipe)
        outputTable << cplex.getTime() <<  ";"; // tempo execucao tfp (cplex)
        outputTable << timeTotal <<  ";"; // tempo execucao tfp (cplex)
        outputTable << " \n ";


    }

    outputTable.close();




}


void TFP_SCP_SIMP::printInstance(){
    rd->show();
}


void TFP_SCP_SIMP::exportILP(IloCplex& cplex,const char* method_SPP){
    char ilpname[50];
    for(int i=0 ; i < 50 ; i++)
        ilpname[i] = 0;
    // string ilpname;
    strcat(ilpname, "TFP_SCP_simp-");
    strcat(ilpname, method_SPP);
    strcat(ilpname, ".lp");
    cout<< "[INFO]: Creating the the file lp: ";
    cout << ilpname << endl;
    cplex.exportModel(ilpname);
}

void TFP_SCP_SIMP::solveILP(){
    
    double cpu0, cpu1;
    cpu0 = get_wall_time(); 
    if (!cplex.solve()){
        env.error() << "Failed to optimize LP." << endl;
        // throw(-1);
    }
    // else{
        // cout << "--------------------------------------------------------" << end;
        // cout << "Solution status = " << cplex.getStatus()   << endl;
        // cout << "Solution value  = " << cplex.getObjValue() << endl;
        // cout << "--------------------------------------------------------" << endl;
    // }
    cpu1 = get_wall_time();
    this->timeTFP = cpu1 - cpu0; 
}

