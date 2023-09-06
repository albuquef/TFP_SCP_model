//#include <iostream>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <time.h>
#include <sys/time.h>
#include <string>
#include <filesystem>
#include <list>
#include<queue>

#include<fstream>
using namespace std;

// ********** VARIAVEIS GLOBAIS **********
int num_vertices;/* numero de vertices/individuos do grafo*/
int num_skills;/*numero de habilidades*/
int num_teams; /* quantidade de total de equipes/projetos*/

int **G; /*matriz |n|x|n| grafo*/
int **K; /*matriz |n|x|f| de pessoa - habilidade*/
int **R; /* matriz demanda |m|x|f|: quantidade de cada pessoa com habilidade k_a em cada equipe*/

char instanceG[50];
char instanceKR[50];

char CURRENT_DIR[500];

int current_day, current_month, current_year;

//double runTime;
#define BILLION 1000000000L


//typedef
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3Matrix;


//--------------------------------------------------------------------------------------------------------------
double get_wall_time(){
    struct timeval time;
    if(gettimeofday(&time,nullptr)){
        // HANDLE ERROR
        return 0;
    }else{
        return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_usec*0.000001); //microsegundos
    }
}

//--------------------------------------------read data------------------------------------------------------------------
//load matrix G
void load_Graph(char arq[], bool show, string instance_type){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int edge_weight;
    fscanf(inst,"%d",&num_vertices);
    G = (int**)(malloc(num_vertices*sizeof(int*)));

    for (int u=0;u<num_vertices;u++)
    {
        G[u] = (int*)malloc(num_vertices*sizeof(int));

        for (int v=0;v<num_vertices;v++)
        {
            fscanf(inst,"%d",&edge_weight);
            // G[u][v] = edge_weight;
            if(u == v){
                G[u][v] = 1;
            }else{
                if(edge_weight>0)
                    G[u][v] = 1;
                else if(edge_weight<0)
                    G[u][v] = -1;
                else if(edge_weight==0)
                    G[u][v] = 0;
            }
        }

    }
    fclose(inst);


    if (instance_type == "random"){
        // cout << "random" << endl;
        for (int u = 0; u < num_vertices; u++){
            for (int v = u+1; v < num_vertices; v++){
                if (G[u][v] == 1 && G[v][u] != 1){
                    G[v][u] = 1;
                }else if (G[v][u] == 1 && G[u][v] != 1){
                    G[u][v] = 1;
                }else if (G[u][v] == -1 && G[v][u] != -1){
                    G[v][u] = -1;
                }else if (G[v][u] == -1 && G[u][v] != -1){
                    G[u][v] = -1;
                }
            }
        }
    }


    if(show){
        cout << endl << "Numero de individuos: " << num_vertices << endl;
        cout << " Adjacency Matrix \n";
        for (int u = 0; u < num_vertices; u++){
            for (int v = 0; v < num_vertices; v++)
                cout <<  G[u][v] << " ";
            cout<< "\n";
        }
    }
}

//load individuals skills
void load_K(char arq[], bool show){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int value;
    fscanf(inst,"%d",&num_skills);
    K = (int**)(malloc(num_vertices*sizeof(int*)));
    for (int u=0;u<num_vertices;u++){
        K[u] = (int*)malloc(num_skills*sizeof(int));

        for (int s=0;s<num_skills;s++){
            fscanf(inst,"%d",&value);
            K[u][s] = value;
        }
    }
    if(show){
        printf("\nHabilidades = %d \n", num_skills);
        for (int u=0;u<num_vertices;u++){
            for (int s=0;s<num_skills;s++)
                printf("%d ", K[u][s]);
            printf("\n");
        }
    }

    fclose(inst);
}

//load the skills necessary in projects/teams
void load_R(char arq[], bool show){

    FILE *inst = fopen(arq,"r");
    if(!inst)
        printf("erro na leitura do arquivo!");

    int value;
    fscanf(inst,"%d",&num_teams);
    R = (int**)(malloc(num_teams*sizeof(int*)));

    for (int j=0;j<num_teams;j++){
        R[j] = (int*)malloc(num_skills*sizeof(int));

        for (int s=0;s<num_skills;s++){
            fscanf(inst,"%d",&value);
            R[j][s] = value;
        }
    }
    if(show){
        printf("\nEquipes = %d\n", num_teams);

        for (int j=0;j<num_teams;j++){
            for (int s=0;s<num_skills;s++)
                printf("%d ", R[j][s]);
            printf("\n");
        }
    }

    fclose(inst);
}


static void
instanceReader(bool show, int vert, int graph_class, int class_type, string instance_type){
    char arq1[600]; char arq2[600]; char arq3[600];
    char arq_base[500];
    sprintf(arq_base, "%s/instances_COR", CURRENT_DIR);

    if (instance_type == "random"){
        sprintf(instanceG, "%dverticesS%d", vert, graph_class);
    }
    else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
        sprintf(instanceG, "%dvertices_%s_S%d", vert, instance_type.c_str(), graph_class);

    }else{
        cout << "[INFO] not a valid instance type" << endl;
    }


    sprintf(arq1, "%s/%dVertices/%s.txt", arq_base, vert, instanceG);
    cout << "[INFO] Read graph: "<< instanceG << endl;

    sprintf(instanceKR, "%dVertices/class1/%d", vert, class_type);
    sprintf(arq2, "%s/%s/K.txt", arq_base, instanceKR);
    sprintf(arq3, "%s/%s/R.txt", arq_base, instanceKR);
    cout << "[INFO] Instance class type: "<< instanceKR << endl;


    load_Graph(arq1, show, instance_type); load_K(arq2, show); load_R(arq3, show);

}

static void
graphReader(bool show, int vert, int graph_class, string instance_type){
    char arq1[600]; 
    char arq_base[500];
    sprintf(arq_base, "%s/instances_COR", CURRENT_DIR);

    if (instance_type == "random"){
        sprintf(instanceG, "%dverticesS%d", vert, graph_class);
    }
    else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
        sprintf(instanceG, "%dvertices_%s_S%d", vert, instance_type.c_str(), graph_class);

    }else{
        cout << "[INFO] not a valid instance type" << endl;
    }


    sprintf(arq1, "%s/%dVertices/%s.txt", arq_base, vert, instanceG);

    load_Graph(arq1, show, instance_type); 
}

static void
    instanceReader (bool show, int vert, int graph_class, int class_type, string instance_type);

static void
   allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda);

static void
   createModel (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda);

static void
   createModel_MTZ (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda, NumVar3Matrix r);

static void
   objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest_OneTeam (IloModel model, BoolVar3Matrix x);

static void
   rest_MinSkill (IloModel model, BoolVar3Matrix x);

static void
   rest_LinY (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   rest_Flow (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest_flowSameTeam (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

static void
   rest_PathComp (IloModel model, BoolVar4Matrix f, IntVarMatrix lambda);

static void
   rest_NegEdge (IloModel model, BoolVar3Matrix y);

static void
   restCuts_1 (IloModel model, BoolVar4Matrix f);

static void
   restCuts_2 (IloModel model, BoolVar4Matrix f);

static void
   restCycle_MTZ (IloModel model, BoolVar4Matrix f, NumVar3Matrix r);


static void
    displaySolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y,
                   BoolVar4Matrix f,
                   IntVarMatrix lambda);
static void
    saveSolution(IloCplex& cplex,
                   BoolVar3Matrix X,
                   BoolVar3Matrix Y,
                   BoolVar4Matrix F,
                   int class_type);
static void
    saveResults(IloCplex& cplex,
                    double timeF);


static void
    savePath(IloCplex& cplex,
                   BoolVar4Matrix F,
                   int class_type);

static void
    testeGrafos();

void traverse(int u, bool visited[]) {
   visited[u] = true; //mark v as visited
   for(int v = 0; v<num_vertices; v++) {
      if(G[u][v]) {
         if(!visited[v])
            traverse(v, visited);
      }
   }
}
bool isConnected() {
   bool *vis = new bool[num_vertices];
   //for all vertex u as start point, check whether all nodes are visible or not
   for(int u = 0; u < num_vertices; u++) {
      for(int i = 0; i<num_vertices; i++)
         vis[i] = false; //initialize as no node is visited
         traverse(u, vis);
      for(int i = 0; i<num_vertices; i++) {
         if(!vis[i]) //if there is a node, not visited by traversal, graph is not connected
            return false;
      }
   }
   return true;
}

static void
    runTests(string SEC_type);

int main()
{
    getcwd(CURRENT_DIR, 500);

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    current_day = tm.tm_mday; current_month = tm.tm_mon + 1; current_year = tm.tm_year + 1900;

    // testeGrafos();  
    // dt nao ta funcionando
    string sec_type = "mtz";
    runTests(sec_type);


    return 0;
}


static void
runTests(string SEC_type){

    double cpu0_exec, cpu1_exec;
    cpu0_exec = get_wall_time();

    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    // graph 25 vertices epinions_S2: disconnected 
    // graph 25 vertices epinions_S3: disconnected (9-18-22)   
        

    // run all
    // int vert = 25; int type_min = 0, type_max = 3; int gclass_min = 1, glcass_max = 4; int ctype_min = 1, ctype_max = 7;

    // set instances
    int vert = 5;
    int type_min = 0, type_max = type_min+1;
    int gclass_min = 3, glcass_max = gclass_min+1;
    int ctype_min = 5, ctype_max = ctype_min+1;


    for(int i=type_min; i<type_max; i++)
        for(int gclass = gclass_min; gclass<glcass_max; gclass++) {// 1-3

            graphReader(0,vert, gclass, intances_types[i]);
            // graphReader(0,8, 1, "random");
            // graphReader(0,8, 2, "random");
            // graphReader(0,14, 1, "random");

            if(isConnected()){
                for(int ctype = ctype_min; ctype<ctype_max; ctype++){  // 1-6
                    IloEnv env;

                    try {

                        // instanceReader(0,vert, gclass,ctype, intances_types[i]);
                        // instanceReader(0,5, 3,5, "random");
                        instanceReader(0,14, 1,1, "random");

                        // create ILP problem
                        IloModel model(env);
                        cout << "[INFO]: Create variables" << endl;
                        BoolVar3Matrix x(env,num_vertices);
                        BoolVar3Matrix y(env,num_vertices);
                        BoolVar4Matrix f(env,num_vertices);
                        IntVarMatrix lambda(env,num_vertices);
                        allocVars(env, x, y, f, lambda);
                        
                        cout << "[INFO]: Create model" << endl;
                        if(SEC_type == "MTZ" || SEC_type == "mtz"){
                            cout << "[INFO]: MTZ constraint" << endl;
                            NumVar3Matrix r(env,num_vertices);
                            createModel_MTZ(model, x, y, f, lambda,r);
                        }
                        else if(SEC_type == "" || SEC_type == "NONE" || SEC_type == "none"){
                            cout << "[INFO]: No SECs constraints" << endl, createModel(model, x, y, f, lambda);}
                        else{
                            cout << "Not a valid SEC constraint" << endl, throw(-1);}
                        
                        // // createModel(model, x, y, f, lambda);
                        IloCplex cplex;
                        cplex = IloCplex(model);
                        // cplex.exportModel("TFP_SPC_model.lp");
                        
                        double cpu0, cpu1;
                        cplex.setParam(IloCplex::TiLim, 3600); // time limit 1h
                        cplex.setParam(IloCplex::TreLim, 7000); // memory limit 7GB
                        cplex.setParam(IloCplex::WorkMem, 10000);

                        // cplex.setOut(env.getNullStream()); // remove Logging output

                        cpu0 = get_wall_time();
                        if (!cplex.solve()){
                            env.error() << "[INFO]: Failed to optimize ILP." << endl;
                            cout << "Solution status = " << cplex.getStatus()   << endl;
                            throw(-1);
                        }
                        cpu1 = get_wall_time();
                        double runTime = cpu1 - cpu0;

                        cout << "Solution status = " << cplex.getStatus()   << endl;
                        cout << "Solution value  = " << cplex.getObjValue() << endl;
                        cout << "cplex time: " << cplex.getTime() << endl;
                        cout << "run time: " << runTime << endl;

                        // displaySolution(cplex,x,y,f,lambda);
                        // saveSolution(cplex,x,y,f,ctype);
                        // saveResults(cplex,runTime);

                    }catch (IloException& e) {
                        cerr << "Concert exception caught: " << e << endl;
                    }
                    catch (...) {
                        cerr << "Unknown exception caught" << endl;
                    }

                    env.end();
                }
            }else{cout << "[INFO]: Graph "<<  vert << "verticesS" << gclass << " ("<< intances_types[i] <<") is disconnected"   << endl;}
        }


    cout << "--------------------------------------------------------" << endl;
    cpu1_exec = get_wall_time();
    cout << "Total time: " << cpu1_exec - cpu0_exec << endl;

}

static void
allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda){

    // var x[u][j][s] : if worker u is in projecet j with skill s
    for (int u=0; u<num_vertices; u++){
        x[u] = BoolVarMatrix(env, num_teams);
        for (int j = 0; j < num_teams; j++)
            x[u][j] = IloBoolVarArray(env, num_skills);
    }

    // var y[u][v][j] : if worker u and v is in the same project j
    for (int u=0; u<num_vertices; u++){
        y[u] = BoolVarMatrix(env, num_vertices);
        for (int v = u+1; v < num_vertices; v++) // u<v
            y[u][v] = IloBoolVarArray(env, num_teams);
    }

    // var f[u][v][p][q]: if arc pq is used in the flow associate to path(u,v)
    for (int u=0; u<num_vertices; u++){
        f[u] = BoolVar3Matrix(env, num_vertices);
        for (int v=u+1; v<num_vertices; v++){ // u<v
            f[u][v] = BoolVarMatrix(env, num_vertices);
            for (int p=0; p<num_vertices; p++)
                f[u][v][p] = IloBoolVarArray(env, num_vertices);
        }
    }


    // var lambda[u][v]
    for (int u = 0; u < num_vertices; u++){
        lambda[u] = IloIntVarArray(env, num_vertices, 0, num_vertices*(num_vertices-1)); // must put the bounds
    }


}

static void
createModel (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda){

    // add var x[u][j][s] in model
    for(int j=0; j<num_teams; j++) // for all project j 
        for(int u=0; u< num_vertices; u++) // for all worker u
            for(int s=0; s<num_skills; s++){ 
                if(K[u][s] > 0 && R[j][s]>0){ // for all skills: s(u,j) skills of u and project j need skill u
                    char name[50];
                    sprintf(name, "x(%d,%d,%d)", u+1, j+1,s+1);
                    x[u][j][s].setName(name);
                    model.add(x[u][j][s]);
                }
            }

    // add var y[u][v][j]
    for(int u=0; u<num_vertices; u++) // for all u 
        for(int v = u+1; v<num_vertices; v++) // for all v: u<v
            for(int j = 0; j<num_teams; j++){ // for all project j
                char name[30];
                sprintf(name, "y%d%d%d", u+1, v+1,j+1);
                y[u][v][j].setName(name);
                model.add(y[u][v][j]);
            }


    // add var f[u][v][p][q]
    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++) // for all u,v: u<v
            if(G[u][v] == 0){ // and u,v not in E
                for(int p=0; p<num_vertices; p++)
                    for(int q=0; q<num_vertices; q++) // for all arc (p,q)
                        if(p!=q && G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                            char name[30];
                            sprintf(name, "F%d%d%d%d", u+1, v+1,p+1, q+1);
                            f[u][v][p][q].setName(name);
                            model.add(f[u][v][p][q]);
                        }
            }


    // add var lambda[u][v]
    for(int u = 0; u<num_vertices; u++) //
        for(int v = u+1; v<num_vertices; v++) // for all u,v : u<v
            if(G[u][v] == 0){ // and u,v not in E
                char name[20];
                sprintf(name, "lambd%d%d", u+1, v+1);
                lambda[u][v].setName(name);
                model.add(lambda[u][v]);
            }

    //object function
    objFunction(model,y,f);
    //constraints
    rest_OneTeam(model, x); // max one team with one skill
    rest_MinSkill(model, x); // min skill s per team j
    rest_LinY(model,x,y); // linearization y with x
    rest_NegEdge(model,y); // negative edge incompatibility
    rest_Flow(model,y,f); // flow const
    rest_flowSameTeam(model,y,f); // flow only if u,v is in the same team
    rest_PathComp(model,f,lambda); // every path with neg edges is pair
    //cuts
    // restCuts_2(model,f); // not repeat vertice in path

}

static void
createModel_MTZ (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y, BoolVar4Matrix f, IntVarMatrix lambda, NumVar3Matrix r){

    createModel(model, x, y, f, lambda);
//-------------------------MTZ VARIABLE---------------------------------------------
    IloEnv env = model.getEnv();
    
    // aloc var MTZ r[u][v][i]
    for (int u=0; u<num_vertices; u++){
        r[u] = NumVarMatrix(env, num_vertices);
        for (int v = u+1; v < num_vertices; v++) // u<v
            r[u][v] = IloNumVarArray(env, num_vertices,0,num_vertices, ILOFLOAT);
    }

    // add var MTZ r[u][v][i]: position of vertice i in u,v-flow
    for(int u=0; u<num_vertices; u++) // for all u 
        for(int v = u+1; v<num_vertices; v++) // for all v: u<v
            for(int i = 0; i<num_vertices; i++){ // for all vertice i
                char name[30];
                sprintf(name, "r%d%d%d", u+1, v+1,i+1);
                r[u][v][i].setName(name);
                model.add(r[u][v][i]);
            }


    //constraints
    restCycle_MTZ(model,f,r); // break cycles using MTZ constraint

    //valid ineq
    // restCuts_2(model,f); // not repeat vertice in path

}

static void
objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);
    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++) // for all u,v : u<v
            for(int p=0; p<num_vertices; p++)
                for(int q=0; q<num_vertices; q++)
                    if(G[u][v] == 0 && p!= q && G[p][q] != 0){ // (u,v) not in E+ and (p,q) is a possible arc between (u,v)-path
                       objExpr += f[u][v][p][q];
                    }

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j = 0; j<num_teams; j++){
                if(G[u][v] >= 1){objExpr += y[u][v][j];}  // (u,v) in E+
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();
}


static void
rest_OneTeam (IloModel model, BoolVar3Matrix x){

    IloEnv env = model.getEnv();

    // each worker uses at most one skill 
    for (int u=0; u<num_vertices; u++) {
        IloExpr expr(env);
        for (int j = 0; j <num_teams; j++) {
            for(int s = 0; s<num_skills; s++){
                if(R[j][s]>0 && K[u][s]>0) // team j need skill s and indidual u have skill s
                    expr += x[u][j][s];
            }
        }
        model.add(expr <= 1);   //var empty_sum to know if expr is empty or not, only add if not
        expr.end();
    }

}

static void
rest_MinSkill (IloModel model, BoolVar3Matrix x){
    IloEnv env = model.getEnv();

    // minimum number of individuals in project j with skill s
    for (int j=0; j <num_teams; j++) {
        for (int s=0; s <num_skills; s++) {
            IloExpr expr(env);
            if (R[j][s] > 0) { //for all skill s in project j
                for (int u=0; u < num_vertices; u++) {
                    if(K[u][s] > 0) // if worker u have skill s
                        expr += x[u][j][s];
                }
                model.add(expr >= R[j][s]);
                expr.end();
            }
        }
    }

}

static void
rest_LinY (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){
    IloEnv env = model.getEnv();

    // set var y with x (linearization)
    for (int u=0; u<num_vertices; u++)
        for (int v = u+1; v < num_vertices; v++) // for all u,v: u<v
            for (int j = 0; j <num_teams; j++){
                IloExpr exprU(env);
                IloExpr exprV(env);
                for (int s = 0; s < num_skills; s++) {
                    if(K[u][s] > 0 && R[j][s] > 0) // individual u have skill s and team j need skill s
                        exprU += x[u][j][s];
                    if(K[v][s] > 0 && R[j][s] > 0)
                        exprV += x[v][j][s];
                }
                model.add(exprU + exprV - y[u][v][j] <= 1);
                model.add(y[u][v][j]  <= exprU);
                model.add(y[u][v][j] <= exprV);
                exprU.end(); exprV.end();
            }

}

static void
rest_Flow (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v<num_vertices; v++) 
            if(G[u][v] == 0){                   // for all u<v and (u,v) not in E+
                for(int q = 0; q<num_vertices; q++){ // for fixed vertice q
                    IloExpr sumArcs(env);
                    for(int p = 0; p<num_vertices; p++) // for every arc (p,q) <=> (p,q) or (q,p) in E 
                        if(p!=q && G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                            sumArcs += f[u][v][p][q];
                            sumArcs += -f[u][v][q][p];
                        }
                    if(q == u){
                        IloExpr expr(env);
                        for (int j = 0; j < num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == -expr); expr.end();
                    }else if(q == v){
                        IloExpr expr(env);
                        for (int j = 0; j<num_teams; j++) {
                            expr += y[u][v][j];
                        }
                        model.add(sumArcs == expr); expr.end();
                    }else{
                        model.add(sumArcs == 0);
                    }
                }
            }

}


static void
rest_NegEdge (IloModel model, BoolVar3Matrix y){

    IloEnv env = model.getEnv();


    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v<num_vertices; v++)
            if(G[u][v] == -1){
                IloExpr expr(env);
                for(int j =0; j<num_teams; j++)
                    expr += y[u][v][j];
                model.add(expr == 0);
                expr.end();
            }


}


static void
rest_flowSameTeam (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    // only use flow between u,v if u,v work in the same team
    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v<num_vertices; v++)
            if(G[u][v] == 0){  // for all u,v in V: u<v and (u,v) not in E+
                for(int p=0; p<num_vertices; p++)
                    for(int q=0; q<num_vertices; q++)
                        if(p!=q && G[p][q] != 0){ // (p,q) in A
                            IloExpr expr(env);
                            for (int j=0; j <num_teams; j++) {
                                expr += y[u][v][j];
                            }
                            model.add(f[u][v][p][q] <= expr);
                            expr.end();
                        }
        }

}


static void
rest_PathComp (IloModel model, BoolVar4Matrix f, IntVarMatrix lambda){
    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v< num_vertices; v++){
            if(G[u][v] == 0 ){ // (u,v) not in E+
                IloExpr expr(env);
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!=q && G[p][q] < 0){ // sum of negatives arcs(A-)
                            expr += f[u][v][p][q];
                        }
                model.add(expr - 2*lambda[u][v] == 0); // pair number of negatives edges
                expr.end();
            }
        }
}

static void
restCycle_MTZ (IloModel model, BoolVar4Matrix f, NumVar3Matrix r){

    IloEnv env = model.getEnv();

    for(int u = 0; u<num_vertices;u++)
        for(int v = u+1; v<num_vertices;v++){
            for(int p=0;p<num_vertices;p++)
                for(int q=0;q<num_vertices;q++){
                    if(p!=q && abs(G[p][q])){
                        IloExpr expr(env);
                        expr = r[u][v][p] - r[u][v][q] + num_vertices*f[u][v][p][q];
                        model.add(expr <= (num_vertices-1));
                        expr.end();
                    }
                }
        }



}

static void
restCycle_DT(IloModel model, BoolVar4Matrix f){
    
    
    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v< num_vertices; v++)
            for(int p=0; p<num_vertices;p++)
                for(int q=0; q<num_vertices;q++){
                    if(p!=q && G[u][v] == 0){
                        model.add(f[u][v][p][q] + f[u][v][q][p] == 1);
                    }
                    for(int w = 0; w<num_vertices;w++)
                        if(p!=q && p!=w && w!=q)
                            model.add(f[u][v][p][q] + f[u][v][p][w] + f[u][v][w][q] <= 2);
                    
                }

}

static void
restCuts_1 (IloModel model, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    // pass only one time per vertice in a  u,v-flow
    // for(int u = 0; u<num_vertices; u++)
    //     for(int v = u+1; v< num_vertices; v++)
    //         for(int p=0; p<num_vertices;p++)
    //             for(int q=0; q<num_vertices;q++)
    //                 if(p!=q && G[p][q] = -1 && p == u && q==v){
    //                     model.add(f[u][v][p][q]==0);
    //                 }
        
}

static void
restCuts_2 (IloModel model, BoolVar4Matrix f){
    IloEnv env = model.getEnv();

    // pass only one time per vertice in a  u,v-flow
    for(int u = 0; u<num_vertices; u++)
        for(int v = u+1; v< num_vertices; v++){
            for(int q = 0; q<num_vertices; q++){
                // if(u!=q){
                    IloExpr expr(env);
                    for(int p=0; p<num_vertices;p++)
                        expr += f[u][v][p][q]; 
                
                   model.add(expr<=1);
                   expr.end(); 
                // }

            }
        }
}

static void
displaySolution(IloCplex& cplex,
               BoolVar3Matrix x,
               BoolVar3Matrix y,
               BoolVar4Matrix f,
               IntVarMatrix lambda){



    // cout << "Solution status = " << cplex.getStatus()   << endl;
    // cout << "Solution value  = " << cplex.getObjValue() << endl;
    // cout << "Valores das Variaveis " << endl;

    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;

    for(int j=0; j<num_teams; j++)
        for(int u=0; u<num_vertices; u++)
            for(int s=0; s<num_skills; s++){
                if(K[u][s] > 0 && R[j][s]>0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        cout << "\t " << "x"<< "["<< u+1  << "," << j+1 << "," << s+1 << "]" << " -> " << " 1 \n";
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j=0; j<num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    cout << "\t " << "y"<< "["<< u+1  << "," << v+1 << "," << j+1 << "]" << " -> " << " 1 \n";
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<num_vertices; u++)
        for(int v = 0; v<num_vertices; v++)
            if(u<v && G[u][v] == 0){
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!= q && G[p][q] != 0){
                            if(cplex.getValue(f[u][v][p][q]) > 0){
                                cout << "\t " << "f"<< "["<< u+1  << "," << v+1 << "," << p+1 << ", " << q+1 << "]" << " -> " << " 1 \n";
                            }
                        }
            }

}


static void
    saveSolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y,
                   BoolVar4Matrix f,
                   int class_type){


    char arq[1000];
    char arqv_instance[50];
    sprintf(arqv_instance, "%s_%d", instanceG,class_type);
    // sprintf(arq, "%s/results/%s.txt",CURRENT_DIR, arqv_instance);
    sprintf(arq, "%s/results/%d_Vertices_%d-%d-%d_TSP_SPC.ods",CURRENT_DIR, num_vertices, current_year, current_month, current_day);

    FILE *file = fopen(arq, "w");
    if (file == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(file, "Solution value = %.1f \n", cplex.getObjValue());

    for(int j=0; j<num_teams; j++)
        for(int u=0; u<num_vertices; u++)
            for(int s=0; s<num_skills; s++){
                if(K[u][s] > 0 && R[j][s]>0){
                    if(cplex.getValue(x[u][j][s]) > 0){
                        fprintf(file, "x[%d,%d,%d] = 1 \n", u+1, j+1, s+1);
                    }
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j=0; j<num_teams; j++){
                if(cplex.getValue(y[u][v][j]) > 0){
                    fprintf(file, "y[%d,%d,%d] = 1 \n", u+1, v+1, j+1);
                }
            }

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<num_vertices; u++)
        for(int v = 0; v<num_vertices; v++)
            if(u<v && G[u][v] == 0){
                for(int p = 0; p<num_vertices; p++)
                    for(int q = 0; q<num_vertices; q++)
                        if(p!= q && G[p][q] != 0){
                            if(cplex.getValue(f[u][v][p][q]) > 0){
                                fprintf(file, "f[%d,%d,%d,%d] = 1 \n", u+1, v+1, p+1,q+1);
                            }
                        }
            }


    fclose(file);

}


static void
    saveResults(IloCplex& cplex,
                   double timeF){

    char arq[600];
    sprintf(arq, "%s/results/%d_Vertices_2022-06-27_cycle_mtz.ods",CURRENT_DIR, num_vertices);
    
    ofstream outputTable;
    outputTable.open(arq,ios:: app);
    if(outputTable.is_open()){

        outputTable << instanceG << ";"; // grafo instancia
        outputTable << instanceKR << ";"; // KR instancia
        outputTable << num_vertices << ";";   // numero de vertices
        outputTable << num_skills << ";";   // qtd de hab
        outputTable << num_teams << ";";   // qtd de proj
        outputTable << cplex.getStatus() << ";"; // Status cplex
        outputTable << cplex.getObjValue() << ";"; // valor fo
        outputTable << cplex.getNnodes() << ";"; // num nos
        outputTable << cplex.getMIPRelativeGap() <<";"; // gap
        outputTable << timeF <<  ";"; // tempo execucao flp
        outputTable << cplex.getTime() <<  ";"; // tempo execucao cplex
        outputTable << " \n ";


    }

    outputTable.close();

}


static void
    testeGrafos(){

    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    
    cout << "-------------------------------------------------" << endl;
    int vert = 50;
    for(int i=0; i<3; i++){
        for(int gclass = 1; gclass<4; gclass++){  // 1-3
            
            instanceReader(0,vert, gclass,1, intances_types[i]);
            
            if(isConnected())
                cout << "The Graph is connected." << endl;
            else
                cout << "The Graph is not connected." << endl;


            cout << "-------------------------------------------" << endl;
            
        }
    }
    
}

