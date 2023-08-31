//#include <iostream>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <time.h>
#include <sys/time.h>
#include <string>
#include <filesystem>
#include <list>
#include<queue>
#include <cmath>        // std::abs

#include<fstream>
using namespace std;

// ********** VARIAVEIS GLOBAIS **********
int num_vertices;/* numero de vertices/individuos do grafo*/
int num_skills;/*numero de habilidades*/
int num_teams; /* quantidade de total de equipes/projetos*/

int **G; /*matriz |n|x|n| grafo*/
int **G_Zuv; /*matriz |n|x|n| grafo ponderado pelos caminhos minimos positivos*/
int **K; /*matriz |n|x|f| de pessoa - habilidade*/
int **R; /* matriz demanda |m|x|f|: quantidade de cada pessoa com habilidade k_a em cada equipe*/

char instanceG[50];
char instanceKR[50];

char CURRENT_DIR[500];

int current_day, current_month, current_year;

int total=25;

//double runTime;
#define BILLION 1000000000L

#define INF 1000000

//typedef
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;
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
void print_matrix(int **matrix){
    
    cout << " Adjacency Matrix \n";
    for (int u = 0; u < num_vertices; u++){
        for (int v = 0; v < num_vertices; v++)
            cout <<  matrix[u][v] << " ";
        cout<< "\n";
    }
}

//load matrix
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
        print_matrix(G);
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
instanceReader(bool show, int vert, int G_class, int class_type, string instance_type){
    char arq1[600]; char arq2[600]; char arq3[600];
    char arq_base[500];
    sprintf(arq_base, "%s/instances_COR", CURRENT_DIR);

    if (instance_type == "random"){
        sprintf(instanceG, "%dverticesS%d", vert, G_class);
    }
    else if (instance_type == "bitcoinotc" || instance_type == "epinions"){
        sprintf(instanceG, "%dvertices_%s_S%d", vert, instance_type.c_str(), G_class);

    }else{
        cout << "[INFO] not a valid instance type" << endl;
    }


    sprintf(arq1, "%s/%dVertices/%s.txt", arq_base, vert, instanceG);
    cout << "[INFO] Read Graph: "<< instanceG << endl;

    sprintf(instanceKR, "%dVertices/class1/%d", vert, class_type);
    sprintf(arq2, "%s/%s/K.txt", arq_base, instanceKR);
    sprintf(arq3, "%s/%s/R.txt", arq_base, instanceKR);
    cout << "[INFO] Instance class type: "<< instanceKR << endl;


    load_Graph(arq1, show, instance_type); 
    load_K(arq2, show); load_R(arq3, show);

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
    instanceReader (bool show, int vert, int G_class, int class_type, string instance_type);

static void
    allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y);

static void
    allocVars_uv (IloEnv env, BoolVarMatrix f, IloIntVar lambda, IloNumVarArray r);

static void
    createModel_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, IloNumVarArray r, int u, int v);

static void
   createModel_Simplifed (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   objFunction_uv (IloModel model, BoolVarMatrix f, int u, int v);

static void
   objFunction_Simplifed (IloModel model, BoolVar3Matrix y);

static void
   rest_OneTeam (IloModel model, BoolVar3Matrix x);

static void
   rest_MinSkill (IloModel model, BoolVar3Matrix x);

static void
   rest_LinY (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

static void
   rest_Flow_uv (IloModel model, BoolVarMatrix f, int u, int v);

static void
   rest_PathComp_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v);

static void    
    rest_NegEdge (IloModel model, BoolVar3Matrix y);

static void
   restCycleMTZ_uv (IloModel model, BoolVarMatrix f, IloNumVarArray r, int u, int v);

static void
   restCut_uv_1 (IloModel model, BoolVarMatrix f, int u, int v);

static void
    displaySolution(IloCplex& cplex,
                   BoolVar3Matrix x,
                   BoolVar3Matrix y);

static void
    displaySolution_uv(IloCplex& cplex,
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
    createG_Zuv();

static void
    initG_Zuv();

static void
    createG_Zuv_algorithm();

static void
    testeGrafos();
static void 
    traverse(int u, bool visited[]);
static bool 
    isConnected(); 

static void
    compare_Zmatrices();


//--------------------------------------------------------------- DIJKSTRA COMPATIVEL ---------------------------------------------------
bool all_visited(bool vSet[]){

    for(int i =0; i<num_vertices;i++){
        if(vSet[i] == false)
            return false;
    }

    return true;
}
bool isCycle(int ***path, int u, int v, int j){
    
    for (int k = 0; k < num_vertices; k++)
        if(path[u][j][0] != -1)
            if(path[u][j][k] == v){
                return true;
            }
    return false;
}

void add_VerticePath(int ***path_u, int u,int ***path_v, int v){
    
    int index;
    for (int j = 0; j < total; j++)
        if(!isCycle(path_u,u,v,j))
            for (int k = 0; k<num_vertices;k++)
                path_v[v][j][k] = path_u[u][j][k];

    // for (int j = 0; j < num_vertices; j++){
    //     if(!isCycle(path_u,u,v,j)){
    //         for (int k = 0; k<num_vertices;k++)
    //             path_v[v][j][k] = path_u[u][j][k];}
    //     else{
    //         for (int k = 0; k<num_vertices;k++)
    //             path_v[v][j][k] = -1;
    //     }
    // }

    for (int j = 0; j < total; j++){
        if(path_v[v][j][0] != -1 && !isCycle(path_u,u,v,j)){
            for (int k = num_vertices-1; k>=0;k--)
                if(path_v[v][j][k] == -1)
                    index = k;
            path_v[v][j][index] = v;      
        }
    }

}

int empty_index(int ***path,int u){

    int index = -1;
    for (int j = total-1; j >= 0; j--){
        if(path[u][j][0] == -1)
            index = j;
    }
    return index;
}

void zero_Paths(int ***path, int v){
    // int total = 1000;
    for (int j = 0; j < total; j++)
        for (int k = 0; k<num_vertices;k++)
            path[v][j][k] = -1;

}

void add_newPath(int ***path_u, int u,int ***path_v, int v){
    // int total = 1000;
    int index;
    int last;
    for (int j = 0; j < total; j++){        
        if(path_u[u][j][0] != -1 && !isCycle(path_u,u,v,j) && index!=-1){
            index = empty_index(path_v,v);
            // cout << "aqui" << endl;
            // cout << j << endl;
            for (int k = 0; k<num_vertices;k++){
                path_v[v][index][k] = path_u[u][j][k];
            }
            for (int k = num_vertices-1; k>=0;k--)
                if(path_v[v][index][k] == -1)
                    last = k;
            path_v[v][index][last] = v;  
        }
        
    }


}

bool equal_Vectors(int *vector1, int *vector2){

    for(int i=0;i<num_vertices;i++)
        if(vector1[i] != vector2[i])
            return false;

    return true;


}


bool isVector_inMatrix(int *vector, int **matrix){

    // int total = 1000;

    int *vector2;
    vector2 = new int[num_vertices];

    for(int i=0; i<num_vertices;i++){
        for(int j=0;j<total;j++)
            vector2[j] = matrix[i][j];
        if(equal_Vectors(vector,vector2))
            return true;
    }

    return false;
}


bool is_newPath(int ***path_u, int u,int ***path_v, int v){

    // int total = 1000;

    int *vector_path;
    vector_path = new int[num_vertices];
 
    int **matrix_paths;
    matrix_paths = new int*[total];
    for(int j = 0; j <total; j++)
        matrix_paths[j] = new int[num_vertices];

    for(int j = 0; j < total; j++)
        for(int k = 0; k < num_vertices; k++){
            if(path_v[v][j][k] != v)
                matrix_paths[j][k] = path_v[v][j][k];
            if(path_v[v][j][k] == v)
                matrix_paths[j][k] = -1;
        }

    for (int j_u = 0; j_u < total; j_u++){
        if(path_u[u][j_u][0] != -1){
            for(int k = 0; k < num_vertices; k++)
                vector_path[k] = path_u[u][j_u][k];
            if(!isVector_inMatrix(vector_path,matrix_paths)){
                // for (int i = 0; i < num_vertices; ++i)
                //     delete [] matrix_paths[i];
                // delete [] matrix_paths;
             
                return true;
            }  
        }
    }


    // for (int i = 0; i < num_vertices; ++i)
    //     delete [] matrix_paths[i];
    // delete [] matrix_paths;

    return false;

}

bool isVertice_Cycle(int ***path_u, int u, int v){

    // int total = 1000;

    for (int j = 0; j < total; j++){
        if(path_u[u][j][0] != -1 && !isCycle(path_u,u,v,j)){
            return false;
        }
        
    }

    return true;

}
bool isPaths_void(int ***path_v, int v){

    for(int j =0;j<total;j++)
        for(int k =0; k < num_vertices;k++)
            if (path_v[v][j][k] != -1)
                return false;
    return true;

}

void print_paths(int ***path, int u){
    
    // int total = 1000;

    int cont = 0;
    cout << "--------------" << endl;
    cout << "caminho "<< u+1 << " :" << endl;
    for(int j = 0; j<total;j++){
        if(path[u][j][0]!= -1){
            cont+=1;
            for(int k = 0; k<num_vertices;k++){
                cout << path[u][j][k]+1 << ", ";}
        cout << endl;}
    }
    cout << "Qte de caminhos: " << cont << endl; 
    cout << "--------------" << endl;

}

void SCP_algorithm(int src)
{

    // int total = 1000;
    int ***PosPath; int ***NegPath;


    PosPath = (int***)(malloc(num_vertices*sizeof(int**)));
    NegPath = (int***)(malloc(num_vertices*sizeof(int**)));

    for (int i=0;i<num_vertices;i++){
        PosPath[i] = (int**)malloc(total*sizeof(int*)), NegPath[i] = (int**)malloc(total*sizeof(int*));
        for(int j = 0; j <total; j++){
            PosPath[i][j] = (int*)malloc(num_vertices*sizeof(int)), NegPath[i][j] = (int*)malloc(num_vertices*sizeof(int));
        }
    }


    // int ***PosPath; int ***NegPath;
    // PosPath = new int **[num_vertices]; NegPath = new int **[num_vertices];
    // for(int i = 0; i <num_vertices; i++){
    //     PosPath[i] = new int *[total], NegPath[i] = new int *[total];
    //     for(int j = 0; j <total; j++){
    //         PosPath[i][j] = new int[num_vertices], NegPath[i][j] = new int[num_vertices];
    //     }
    // }

    int dist[num_vertices][2];  // [|neg_path|, |pos_path|]

    bool vSet[num_vertices]; // vSet[i] list vertex to be processed

 
    // Initialize all distances as INFINITE and vSet[] as false
    for (int i = 0; i < num_vertices; i++)
        vSet[i] = false, dist[i][0] = INF, dist[i][1] = INF; //dist[i][0] = INT_MAX, dist[i][1] = INT_MAX;
    
    for (int i = 0; i < num_vertices; i++)
        for (int j = 0; j < total; j++)
            for (int k = 0; k < num_vertices; k++)
                PosPath[i][j][k] = -1, NegPath[i][j][k] = -1;

    // Distance of source vertex from itself is always 0
    dist[src][0] = 0; dist[src][1] = 0; 
    // for (int k = 0; k < num_vertices; k++)
    PosPath[src][0][0] = src, NegPath[src][0][0] = src;


    // NegPath[0][1][0] = 6;
    // NegPath[0][1][1] = 5; 
    // NegPath[0][1][2] = 1;
    // NegPath[0][1][3] = 0;

    // NegPath[2][0][0] = 6;
    // NegPath[2][0][1] = 5; 
    // NegPath[2][0][3] = 3;
    // NegPath[2][0][4] = 2;

    // NegPath[2][1][0] = 6;
    // NegPath[2][1][1] = 5; 
    // NegPath[2][1][2] = 1;
    // NegPath[2][1][4] = 2;


    // int v = 1;
    // if(isVertice_Cycle(NegPath,2,v))
    //     cout << "ciclo" << endl;
    // if(!isVertice_Cycle(NegPath,2,v))
    //     cout << "sem ciclo" << endl;



    // cout << "--------------" << endl;
    // cout << "caminho negativo 3:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(NegPath[2][j][0]!= -1){
    //         for(int k = 0; k<num_vertices;k++)
    //             cout << NegPath[2][j][k]+1 << ", ";
    //     cout << endl;}
    // }
    // cout << "custo = " << dist[2][0] << endl;
    // cout << "--------------" << endl;


    // NegPath[0][0][0] = 6;
    // NegPath[0][0][1] = 5; 
    // NegPath[0][0][2] = 1;
    // NegPath[0][0][3] = 0;
    

    // add_newPath(NegPath, 0,NegPath, 2);
    // // // NegPath[2][1][0] = 6;
    // // // NegPath[2][1][1] = 5; 
    // // // NegPath[2][1][2] = 1;
    // // // NegPath[2][1][3] = 0;
    // // // NegPath[2][1][4] = 2;


    // cout << "--------------" << endl;
    // cout << "caminho negativo 3:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(NegPath[2][j][0]!= -1){
    //         for(int k = 0; k<num_vertices;k++)
    //             cout << NegPath[2][j][k]+1 << ", ";
    //     cout << endl;}
    // }
    // cout << "custo = " << dist[2][0] << endl;
    // cout << "--------------" << endl;


    // // add_VerticePath(NegPath, 2,PosPath, 4);
    
    // add_newPath(NegPath, 2,PosPath, 4);

    // cout << "--------------" << endl;
    // cout << "caminho positivos 5:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(PosPath[4][j][0]!= -1){
    //         for(int k = 0; k<num_vertices;k++)
    //             cout << PosPath[4][j][k]+1 << ", ";
    //     cout << endl;}
    // }
    // cout << "custo = " << dist[4][0] << endl;
    // cout << "--------------" << endl;




    // if(!isVertice_Cycle(PosPath,4,3))
    //     cout << "sem ciclo" << endl;
    // else
    //     cout << " ciclo" << endl;
    
    
    // // add_VerticePath(PosPath, 4,NegPath,3);
    // add_newPath(PosPath, 4,NegPath,3);

    // cout << "--------------" << endl;
    // cout << "caminho negativo 4:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(NegPath[3][j][0]!= -1){
    //         for(int k = 0; k<num_vertices;k++)
    //         cout << NegPath[3][j][k]+1 << ", ";
    //     cout << endl;}
    // }
    // cout << "custo = " << dist[3][0] << endl;
    // cout << "--------------" << endl;


    // exit(0);

    // while have a vertex to be processed 
    while(!all_visited(vSet)){

        int min = INF, min_index;
        // for (int v = 0; v < num_vertices; v++){
        //     if (vSet[v] == false && dist[v][1] < min)
        //         min = dist[v][1], min_index = v;
        //     if (vSet[v] == false && dist[v][0] < min)
        //         min = dist[v][0], min_index = v;
        // }
        // for (int v = num_vertices-1; v >= 0; v--){
        //     if (vSet[v] == false && dist[v][0] < min)
        //         min = dist[v][0], min_index = v;
        //     if (vSet[v] == false && dist[v][1] < min)
        //         min = dist[v][1], min_index = v;
        // }
        
        for (int v = 0; v < num_vertices; v++){
            if (vSet[v] == false && dist[v][1] < min)
                min = dist[v][1], min_index = v;}
         for (int v = 0; v < num_vertices; v++){
            if (vSet[v] == false && dist[v][0] < min)
                min = dist[v][0], min_index = v;}



        int u = min_index;

        vSet[u] = true;
        // cout << "-----------------------------------" << endl;
        // cout << "vertice u=" << u+1 << endl;
        // cout << "custos: " << dist[u][0] << ", " << dist[u][1] << endl;
        for (int v = 0; v < num_vertices; v++){
                

                // --------------------------------------------Aresta Positiva----------------------------------------------------------------
                if(G[u][v] > 0 && u!=v){ 
                    // cout << "--------------" << endl;
                    // cout << "aresta positiva " << u+1 << ", " << v+1 << endl;
                    // print_paths(NegPath,u);
                    // print_paths(PosPath,u);
                    // print_paths(NegPath,v);
                    // print_paths(PosPath,v);
                    if(dist[u][1]!= INF && abs(dist[u][1]) + abs(G[u][v]) <= abs(dist[v][1]) && !isVertice_Cycle(PosPath,u,v)){
                         
                        if(abs(dist[u][1]) + abs(G[u][v]) == abs(dist[v][1]) && is_newPath(PosPath,u,PosPath,v)){
                            // if(v == 13){
                                // cout << "valores iguais "<< endl;
                            // print_paths(PosPath,u);
                            // print_paths(PosPath,v);
                            add_newPath(PosPath,u,PosPath,v);
                        }
                        // add_VerticePath(PosPath,u,PosPath,v);
                        // add_newPath(NegPath,u,NegPath,v);
                        if(abs(dist[u][1]) + abs(G[u][v]) < abs(dist[v][1])){
                            zero_Paths(PosPath,v);
                            add_newPath(PosPath,u,PosPath,v);
                            dist[v][1] = dist[u][1] + abs(G[u][v]), vSet[v] = false;
                            // if(v == 13)
                                // cout << "reset 13 "<< endl;
                            
                        }
                        // cout << "custos v: " << dist[v][0] << ", " << dist[v][1] << endl;

                        
                    }else if (dist[v][1] == INF && !isVertice_Cycle(PosPath,u,v) && isPaths_void(PosPath,v)){
                        // add_VerticePath(PosPath,u,PosPath,v);
                        add_newPath(PosPath,u,PosPath,v);
                    }
                    
                    if(u!= src && dist[u][0] != INF && abs(dist[u][0]) + abs(G[u][v]) <= abs(dist[v][0]) && !isVertice_Cycle(NegPath,u,v)){
                        
                        if(abs(dist[u][0]) + abs(G[u][v]) == abs(dist[v][0])){// && is_newPath(NegPath,u,NegPath,v)){
                            
                            // add_newPath(NegPath,u,NegPath,v);
                            // cout << "valores iguais "<< endl;
                        }
                        // add_VerticePath(NegPath,u,NegPath,v);
                        // add_newPath(NegPath,u,NegPath,v);
                        if(abs(dist[u][0]) + abs(G[u][v]) < abs(dist[v][0])){
                            zero_Paths(NegPath,v);
                            add_newPath(NegPath,u,NegPath,v);
                            dist[v][0] = dist[u][0] + abs(G[u][v]), vSet[v] = false; 
                        } 
                        // cout << "custos v: " << dist[v][0] << ", " << dist[v][1] << endl;

                    }
                    else if (dist[v][0] == INF && !isVertice_Cycle(NegPath,u,v) && isPaths_void(NegPath,v)){
                        // add_VerticePath(NegPath,u,NegPath,v);
                        add_newPath(NegPath,u,NegPath,v);
                    }

                }
                // --------------------------------------------Aresta Negativa----------------------------------------------------------------
                if(G[u][v] < 0 && u!=v){

                    // cout << "--------------" << endl;
                    // cout << "aresta negativa " << u+1 << ", " << v+1 << endl;
                    // print_paths(NegPath,u);
                    // print_paths(PosPath,u);
                    // print_paths(NegPath,v);
                    // print_paths(PosPath,v);

                    if(dist[u][1] != INF && abs(dist[u][1]) + abs(G[u][v]) <= abs(dist[v][0]) && !isVertice_Cycle(PosPath,u,v)){
                         

                        if(abs(dist[u][1]) + abs(G[u][v]) == abs(dist[v][0])){//  && is_newPath(PosPath,u,NegPath,v)){
                         
                            // add_newPath(PosPath,u,NegPath,v);
                            // cout << "valores iguais "<< endl;
                        }
                        // add_VerticePath(PosPath,u,NegPath,v);
                        // add_newPath(PosPath,u,NegPath,v);
                        if(abs(dist[u][1]) + abs(G[u][v]) < abs(dist[v][0])){   
                            zero_Paths(NegPath,v);         
                            add_newPath(PosPath,u,NegPath,v);
                            dist[v][0] = dist[u][1] + abs(G[u][v]), vSet[v] = false;
                        }
                        // cout << "custos v: " << dist[v][0] << ", " << dist[v][1] << endl;

                        
                    }else if (dist[v][0] == INF && !isVertice_Cycle(PosPath,u,v) && isPaths_void(NegPath,v)){
                        // add_VerticePath(PosPath,u,NegPath,v);
                        add_newPath(PosPath,u,NegPath,v);
                    }
                    
                    if(u!= src && dist[u][0] != INF && abs(dist[u][0]) + abs(G[u][v]) <= abs(dist[v][1]) && !isVertice_Cycle(NegPath,u,v)){
                        
                        if(abs(dist[u][0]) + abs(G[u][v]) == abs(dist[v][1])){//  && is_newPath(NegPath,u,PosPath,v)){
                            // add_newPath(NegPath,u,PosPath,v);
                            // cout << "valores iguais "<< endl;
                        }

                        // add_VerticePath(NegPath,u,PosPath,v);
                        // add_newPath(NegPath,u,PosPath,v);
                        if(abs(dist[u][0]) + abs(G[u][v]) < abs(dist[v][1])){
                            zero_Paths(PosPath,v);
                            add_newPath(NegPath,u,PosPath,v);
                            dist[v][1] = dist[u][0] + abs(G[u][v]), vSet[v] = false;}  
                        // cout << "custos v: " << dist[v][0] << ", " << dist[v][1] << endl;



                    }else if (dist[v][1] == INF && !isVertice_Cycle(NegPath,u,v) && isPaths_void(PosPath,v)){
                        // add_VerticePath(NegPath,u,PosPath,v);
                        add_newPath(NegPath,u,PosPath,v);
                    }
            
                }
        }
        
    }
    // cout << "--------------" << endl;
    // cout << "caminho negativo 4:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(NegPath[3][j][0]!= -1)
    //         for(int k = 0; k<num_vertices;k++)
    //         cout << NegPath[3][j][k]+1 << ", ";
    //     cout << endl;
    // }
    // cout << "custo = " << dist[3][0] << endl;
    // cout << "--------------" << endl;

    // cout << "--------------" << endl;
    // cout << "caminho negativo 3:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(NegPath[2][j][0]!= -1)
    //         for(int k = 0; k<num_vertices;k++)
    //         cout << NegPath[2][j][k]+1 << ", ";
    //     cout << endl;
    // }
    // cout << "custo = " << dist[2][0] << endl;
    // cout << "--------------" << endl;

    // cout << "--------------" << endl;
    // cout << "caminho positivo 5:" << endl;
    // for(int j = 0; j<num_vertices;j++){
    //     if(PosPath[4][j][0]!= -1)
    //         for(int k = 0; k<num_vertices;k++)
    //         cout << PosPath[4][j][k]+1 << ", ";
    //     cout << endl;
    // }
    // cout << "custo = " << dist[4][1] << endl;
    // cout << "--------------" << endl;


    for (int v=0;v<num_vertices;v++){
        // G_Zuv[src][v] = dist[v][1];
        // G_Zuv[v][src] = dist[v][1];
        if (G_Zuv[v][src] != 0){G_Zuv[src][v] = G_Zuv[v][src];} //min(G_Zuv[v][src],dist[v][1]);}
        else {G_Zuv[v][src] = dist[v][1];}
        if (G_Zuv[src][v] != 0){G_Zuv[src][v] =  G_Zuv[src][v];} // min(G_Zuv[src][v],dist[v][1]);}
        else {G_Zuv[src][v] = dist[v][1];}
        // G_Zuv[src][v] = min(G_Zuv[src][v],dist[v][1]);
    }

    /* DESALOCAR MEMORIA */

    // for(int i=0;i<num_vertices;i++)
    // {
    //     for(int j=0;j<total;j++)
    //     {
    //         delete [] NegPath[i][j];
    //         delete [] PosPath[i][j];
    //     }
    //     delete [] NegPath[i];
    //     delete [] PosPath[i];
    //     delete [] dist[i];
    // }
    // delete [] NegPath;
    // delete [] PosPath;
    // delete [] dist;
    // delete [] vSet;


}
//---------------------------------------------------------------------------------------------------------------------------------------------------
static void
    runTests(string Z_mode);

int main()
{
    getcwd(CURRENT_DIR, 500);

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    current_day = tm.tm_mday; current_month = tm.tm_mon + 1; current_year = tm.tm_year + 1900;


    // ps. SCP_algorithm have a counterexample, not work in all cases //

    runTests("MTZ");
    // runTests("SCP");


    return 0;
}


static void
runTests(string Z_mode){


    double cpu0_exec, cpu1_exec;
    cpu0_exec = get_wall_time();

    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    // G 25 vertices epinions_S3: disconnected (9-18-22)


// run all
    // int vert = 25; int type_min = 0, type_max = 3; int gclass_min = 1, glcass_max = 4; int ctype_min = 1, ctype_max = 7;

// set instances
    int vert = 14;
    int type_min = 0, type_max = type_min+1;
    int gclass_min = 1, glcass_max = gclass_min+1;
    int ctype_min = 1, ctype_max = ctype_min+1;

    for(int i=type_min; i<type_max; i++)
        for(int gclass = gclass_min; gclass<glcass_max; gclass++) {// 1-3

            graphReader(0,vert, gclass, intances_types[i]);

            if(isConnected()){
                for(int ctype = ctype_min; ctype<ctype_max; ctype++){  // 1-6

                    IloElambda{

                        instanceReader(0,vert, gclass,ctype, intances_types[i]);
                        // instanceReader(0,5, 3,5, "random");
                        // instanceReader(0,6, 2,3, "random");
                        // instanceReader(0,8, 1,1, "random");
                        // instanceReader(0,14, 1,1, "random");


                        cout << "[INFO]: Create Graph Zuv" << endl;
                        double timeZ_start = get_wall_time();
                        
                        // Subproblems
                        if(Z_mode == "mtz" || Z_mode == "MTZ"){
                            createG_Zuv();
                        }
                        // Algorithm 
                        if(Z_mode == "scp" || Z_mode == "SCP"){
                            initG_Zuv();
                            for(int u=0; u<num_vertices; u++)
                                SCP_algorithm(u);
                        }
                            
                        double timeZ_final = get_wall_time();
                        double timeG = timeZ_final - timeZ_start;

                        // create ILP problem
                        IloModel model(env);
                        cout << "[INFO]: Create variables simplified" << endl;
                        BoolVar3Matrix x(env,num_vertices);
                        BoolVar3Matrix y(env,num_vertices);
                        allocVars(env, x, y);
                        cout << "[INFO]: Create model simplified" << endl;
                        createModel_Simplifed(model, x, y);
                        IloCplex cplex;
                        cplex = IloCplex(model);
                        // cplex.exportModel("TFP_SPC_Simpl.lp");
                        
                        double cpu0, cpu1;
                        cplex.setParam(IloCplex::TiLim, 3600); // time limit 2h (7200)
                        cplex.setParam(IloCplex::TreLim, 30000); // memory limit 7GB
                        // cplex.setParam(IloCplex::WorkMem, 2000);

                        // cplex.setOut(env.getNullStream()); // remove Logging output


                        cpu0 = get_wall_time();
                        if (!cplex.solve()){
                            env.error() << "[INFO]: Failed to optimize ILP." << endl;
                            cout << "Solution status = " << cplex.getStatus()   << endl;
                            throw(-1);
                        }
                        cpu1 = get_wall_time();
                        double runTime = (cpu1 - cpu0) + timeG;
                        

                        cout << "Solution status = " << cplex.getStatus()   << endl;
                        cout << "Solution value  = " << cplex.getObjValue() << endl;
                        cout << "cplex time: " << cplex.getTime() << endl;
                        cout << "run time: " << runTime << endl;

                        displaySolution(cplex,x,y);
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
initG_Zuv(){

    G_Zuv = (int**)(malloc(num_vertices*sizeof(int*)));
    for (int u=0;u<num_vertices;u++)
        G_Zuv[u] = (int*)malloc(num_vertices*sizeof(int));


    for(int u=0; u<num_vertices;u++)
        for(int v=0; v<num_vertices; v++)
            if(u==v || G[u][v] == 1)
                G_Zuv[u][v] = 1;
            else
                G_Zuv[u][v] = 0;


}


static void
allocVars_uv (IloEnv env, BoolVarMatrix f, IloIntVar lambda, IloNumVarArray r){

    // var f(u,v)[p][q]: if arc pq is used in the flow associate to path(u,v)
    for (int p=0; p<num_vertices; p++)
        f[p] = IloBoolVarArray(env, num_vertices);

}

static void
allocVars (IloEnv env, BoolVar3Matrix x, BoolVar3Matrix y){

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


}

static void
createG_Zuv(){

    initG_Zuv();

    for(int u=0; u<num_vertices;u++)
        for(int v=u+1; v<num_vertices; v++){
            if(G[u][v] == 0){

                IloEnv env;
                try{
                    IloModel model(env);
                    
                    BoolVarMatrix f(env,num_vertices);
                    IloIntVar lambda(env);
                    IloNumVarArray r(env, num_vertices,0,num_vertices, ILOFLOAT);
                    allocVars_uv(env,f,lambda,r);
                    createModel_uv(model,f,lambda,r,u,v);       // can optmize var f, put all dif u,v = zero or better use 2 indices, only pq
                    IloCplex cplex;
                    cplex = IloCplex(model);
                    cplex.setOut(env.getNullStream());
                    
                    if (!cplex.solve()){
                        
                        // infeasible -> infinite
                        G_Zuv[u][v] = INF; 
                        G_Zuv[v][u] = INF;
                        throw(-1);
                    }


                    G_Zuv[u][v] = int(cplex.getObjValue());
                    G_Zuv[v][u] = int(cplex.getObjValue());

                    // cout << "-------------------------------------------------- \n\n" << endl;
                }catch (IloException& e) {
                        cerr << "Concert exception caught: " << e << endl;
                    }
                    catch (...) {
                        // cerr << "Unknown exception caught" << endl;
                    }

                env.end();


            }if(G[u][v] == 1){
                G_Zuv[u][v] = 1;
                G_Zuv[v][u] = 1;
            }if(G[u][v] == -1){
                // incompatible -> infinite
                G_Zuv[u][v] = INF;
                G_Zuv[v][u] = INF;
            }

        }
}

static void
createModel_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, IloNumVarArray r, int u, int v){

    // add var f(u,v)[p][q]
    if(G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<num_vertices; p++)
            for(int q=0; q<num_vertices; q++) // for all arc (p,q)
                if(p!=q && G[p][q] != 0){ // must have an edge(p,q) to have an arc (p,q)
                    char name[50];
                    sprintf(name, "F(%d,%d)%d%d", u+1, v+1,p+1, q+1);
                    f[p][q].setName(name);
                    model.add(f[p][q]);
                }
    }

    // add var lambda(u,v)
    if(G[u][v] < 1){ // and u,v not in E{+}
        char name[20];
        sprintf(name, "lambd(%d,%d)", u+1, v+1);
        lambda.setName(name);
        model.add(lambda);
    }


    // add var r(u,v)[p]
    if(G[u][v] < 1){ // and u,v not in E{+}
        for(int p=0; p<num_vertices; p++){
            char name[50];
            sprintf(name, "r(%d,%d)%d", u+1, v+1,p+1);
            r[p].setName(name);
            model.add(r[p]);
        }
    }

    objFunction_uv(model,f,u,v);
    rest_Flow_uv(model,f, u, v); // flow const 
    rest_PathComp_uv(model,f,lambda, u,v); // every path with neg edges is pair
    restCycleMTZ_uv(model,f,r,u,v); // break cycle using MTZ
    
    // restCut_uv_1(model,f,u,v); // not repeat vertice in path



}

static void
createModel_Simplifed (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y){

    // add var x[u][j][s] in model
    for(int j=0; j<num_teams; j++) // for all project j 
        for(int u=0; u< num_vertices; u++) // for all worker u
            for(int s=0; s<num_skills; s++){ // 
                if(K[u][s] > 0 && R[j][s]>0){ // for all skills: s(u)
                    char name[30];
                    sprintf(name, "x%d%d%d", u+1, j+1,s+1);
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

    
    // print_matrix(G_Zuv);
    objFunction_Simplifed(model,y);

    rest_OneTeam(model, x); // max one team with one skill
    rest_MinSkill(model, x); // min skill s per team j
    rest_LinY(model,x,y); // linearization y with x
    rest_NegEdge(model,y); // negative edge incompatibility

}

static void
objFunction_uv (IloModel model, BoolVarMatrix f, int u, int v){
    IloEnv env = model.getEnv();
    IloExpr objExpr(env);
    for(int p=0; p<num_vertices; p++)
        for(int q=0; q<num_vertices; q++)
            if(G[u][v] < 1 && p!= q && G[p][q] != 0){ // (u,v) not in E+ and (p,q) is a possible arc between (u,v)-path
                objExpr += f[p][q];
            }

    IloObjective obj = IloMinimize(env, objExpr);
    model.add(obj);
    objExpr.end();

}

static void
objFunction_Simplifed (IloModel model, BoolVar3Matrix y){
    IloEnv env = model.getEnv();

    IloExpr objExpr(env);

    for(int u=0; u<num_vertices; u++)
        for(int v=u+1; v<num_vertices; v++)
            for(int j = 0; j<num_teams; j++){
                objExpr += G_Zuv[u][v]*y[u][v][j];  // (u,v) in E+
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
        model.add(expr <= 1);   // posso adcionar variavael boleana para saber sse sumConst está empty ou nao
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
                // model.add(y[u][v][j]  <= exprU);
                // model.add(y[u][v][j] <= exprV);
                exprU.end(); exprV.end();
            }

}

static void
rest_Flow_uv (IloModel model, BoolVarMatrix f, int u, int v){
    IloEnv env = model.getEnv();

    if(G[u][v] < 1){  // for all u<v and (u,v) not in E+
        for(int q = 0; q<num_vertices; q++){ // for fixed vertice q
            IloExpr sumArcs(env);
            for(int p = 0; p<num_vertices; p++) // for every arc (p,q) <=> (p,q) or (q,p) in E 
                if(p!=q && G[p][q] != 0){ //  (p,q) in E <=> (q,p) in E
                    sumArcs += f[p][q];
                    sumArcs += -f[q][p];
                }
            if(q == u){
                IloExpr expr(env);
                model.add(sumArcs == -1); expr.end();
            }else if(q == v){
                IloExpr expr(env);
                model.add(sumArcs == 1); expr.end();
            }else{
                model.add(sumArcs == 0);
            }
        }
    }

}

static void
rest_PathComp_uv (IloModel model, BoolVarMatrix f, IloIntVar lambda, int u, int v){
    IloEnv env = model.getEnv();
    
    if(G[u][v] < 1 ){ // (u,v) not in E+
        IloExpr expr(env);
        for(int p = 0; p<num_vertices; p++)
            for(int q = 0; q<num_vertices; q++)
                if(p!=q && G[p][q] < 0){ // sum of negatives arcs(A-)
                    expr += f[p][q];
                }
        model.add(expr - 2*lambda == 0); // pair number of negatives edges
        expr.end();
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
restCycleMTZ_uv (IloModel model, BoolVarMatrix f, IloNumVarArray r, int u, int v){

    IloEnv env = model.getEnv();
    if(G[u][v] < 1 ){ // (u,v) not in E+

        for(int p=0;p<num_vertices;p++)
            for(int q=0;q<num_vertices;q++){
                if(p!=q && abs(G[p][q])){
                    IloExpr expr(env);
                    expr = r[p] - r[q] + num_vertices*f[p][q];
                    model.add(expr <= (num_vertices-1));
                    expr.end();
                }
            }
    }
}

static void
restCut_uv_1 (IloModel model, BoolVarMatrix f, int u, int v){
    IloEnv env = model.getEnv();

    // pass only one time per vertice in a  u,v-flow
    for(int q = 0; q<num_vertices; q++){
        IloExpr expr(env);
        for(int p=0; p<num_vertices;p++)
            expr += f[p][q]; 
    
        model.add(expr<=1);
        expr.end(); 

    }
        
}

static void
displaySolution(IloCplex& cplex,
               BoolVar3Matrix x,
               BoolVar3Matrix y){



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

}

static void
    displaySolution_uv(IloCplex& cplex,
                   BoolVar4Matrix f,
                   IntVarMatrix lambda){

    cout << "--------------------------------------------------------" << endl;

    for(int u = 0; u<num_vertices; u++)
        for(int v = 0; v<num_vertices; v++)
            if(u<v && G[u][v] < 1){
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
    sprintf(arq, "%s/results/%d_Vertices_%d-%d-%d_TSP_SPC_simplified.ods",CURRENT_DIR, num_vertices, current_year, current_month, current_day);

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
            if(u<v && G[u][v] < 1){
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
    // sprintf(arq, "%s/results/%d_Vertices_2022-06-27_AL.ods",CURRENT_DIR, num_vertices);
    sprintf(arq, "%s/results/%d_Vertices_%d-%d-%d_TSP_SPC_subproblems.ods",CURRENT_DIR, num_vertices, current_year, current_month, current_day);
    // sprintf(arq, "%s/results/%d_Vertices_%d-%d-%d_TSP_SPC_algorithm.ods",CURRENT_DIR, num_vertices, current_year, current_month, current_day);
    
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
    compare_Zmatrices(){

        char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};


        instanceReader(0,7, 1,1, intances_types[0]);
        cout << "Create graph Algorithm" << endl;
        initG_Zuv();
        for(int u=0; u<num_vertices; u++)
            SCP_algorithm(u);


        int **G_Zuv_A;
        G_Zuv_A = (int**)(malloc(num_vertices*sizeof(int*)));
        for (int u=0;u<num_vertices;u++)
            G_Zuv_A[u] = (int*)malloc(num_vertices*sizeof(int));

        for (int u = 0; u < num_vertices; u++)
            for (int v = 0; v < num_vertices; v++){
                if(u==v)
                    G_Zuv_A[u][v] = 1;
                else
                    G_Zuv_A[u][v] = G_Zuv[u][v];
            }
                
        cout << "Create graph Subproblems " << endl;
        createG_Zuv();
        
        bool iguais = true;
        for (int u = 0; u < num_vertices; u++)
            for (int v = 0; v < num_vertices; v++){
                if(G_Zuv_A[u][v] != G_Zuv[u][v] ){
                    // cout << "Grafos diferentes pos: " << u << " , " <<  v << endl;
                    iguais = false;
                }
                
            }

        if(iguais)
            cout << "Grafos são iguais \n\n " << endl;
        else
            cout << "Grafos são diferentes \n\n " << endl;
        

        cout << "-------------------------------" << endl;

        for(int i=0; i<3; i++)
            for(int gclass = 1; gclass<4; gclass++){  // 1-3

                instanceReader(0,25, gclass,1, intances_types[i]);
                // instanceReader(0,7, 1,1, intances_types[i]);
                if(isConnected()){

                    cout << "Create graph Algorithm" << endl;
                    initG_Zuv();
                    for(int u=0; u<num_vertices; u++)
                        SCP_algorithm(u);


                    int **G_Zuv_A;
                    G_Zuv_A = (int**)(malloc(num_vertices*sizeof(int*)));
                    for (int u=0;u<num_vertices;u++)
                        G_Zuv_A[u] = (int*)malloc(num_vertices*sizeof(int));

                    for (int u = 0; u < num_vertices; u++)
                        for (int v = 0; v < num_vertices; v++){
                            if(u==v)
                                G_Zuv_A[u][v] = 1;
                            else
                                G_Zuv_A[u][v] = G_Zuv[u][v];
                        }
                            
                    cout << "Create graph Subproblems " << endl;
                    createG_Zuv();
                    
                    bool iguais = true;
                    for (int u = 0; u < num_vertices; u++)
                        for (int v = 0; v < num_vertices; v++){
                            if(G_Zuv_A[u][v] != G_Zuv[u][v] ){
                                // cout << "Grafos diferentes pos: " << u << " , " <<  v << endl;
                                iguais = false;
                            }
                            
                        }

                    if(iguais)
                        cout << "Grafos são iguais \n\n " << endl;
                    else
                        cout << "Grafos são diferentes \n\n " << endl;
                    

                    cout << "-------------------------------" << endl;
                    // print_matrix(G_Zuv_A);
                    // cout << "-------------------------------" << endl;
                    // print_matrix(G_Zuv);

                }else{
                    cout << "Graph is not connected "<< endl;
                    cout << "-------------------------------" << endl;
                }

        }



}

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
         if(!vis[i]) //if there is a node, not visited by traversal, G is not connected
            return false;
      }
   }
   return true;
}
static void
    testeGrafos(){

    char intances_types[3][12] =   {"random", "bitcoinotc", "epinions"};
    
    cout << "-------------------------------------------------" << endl;
    int vert = 25;
    for(int i=0; i<3; i++){
        for(int gclass = 1; gclass<4; gclass++){  // 1-3
            
            instanceReader(0,vert, gclass,1, intances_types[i]);
            
            if(isConnected())
                cout << "The G is connected." << endl;
            else
                cout << "The G is not connected." << endl;


            cout << "-------------------------------------------" << endl;
            
        }
    }
    
}

// make test PROB=TFP_SCP_simp METHOD=SIGN NUM_VERT=25 GCLASS=2
// make test PROB=TFP_SCP_simp METHOD=MTZ NUM_VERT=50
//  make test PROB=TFP_SCP_simp METHOD=SIGN NUM_VERT=50
//  make test PROB=TFP_SCP_simp METHOD=MinMatching NUM_VERT=50