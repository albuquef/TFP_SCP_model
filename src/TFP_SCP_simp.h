#ifndef TFP_SCP_SIMP_H
#define TFP_SCP_SIMP_H

#include <iostream>
#include <filesystem>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


#include "reader.h"
#include "graph.h"

using namespace std;

//typedef
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3Matrix;

#include <sys/time.h>


class TFP_SCP_SIMP{

public:
    TFP_SCP_SIMP(Reader *r, const Graph & Grph,const char* method_SPP);
    ~TFP_SCP_SIMP();
    void exportILP      (IloCplex& cplex,const char* method_SPP);
    void solveILP       (void);
    void printSolution  (IloCplex& cplex,
                        BoolVar3Matrix x,
                        BoolVar3Matrix y);
    void saveSolution   (IloCplex& cplex,
                        BoolVar3Matrix x,
                        BoolVar3Matrix y,
                        int class_type);
    // void saveResults    (IloCplex& cplex,
    //                     double timeCplex, double time_WeightedGraph);
    void saveResults    (double timeTotal);
    void printInstance  (void);
    char instanceG[50];
    char instanceKR[50];
    char CURRENT_DIR[500];
    int current_day, current_month, current_year;
    const char* method_SPP;
    double timeTFP;
    double time_GraphSPP;
    

private:
    Reader *rd;
    Graph GRAPH;
    // int** Weighted_Graph;
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    BoolVar3Matrix x;
    BoolVar3Matrix y;

    void initModel      (const char* method_SPP);
    void initVariables  (void);
    void createModel    (IloModel model,
                        BoolVar3Matrix x,
                        BoolVar3Matrix y);
      
    void
    objFunction (IloModel model, BoolVar3Matrix y);

    void
    constr_OneTeam (IloModel model, BoolVar3Matrix x);

    void
    constr_MinSkill (IloModel model, BoolVar3Matrix x);

    void
    constr_LinY (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

    void
    constr_Incomp (IloModel model, BoolVar3Matrix x); // u and v not comp


    void create_GraphSPP_ILP();
    void allocVars_SPP(IloEnv env_MTZ,BoolVarMatrix f);

    void createModel_MTZ(IloModel model_MTZ,BoolVarMatrix f,IloIntVar lambda,IloNumVarArray pi,int u, int v);  
    void createModel_SIGN(IloModel model_MTZ,BoolVarMatrix f,IloIntVar lambda,IloBoolVarArray mu,int u, int v);  


    void objFunction_SPP_uv (IloModel model_SPP,BoolVarMatrix f, int u, int v);
    void constr_Flow_uv(IloModel model_SPP,BoolVarMatrix f, int u, int v);
    void constr_PathComp_uv (IloModel model_SPP, BoolVarMatrix f, IloIntVar lambda, int u, int v);
    void constr_BreakCycle_MTZ_uv(IloModel model_SPP, BoolVarMatrix f, IloNumVarArray pi, int u, int v);
    void constr_BreakCycle_SIGN_uv(IloModel model_SPP, BoolVarMatrix f, IloBoolVarArray mu, int u, int v);

};



#endif // TFP_SCP_SIMP_H