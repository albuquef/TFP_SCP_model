#ifndef TFP_SCP_H
#define TFP_SCP_H
#include <iostream>
#include <filesystem>

#include "reader.h"
#include <ilcplex/ilocplex.h>

using namespace std;

//double runTime;
#define BILLION 1000000000L

ILOSTLBEGIN

//typedef
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix> BoolVar4Matrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVar3Matrix;

#include <sys/time.h>

class TFP_SCP
{

    public: 
        // TFP_SCP(Reader *r, const char* typeSEC):rd(r){rd->show();};
        TFP_SCP(Reader *r, const char* typeSEC);
        ~TFP_SCP();
        void exportILP      (IloCplex& cplex,const char *typeSEC);
        void solveILP       (void);
        void printSolution  (IloCplex& cplex,
                            BoolVar3Matrix x,
                            BoolVar3Matrix y,
                            BoolVar4Matrix f,
                            IntVarMatrix lambda);
        void saveSolution   (IloCplex& cplex,
                            BoolVar3Matrix x,
                            BoolVar3Matrix y,
                            BoolVar4Matrix f,
                            int class_type);
        // void saveResults    (IloCplex& cplex=this->cplex,
        //                     double timeF=0.0);
        void saveResults    (double timeF=0.0);
        void printInstance  (void);
        char instanceG[50];
        char instanceKR[50];
        char CURRENT_DIR[500];
        int current_day, current_month, current_year;
        const char* typeSEC;  
        double timeTFP_SCP;

    private:
        Reader *rd;
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        BoolVar3Matrix x;
        BoolVar3Matrix y;
        BoolVar4Matrix f;
        IntVarMatrix lambda;
        // void allocVars(env, x, y, f, lambda);

        void initILP        (const char* typeSEC);
        void initVariables  (void);
        void createModel    (IloModel model,
                            BoolVar3Matrix x,
                            BoolVar3Matrix y,
                            BoolVar4Matrix f,
                            IntVarMatrix lambda);
        
        NumVar3Matrix pi; // var mtz
        void createModel_MTZ(IloModel model, 
                            BoolVar3Matrix x,
                            BoolVar3Matrix y,
                            BoolVar4Matrix f, 
                            IntVarMatrix lambda, 
                            NumVar3Matrix pi);

        BoolVar3Matrix mu; // var sign
        void createModel_SIGN(IloModel model, 
                            BoolVar3Matrix x,
                            BoolVar3Matrix y,
                            BoolVar4Matrix f, 
                            IntVarMatrix lambda, 
                            BoolVar3Matrix mu);
        
        void
        objFunction (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

        void
        constr_OneTeam (IloModel model, BoolVar3Matrix x);

        void
        constr_MinSkill (IloModel model, BoolVar3Matrix x);

        void
        constr_LinY (IloModel model, BoolVar3Matrix x, BoolVar3Matrix y);

        void
        constr_Flow (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

        void
        constr_flowSameTeam (IloModel model, BoolVar3Matrix y, BoolVar4Matrix f);

        void
        constr_PathComp (IloModel model, BoolVar4Matrix f, IntVarMatrix lambda);

        void
        constr_NegEdge (IloModel model, BoolVar3Matrix y);

        // void
        // constr_Cuts_1 (IloModel model, BoolVar4Matrix f);

        // void
        // constr_Cuts_2 (IloModel model, BoolVar4Matrix f);

        void
        constr_MTZ (IloModel model, BoolVar4Matrix f, NumVar3Matrix pi);

        void
        constr_SIGN (IloModel model, BoolVar4Matrix f, BoolVar3Matrix mu);



};



#endif // TFP_SCP_H