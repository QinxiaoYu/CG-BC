/*
 * BranchAndCut.h
 *
 *  Created on: 01/10/2013
 *      Author: luishbicalho
 */

#ifndef BRANCHANDCUT_H_
#define BRANCHANDCUT_H_

#include <ilcplex/ilocplex.h>
#include "ColumnGeneration.h"
#include "Utils.h"
#include "LocalSearch.h"

class BranchAndCut {

    public:
        IloEnv bc_env_;
        IloModel bc_model_;
        IloObjective bc_objective_;
        IloCplex * bc_cplex_;

        vector<Edge> columns_;
        vector<vector<int> > edges_in_pool;

        // Variaveis que indicam se uma aresta consta ou nao na solucao
        IloNumVarArray bc_x_;

        ////// Restricoes //////
        IloRangeArray constraints_;

        BranchAndCut(ColumnGeneration &cg);

        void createVariables();

        void createObjectiveFunction();

        void createConstraints(ColumnGeneration &cg);

        void setCplexParameters();

        void getCgBasis(ColumnGeneration &cg);

        void solveBc(double time_limit);



};


#endif /* BRANCHANDCUT_H_ */
