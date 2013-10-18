/*
 * VariableFixing.h
 *
 *  Created on: 19/09/2013
 *      Author: luishbicalho
 */

#ifndef VARIABLEFIXING_H_
#define VARIABLEFIXING_H_

#include <vector>
#include <stack>
#include <ilcplex/ilocplex.h>
#include "Utils.h"
#include "Edge.h"
#include "TypedefsAndDefines.h"
#include "Instance.h"
#include "Kruskals.h"

using std::vector;
using std::stack;

class VariableFixing {

    private:
        int num_vertex;
        int root;
        double lb;
        vector<vector<bool> > adj_list;
        vector<int> px;
        vector<int> tx;
        vector<vector<double> > cost_matrix;

        void initializeTree(vector<vector<bool> > &adj_matrix);
        void dfs(stack<int> &S, vector<int> &discovered, vector<vector<bool> > &adj_matrix);

        Edge findLargestCostEdge(int u, int v);

    public:

        VariableFixing(Graph *G, vector<vector<double> > &lagr_cost_matrix, double c);
        int runVariableFixing(double ub, vector<vector<int> > &edges_in_pool, vector<Edge> &all_edges, vector<Edge*> &columns_);
};


#endif /* VARIABLEFIXING_H_ */
