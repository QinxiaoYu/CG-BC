/*
 * PathRelinking.h
 *
 *  Created on: 12/08/2013
 *      Author: luishbicalho
 */

#ifndef PATHRELINKING_H_
#define PATHRELINKING_H_

#include <vector>
#include <stack>
#include <queue>
#include <algorithm>

#include "Instance.h"
#include "Edge.h"
#include "TypedefsAndDefines.h"

using namespace std;

class TreeGloverModPR {
private:
	double tree_cost;
	int num_vertex;
	vector<vector<bool> > adj_matrix;
	vector<int> px;
	vector<int> tx;
	vector<double> ax;
	int root;
	int x1, x2;
	int y1, y2;
	int u, v;
	int p, q;
	int z;
	Edge incoming_edge;
	Edge outgoing_edge;

	void setEdges(Edge ie, Edge oe);
	bool isInTree(Edge ie);

	void initializeTree(vector<vector<bool> > &adj_matrix);
	void dfs(stack<int> &S, vector<int> &discovered, vector<vector<bool> > &adj_matrix);

	Edge findBasisEquivalentPath(vector<vector<int> > &sol, int op = 0);		// op = 1: o métedo irá retornar a aresta incidente à e.endpoint1 no ciclo
	                                                                            // op = 2: o métedo irá retornar a aresta incidente à e.endpoint2 no ciclo
	void decidingT2();
	void updatingOperations();
	void reRootingT2();
	void attachT2T1();

	void print_values(vector<int> p, vector<int> t) const;
	void print_values_vertical(vector<int> p, vector<int> t) const;

public:
	TreeGloverModPR();
	TreeGloverModPR(vector<Edge> &tree, double c);
	Edge performEdgeExchange(Edge ie, vector<vector<int> > &sol, int c = 0);
	vector<Edge> getTree();
	double getTreeCost();
};

class PathRelinking {

public:
	vector<vector<int> > sol_start;
	vector<vector<int> > sol_target;
	vector<vector<int> > sol;
	vector<Edge> start;
	vector<Edge> target;
	TreeGloverModPR *tree;
	vector<int> degree;

	PathRelinking(vector<Edge> &start, vector<Edge> &target);
	~PathRelinking();
	vector<Edge> runPR();

};

#endif /* PATHRELINKING_H_ */
