/*
 * GeracaoColunas.h
 *
 *  Created on: 08/12/2012
 *      Author: luis
 */

#ifndef COLUMN_GENERATION_H_
#define COLUMN_GENERATION_H_

#include <list>
#include <map>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <set>
#include <fstream>
#include <queue>
#include <locale>
#include <iomanip>

#include <ilcplex/ilocplex.h>

#include "TypedefsAndDefines.h"
#include "Instance.h"
#include "Edge.h"
#include "Kruskals.h"
#include "TreeGlover.h"
#include "PathRelinking.h"
#include "Utils.h"
#include "MinCut.h"
#include "Dinic.h"
#include "VariableFixing.h"

using namespace std;

class Solution {

public:
	vector<Edge> spanning_tree;
	double obj_fun_value;

	Solution() {obj_fun_value = INF;}
	Solution(vector<Edge> st, double f) {spanning_tree = st; obj_fun_value = f;}
};

class Sec {

public:
	vector<int> S;
	int cont;

	Sec() { cont = 0; }
	Sec(vector<int> s) { S = s; cont = 0;}


};

class Blossom {

public:
	vector<int> H;
	vector<Edge> T;
	int cont;

	Blossom() { cont = 0; }
	Blossom(vector<int> h, vector<Edge> t) { H = h; T = t; cont = 0;}


};

class ColumnGeneration {

public:
	ColumnGeneration();

	~ColumnGeneration();

	void createVariables();

	void createObjectiveFunction();

	void createConstraints();

	void createModel();

	void generateInitialColumns();

	void getDualPrices(IloNumArray & u, IloNum & m, IloNumArray & p, IloNumArray & b);

	double pricing(const IloNumArray & u, const IloNum & m, const IloNumArray & p, IloNumArray & b, Columns & new_columns);

	void addNewColumns(const Columns & new_columns, int &con_cols);

	void blossomHeuristicEdgeWeight1(IloNumArray &xsol,vector<vector<int> > &Hs, vector<vector<Edge> > &Ts);

	void blossomHeuristicBridge(IloNumArray &xsol, vector<vector<int> > &Hs, vector<vector<Edge> > &Ts);

	bool findNewBlossoms();

	void addNewBlossom(vector<int> &H, vector<Edge> &T);

	void findNewSECs();

	void addNewSECs(vector<vector<int> > subsets);

	void setCplexParameters();

	Solution kruskalx();

	Solution lagrangeanPrimalHeuristic();

	Solution localSearchEdgeExchange(Solution s);		// First Improvement

//	Solution localSearchPrBased(Solution start, Solution target);

	int variableFixingLPR();

	int variableFixingLagr();

	Solution runPathRelinking();

	void primalSolution();

	void manageRows();

	void runCG(int time_limit);

	void printObjectives();

	void solveBc(ColumnGeneration &cg);

//private:
	IloEnv env_;
	IloModel model_;
	IloObjective objective_;
	IloCplex * cplex_;
	Graph* G;

	vector<Edge> all_edges;
	vector<Solution> elite_solutions;

	Solution best_primal;
	double best_dual;

	int num_secs;
	int num_blossoms;
	int cont_secs_rm;
	int cont_var_fixed;
	int fixed;

	double time_cg;

	ofstream debug;
	int cont_debug;

	// Variaveis que indicam se uma aresta consta ou nao na solucao
	IloNumVarArray x_;

	////// Restricoes //////
	// Restricoes de Grau
	IloRangeArray degree_constraints_;

	// Restrição que garante que a árvore é geradora
	IloRange st_constraint_;

	// Restrições SECs
	IloRangeArray secs_;

	// Blossoms
	IloRangeArray blossoms_;

	//Arestas já inseridas no pool de colunas
	vector<vector<int> > edges_in_pool;

	// Colunas
	Columns columns_;

	// Conjuntos S das retrições SEC's
	vector< vector<int> > S;

	vector<Sec*> secs;

	vector<Blossom*> blossoms;

	//int64_t GLOBAL_LABEL_ID_;

	int INITIAL_POOL_SIZE;

	double tot_start, tot_end;
	double start, end, tblossom, tsec, tprimal, tpric, tpl, tplsec;

	vector<int> bfs(int source, char checked[], vector<vector<int> > &lista_adj);
	vector<int> checa_solucao_bfs (vector<Edge> v);


};

#endif /* COLUMN_GENERATION_H_ */
