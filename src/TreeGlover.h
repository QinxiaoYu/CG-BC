/*
 * Estrutura para representar árvores geradoras proposta por Barr, Glover and Klingman (1979)
 * - p[i] representa o nó antecessor ao nó i na árvore (-1 se i for o nó raiz)
 * - t[i] representa o tamanho da subárvore enraizada por i (tam. da subárvore que contém i e todos os seus sucessores)
 * - a[i] é o custo da aresta de maior peso incidente ao vértice i
 * Os vetores s[] e f[] não foram considerados por serem desnecessários pro uso da estrutura neste trabalho
 * Para o entender a implementação dos métodos é necessário estudá-los com o artigo em mãos.
 * */

#ifndef TREEGLOVER_H_
#define TREEGLOVER_H_

#include <vector>
#include <stack>
#include <iostream>

#include "Edge.h"
#include "Instance.h"
#include "TypedefsAndDefines.h"

using namespace std;

class TreeGlover {

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
	bool t1_equal_tq;

	void setEdges(Edge ie, Edge oe);
	bool isInTree(Edge ie);

	void initializeTree(vector<vector<bool> > &adj_matrix);
	void dfs(stack<int> &S, vector<int> &discovered, vector<vector<bool> > &adj_matrix);

	Edge findBasisEquivalentPath(int op = 0);		// op = 1: o métedo irá retornar a aresta incidente à e.endpoint1 no ciclo
	// op = 2: o métedo irá retornar a aresta incidente à e.endpoint2 no ciclo
	void decidingT2();
	void updatingOperations();
	void reRootingT2();
	void attachT2T1();

	void print_values(vector<int> p, vector<int> t) const;
	void print_values_vertical(vector<int> p, vector<int> t) const;

public:

	TreeGlover(vector<Edge> &tree, double c);
	Edge performEdgeExchange(Edge ie, int c = 0);		// Retorna a aresta que saiu da solução. Se nenhuma troca foi feita, a aresta {0, 0} de custo 0 é retornada
	vector<Edge> getTree();
	double getTreeCost();
};

#endif
