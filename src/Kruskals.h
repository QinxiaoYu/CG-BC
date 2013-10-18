/*
 * Kruskals.h
 *
 *  Created on: 31/05/2013
 *      Author: luishbicalho
 */

#ifndef KRUSKALS_H_
#define KRUSKALS_H_

#include <algorithm>
#include <iostream>
#include <vector>
#include <ilcplex/ilocplex.h>

#include "Instance.h"
#include "Utils.h"
#include "Edge.h"
#include "TypedefsAndDefines.h"

int cmpD(double a, double b);

using namespace std;

// Vertex conection
class qUnion {

public:

	vector< int > id; vector< int > depth; vector< int > id2;

	qUnion(int n);

	// O(logn)
	bool verify(int v1, int v2);

	// O(logn)
	void push(int v1, int v2);

};

class KEdge {

public:

	int p;

	int q;

	double w;

	KEdge() { }
	KEdge(int p_, int q_, double w_) : p(p_) , q(q_) , w(w_) { }

};

class Graph {

public:

	vector< vector< pair<int, double> > > adj;

	Graph(int n = 0);

	// O(e)
	void clear(int n = -1);

	// O(e)
	void resize(int n);

	// O(e)
	bool adjacent(int p, int q);

	// O(1)
	void push_edge(int p, int q, double w, bool dir = false);

	// O(e)
	void pop_edge(int p, int q, bool dir = false);

	// O(e)
	double& edge(int p, int q);

	// O(1)
	int E(int p);

	// O(1)
	int V();

	// O(E + V)
	Graph invert(int n = -1);

	// O(e)
	vector< KEdge > get_edge();

	// O(1)
	vector< pair<int, double> >& operator[](const int& a);

};

bool operator<(KEdge p, KEdge q);

/*
 *  Algortimo de kruskal para encontrar a árvore geradora mínima de um grafo.
 *  Input:
 *      G: grafo (representado por uma matriz de adjacência)
 * */
vector< KEdge > kruskal(Graph* G);

/*
 *  Heurística para encontrar uma árvore geradora que não viola as restrições de grau
 *  Input:
 *      G: grafo (representado por uma matriz de adjacência)
 * */
vector< KEdge > kruskalDegreeConstrained(Graph* G);

/*
 *  Heurística para encontrar SEC's violadas através da construção de uma árvore geradora mínima de G
 *  A cada inclusão de uma aresta {i, j} no algoritmo de kruskal é verificado se os vértices contido no componente
 *  conexo que contém {i, j} gera uma desigualdade SEC violada, através dos valores de relaxação linear da solução
 *  corrente no problema mestre.
 *  Input:
 *                G: grafo
 *             xsol: valor de relaxação linear das arestas presente no modelo em um dado momento do problema mestre
 *    edges_in_pool: estrutura auxiliar no qual [i][j] representa a posição que a aresta {i, j} está no vetor com os
 *                   valores de relaxação linear (xsol)
 *  Output:
 *      secs: vetor com os subconjuntos S que geram uma desigualdade SEC violada
 *       int: índice do vetor sec onde se encontra a desigualdade mais violada
 * */
int heuristicKruskalSEC(Graph* G, IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector< vector<int> > &secs);

vector<KEdge> heuristicKruskalPrBased(Graph* G, const vector<Edge> &start, const vector<Edge> &target);

#endif /* KRUSKALS_H_ */
