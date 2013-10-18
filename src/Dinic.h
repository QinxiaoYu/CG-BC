/*
 * Dinic.h
 *
 *  Created on: 09/06/2013
 *      Author: Dilson Guimarães / Luis H. C. Bicalho
 *
 *      Implementação do algoritmo de Dinitz para encontrar o fluxo máximo (corte mínimo) em um grafo direcionado.
 *
 */

#ifndef DINIC_H_
#define DINIC_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

#include "TypedefsAndDefines.h"

using namespace std;

/* --- Classe EdgeD ---
 *  Representa uma aresta no algoritmo de Dinitz.
 *  Obs: o vértice de origem da aresta está implicito no índice do vetor que representa a lista de adjacência do grafo.
 * */
class EdgeD {

public:
	int v;		// vértice de destino da aresta
	int rev;	// índice na lista de adjacência do vértice v da aresta "que retorna" no grafo de suporte
	int cap;	// capacidade da aresta (quantidade máxima de fluxo que pode passar na aresta)
	EdgeD(int v_, int cap_, int rev_) : v(v_), rev(rev_), cap(cap_) {}

};

/*
 *  --- Classe Dinic ---
 *  Classe que implementa o método de Dinitz para encontrar
 *  o fluxo máximo (através deste o corte mínimo) em um grafo direcionado
 */
class Dinic {

public:
	vector<vector<EdgeD> > g_adj_list;		// lista da adjacência que representa o grafo
	vector<int> level;						// level[i] indica o nível do vértice i em relação ao source
	queue<int> q;							// fila necessária para construir o grafo de nível
	int flow, n;

	Dinic(int n_) : g_adj_list(n_), level(n_), n(n_) { flow = 0; }

	/*  Adiciona uma arco ao grafo. Ao adicionar um arco (u, v) um outro arco (v, u) é criado com cap = 0
     *  Input:
     * 	 	u: extremidade de origem da aresta
     *		v: extremidade de destino da aresta
     *	  cap: capacidade suportada pela aresta (u, v)
     */
	void addEdge(int u, int v, int cap);

	/*
	 *  Constrói o grafo de nível através de uma modificação da busca em largura
	 *  Input:
	 *    src: nó fonte de onde o fluxo inicial deve sair
	 *   sink: nó sorvedouro no qual o fluxo deve chegar
	 *
	 *  Output:
	 *    true, se foi possível construir um grafo de nível partindo de src até sink
	 *   false, caso contrário (neste caso o algoritmo para encontrar fluxo se encerra pois não existe nenhum caminho aumentante)
	 * */
	bool buildLevelGraph(int src, int sink);

	/*
	 *  Encontra um ou vários caminhos aumentantes através do grafo de nível construído por buildLevelGraph
	 *  Input:
     *    src: nó fonte de onde o fluxo inicial deve sair
     *   sink: nó sorvedouro no qual o fluxo deve chegar
     *      f: fluxo máximo que pode sair de u até sink através de uma das arestas cuja extremidade inicial é u
     *
     *  Output:
     *      Quantidade máxima de fluxo que pode passar num dado grafo de nível
	 * */
	int blockingFlow(int u, int sink, int f);

	/*
	 *  Encontra o fluxo máximo suportado pelo grafo reprensentado por g_adj_list
	 *  Input:
     *    src: nó fonte de onde o fluxo inicial deve sair
     *   sink: nó sorvedouro no qual o fluxo deve chegar
     *
     *  Output:
     *      Quantidade máxima de fluxo suportado pelo grafo g_adj_list tomando como fonte src e como soverdouro sink
	 * */
	int maxFlow(int src, int sink);

	/*
	 *  Encontra o corte mínimo do grafo g_adj_list cujo valor é igual ao retornado por maxFlow
	 *  Output:
	 *      Vetor com n (número de vértices) posições no qual vec[i] = 1 indica que o
	 *      vértice i está no "lado esquerdo" do corte e vec[i] = 0, c.c.
	 *      Ou seja, vec[i] = 1 indica que o vértice i pode ser alcançado através de src se
	 *      as arestas do corte mínimo forem suprimidas do grafo
	 * */
	vector<int> minCut();

};

#endif /* DINIC_H_ */
