/*
 * MinCut.h
 *  Created on: 05/05/2013
 *      Author: luishbicalho
 */

#ifndef MINCUT_H_
#define MINCUT_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ilcplex/ilocplex.h>

#include "Utils.h"
#include "Instance.h"
#include "Dinic.h"
#include "TypedefsAndDefines.h"

using namespace std;
/*
 *  --- Classe MinCut ---
 *  Classe que implementa o algoritmo proposto em Padberg & Wolsey (1983) para encontrar
 *  cortes mínimos para SECs utilizando (|V| - 2) problemas de fluxo máximo (Dinic.h)
 * */
class MinCut {

public:
	int num_vertex;         // número de vértices
	int source;             // nó fonte (vértice |n|)
	int sink;               // nó sorvedouro (vértice |n|+1)

	MinCut();
	~MinCut();

	/*
	 *  Encontra cortes violados
	 *  Input:
     *             xsol: valor de relaxação linear das arestas presente no modelo em um dado momento do problema mestre
     *    edges_in_pool: estrutura auxiliar no qual [i][j] representa a posição que a aresta {i, j} está no vetor com os
     *                   valores de relaxação linear (xsol)
     *  Output:
     *               int: índice do corte mais violado em min_cuts
     *          min_cuts: vetor com os subconjuntos S que geram uma desigualdade SEC violada
	 * */
	int solveMinCutWolsey(IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector<vector<int> > &min_cuts);

    /*
     *  Encontra cortes violados que são ortogonais a um corte violado (vec).
     *  As arestas presentes em vec são removidas do grafo suporte e um novos problemas de fluxo máximo são
     *  resolvidos a fim de se obter SECs violadas.
     *  Input:
     *             xsol: valor de relaxação linear das arestas presente no modelo em um dado momento do problema mestre
     *    edges_in_pool: estrutura auxiliar no qual [i][j] representa a posição que a aresta {i, j} está no vetor com os
     *                   valores de relaxação linear (xsol)
     *              vec: conjunto S de uma SEC violada
     *  Output:
     *               int: índice do corte mais violado em min_cuts
     *          min_cuts: vetor com os subconjuntos S que geram uma desigualdade SEC violada
     * */
	int solveMinCutWolseyOrth(IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector<int> &vec, vector<vector<int> > &min_cuts);

};



#endif /* FLUXO_H_ */
