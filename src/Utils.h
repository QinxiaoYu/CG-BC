#ifndef UTIL_H
#define UTIL_H

#include <limits>
#include <vector>
#include <list>
#include <queue>
#include <cmath>

#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <sys/unistd.h>
#include <ilcplex/ilocplex.h>

#include "Instance.h"
#include "Edge.h"
#include "TypedefsAndDefines.h"

using std::vector;
using std::queue;
using std::list;

namespace Utils {

/// funcoes

// para comparacoes entre doubles
int cmp(double x, double y = 0, double tol = 1e-6);

double roundPrecision(double value);

// verifica se A eh um subconjunto de B, onde A armazena os seus elementos e B eh um vetor de booleanos indicando elementos pertencentes ao conjunto
bool isSubset(const std::vector<int> & A, const std::vector<int> & B);

double euclidianNorm( vector<int> &v, IloCplex * master_, IloNumArray &x_, vector<vector<int> > &edges_in_pool);

double innerProduct(vector<int> v1, vector<int> v2);

// Entrada: a lista de adjacência do grafo
// Armazena em H os conjuntos Handles obtidos através de uma busca em largura
void findHandles(vector<vector<int> > &adj_list, vector<vector<int> > &H);

// Entrada: lista de adjacência do grafo, um vetor indicando quais vértices foram visitados e o nó raiz pra iniciar a busca
// Saída: um vetor de tamanho num_vertex onde o vértice pertencente ao handle é sinalizado por 1
vector<int> bfsHandle(vector<vector<int> > &adj_list, char checked[], int source);

// -------- http://www.geeksforgeeks.org/bridge-in-a-graph/
// A recursive function that finds and store bridges using DFS traversal
// u --> The vertex to be visited next
// visited[] --> keeps tract of visited vertices
// disc[] --> Stores discovery times of visited vertices
// parent[] --> Stores parent vertices in DFS tree
// adj --> Adjacency list
// bridges --> vector that stores the bridges
void bridgeUtil(int u, list<int> *adj, bool visited[], int disc[], int low[], int parent[], vector<vector<int> > &bridges);

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void findBridges(vector<Edge> &tree, vector<vector<int> > &bridges);		// bridges é uma matriz V x V inicializada com 0 em todas as posições. Ao final, as pontes são os
// as posições [i][j] cujo valor é 1

void bfs(int source, char checked[], vector<vector<int> > &lista_adj);

void checa_solucao_bfs (vector<Edge> v);

// retorna tempo
double getTime();

}

inline double Utils::getTime() {
	struct tms time_;
	clock_t time_tic;
	double cpt = sysconf(_SC_CLK_TCK);
	times(&time_);
	time_tic = time_.tms_utime;
	return (double) (time_tic / cpt);
}

#endif // UTIL_H
