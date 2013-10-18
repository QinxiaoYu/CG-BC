/*
 * Edge.h
 *
 *  Created on: 04/05/2013
 *      Author: luishbicalho
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <ilcplex/ilocplex.h>

class Edge {

public:
	int endpoint1;          // uma das extremidades da aresta (endpoint1 < endpoint2)
	int endpoint2;          // outra extremidade da aresta
	double cost;            // custo da aresta
	double cost_reduced;    // custo reduzido da aresta no problema de pricing
	int cont;               // número de iterações que a aresta ficou fora da base
	bool fixed;             // indica se a aresta está fixada ou não

	Edge();

	Edge(int u, int v, double c);
	Edge(int u, int v, double c, IloNumVar var);
	void setEdge(int u, int v, double c);


	/*
	 * Os operadores foram sobrecarregados tomando como base o atributo cost_reduced
	 * */
	bool operator<(const Edge & e) const;
	bool operator>(const Edge & e) const;
	bool operator<=(const Edge & e) const;
	bool operator>=(const Edge & e) const;
	bool operator==(const Edge & e) const;
	bool operator!=(const Edge & e) const;

};


#endif /* EDGE_H_ */
