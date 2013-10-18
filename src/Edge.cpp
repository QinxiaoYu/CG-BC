/*
 * Edge.cpp
 *
 *  Created on: 04/05/2013
 *      Author: luishbicalho
 */

#include "Edge.h"

Edge::Edge() {
	endpoint1 = endpoint2 = cost = cost_reduced = cont = fixed = 0;
}

Edge::Edge(int u, int v, double c) {
	setEdge(u, v, c);
}

Edge::Edge(int u, int v, double c, IloNumVar var) {
	setEdge(u, v, c);
}

void Edge::setEdge(int u, int v, double c) {
	endpoint1 = u;
	endpoint2 = v;
	cost = c;
	cost_reduced = c;
	cont = 0;
	fixed = false;
}

bool Edge::operator<(const Edge & e) const {
	return cost_reduced < e.cost_reduced;
}

bool Edge::operator<=(const Edge & e) const {
	return cost_reduced <= e.cost_reduced;
}

bool Edge::operator==(const Edge & e) const {
	return endpoint1 == e.endpoint1 && endpoint2 == e.endpoint2 && cost_reduced == e.cost_reduced;
}

bool Edge::operator>(const Edge & e) const {
	return cost_reduced > e.cost_reduced;
}

bool Edge::operator>=(const Edge & e) const {
	return cost_reduced >= e.cost_reduced;
}

bool Edge::operator!=(const Edge & e) const {
	return endpoint1 != e.endpoint1 || endpoint2 != e.endpoint2 || cost_reduced != e.cost_reduced;
}
