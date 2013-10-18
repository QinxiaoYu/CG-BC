/*
 * Dinic.cpp
 *
 *  Created on: 09/06/2013
 *      Author: luishbicalho
 */

#include "Dinic.h"

void Dinic::addEdge(int u, int v, int cap) {
	if(u == v) return;
	EdgeD e(v, cap, int(g_adj_list[v].size()));
	EdgeD r(u, 0, int(g_adj_list[u].size()));
	g_adj_list[u].push_back(e);
	g_adj_list[v].push_back(r);
}

bool Dinic::buildLevelGraph(int src, int sink) {
	fill(ALL(level), -1);
	while(not q.empty()) q.pop();
	level[src] = 0;
	q.push(src);
	while(not q.empty()){
		int u = q.front();
		q.pop();
		FOREACH(e, g_adj_list[u]) {
			if(not e->cap or level[e->v] != -1) continue;
			level[e->v] = level[u] + 1;
			if(e->v == sink) return true;
			q.push(e->v);
		}
	}
	return false;
}

int Dinic::blockingFlow(int u, int sink, int f) {
	if(u == sink or not f) return f;
	int fu = f;
	FOREACH(e, g_adj_list[u]){
		if(not e->cap or level[e->v] != level[u] + 1) continue;
		int mincap = blockingFlow(e->v, sink, min(fu, e->cap));
		if(mincap){
			g_adj_list[e->v][e->rev].cap += mincap;
			e->cap -= mincap;
			fu -= mincap;
		}
	}
	if(f == fu) level[u] = -1;
	return f - fu;
}

int Dinic::maxFlow(int src, int sink) {
	flow = 0;
	while(buildLevelGraph(src, sink))
		flow += blockingFlow(src, sink, numeric_limits<int>::max());
	return flow;
}

vector<int> Dinic::minCut() {

	vector<int> mc(n, 0);		// vetor de tamanho num_vertex. 1 se o vértice pertence à S, 0 caso contrário.
	for(int i=0; i<n; i++)
		if((*this).level[i] != -1)
			mc[i] = 1;
	return mc;

}


