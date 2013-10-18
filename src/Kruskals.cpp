/*
 * Kruskals.cpp
 *
 *  Created on: 31/05/2013
 *      Author: luishbicalho
 */

#include "Kruskals.h"

int cmpD(double a, double b){
	return (a <= b + EPS) ? (a + EPS < b) ? -1 : 0 : 1;
}

qUnion::qUnion(int n){
	id.resize(n);
	depth.resize(n);
	for(int i=0;i<n;i++) {
		id[i] = i;
		depth[i] = 1;
	}
}

// O(logn)
bool qUnion::verify(int v1, int v2) {
	return id[v1] == id[v2];
}

// O(logn)
void qUnion::push(int v1, int v2) {
	int v1id = id[v1];
	for(int i=0;i<Instance::getInstance()->num_vertex;i++)
		if(id[i]==v1id)
			id[i] = id[v2];
}

Graph::Graph(int n){
	resize(n);
}

// O(e)
void Graph::clear(int n) {
	if(n == -1)
		n = V();
	for(int i=0;i<n;i++)
		adj[i].clear();
}

// O(e)
void Graph::resize(int n){
	adj.resize(n);
}

// O(e)
bool Graph::adjacent(int p, int q){
	typename vector< pair<int, double> >::iterator it;
	for(it = adj[p].begin();it != adj[p].end();it++)
		if(it->first == q)
			return true;
	return false;
}

// O(1)
void Graph::push_edge(int p, int q, double w, bool dir){
	pair<int, double> E(q,w);
	adj[p].push_back(E);
	if(!dir){
		E.first = p;
		adj[q].push_back(E);
	}      // Verifica nao-direcionado.
}

// O(e)
void Graph::pop_edge(int p, int q, bool dir){
	typename vector< pair<int, double> >::iterator it;
	for(it = adj[p].begin();it != adj[p].end();it++)
		if(it->first == q){
			adj[p].erase(it);
			break;
		}
	if(dir)
		return; // Verifica nao-direcionado.
	for(it = adj[q].begin();it != adj[q].end();it++)
		if(it->first == p){
			adj[q].erase(it);
			break;
		}
}

// O(e)
double& Graph::edge(int p, int q){
	typename vector< pair<int, double> >::iterator it;
	for(it = adj[p].begin();it != adj[p].end();it++)
		if(it->first == q)
			return it->second;
	return it->second;
}

// O(1)
int Graph::E(int p){
	return adj[p].size();
}

// O(1)
int Graph::V(){
	return adj.size();
}

// O(E + V)
Graph Graph::invert(int n){
	typename vector< pair<int, double> >::iterator it;
	if(n == -1) n = V();
	Graph G(this->V());
	for(int i=0;i<n;i++)
		for(it = adj[i].begin();it != adj[i].end();it++){
			int p = i;
			int q = it->first;
			double w = it->second;
			G.push_edge(q,p,w,true);
		}
	return G;
}

// O(e)
vector< KEdge > Graph::get_edge(){
	int size = V();
	vector< KEdge > ret;
	typename vector< pair<int, double> >::iterator it;
	for(int i=0;i<size;i++)
		for(it = adj[i].begin();it != adj[i].end();it++)
			if(it->first >= i){
				KEdge E;
				E.p = i;
				E.q = it->first;
				E.w = it->second;
				ret.push_back(E);
			}
	return ret;
}

// O(1)
vector< pair<int, double> >& Graph::operator[](const int& a){
	return adj[a];
}

bool operator<(KEdge p, KEdge q){
	return p.w < q.w;
}

vector< KEdge > kruskal(Graph* G) {
	vector< KEdge > E = G->get_edge();
	sort(E.begin(), E.end());
	qUnion qunion(Instance::getInstance()->num_vertex);
	int size = E.size();
	vector< KEdge > mst;
	for(int i=0; i<size; i++){
		if(qunion.verify(E[i].p, E[i].q))
			continue;
		mst.push_back(E[i]);
		qunion.push(E[i].p, E[i].q);
	}

	return mst;
}

vector< KEdge > kruskalDegreeConstrained(Graph* G) {
	vector< KEdge > E = G->get_edge();
	sort(E.begin(), E.end());
	qUnion qunion(Instance::getInstance()->num_vertex);
	int size = E.size();
	vector< KEdge > mst;

	int numv = Instance::getInstance()->degree_constraint.size();
	vector<int> contg(numv, 0);

	for(int i=0; i<size; i++) {
		if(contg[E[i].p] <= Instance::getInstance()->degree_constraint[E[i].p]-1 && contg[E[i].q] <= Instance::getInstance()->degree_constraint[E[i].q] - 1) {
			if(qunion.verify(E[i].p, E[i].q))
				continue;

			if((int)mst.size() == Instance::getInstance()->num_vertex - 2) {
				mst.push_back(E[i]);
				qunion.push(E[i].p, E[i].q);
				++contg[E[i].p];
				++contg[E[i].q];
				break;
			}
			vector<int> deltat(numv, 0);
			vector<int> dt(numv, 0);
			for(int j=0; j<Instance::getInstance()->num_vertex; j++) {
				dt[qunion.id[j]] += Instance::getInstance()->degree_constraint[j];
			}
			for(int j=0; j<(int) mst.size(); j++) {
				deltat[qunion.id[mst[j].p]]++;
				deltat[qunion.id[mst[j].q]]++;
			}
			if(dt[qunion.id[E[i].p]] - deltat[qunion.id[E[i].p]] > 1 || dt[qunion.id[E[i].q]] - deltat[qunion.id[E[i].q]] > 1) {
				mst.push_back(E[i]);
				qunion.push(E[i].p, E[i].q);
				++contg[E[i].p];
				++contg[E[i].q];
			}
		}
	}

	return mst;
}

int heuristicKruskalSEC(Graph* G, IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector<vector<int> > &secs) {

	vector< KEdge > E = G->get_edge();
	sort(E.begin(), E.end());
	qUnion qunion(Instance::getInstance()->num_vertex);
	int size = E.size();
	vector< KEdge > mst;

	double most_violated_value = INF;
	int id_most_violated = -1;
	int cont_edges = 0;
	int vref;

	for(int i=0; i < size; i++) {
		if(qunion.verify(E[i].p, E[i].q))
			continue;
		mst.push_back(E[i]);
		qunion.push(E[i].p, E[i].q);

		++cont_edges;
		if(cont_edges > 1) {
			vref = E[i].q;		// vértice referência para procurar quais estão na mesma componente conexa
			vector<int> S(Instance::getInstance()->num_vertex, 0);
			for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
				if(qunion.verify(i,vref) == true)	// O vértice i está na mesma componente conexa de vref?
					S[i] = 1;						// Inclui i no corte
			}

			double sum = 0;
			int cardS = 0;
			for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
				if(S[i] == 1) {	// O vertice i pertence a S?
					cardS++;		// Adiciona a unidade referente ao elemento i pertencente à S
					for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
						if(S[j] == 1 && edges_in_pool[i][j] != -1) {		// Ambos os vertices, i e j, pertencem à S e {i, j} aresta está presente no pool de colunas?
							sum += xsol[edges_in_pool[i][j]];		// Soma o valor de relaxacao linear referente à aresta {i, j}
						}
					}
				}
			}

			if(sum - double(cardS-1) > 0.00001 && cardS > 1 && cardS < Instance::getInstance()->num_vertex) {		// O corte atual está violado?
				secs.push_back(S);
				if(sum - double(cardS-1) < most_violated_value) {		// O corte atual é mais violado de todos até o momento?
					id_most_violated = (int)secs.size() - 1;
					most_violated_value = sum - double(cardS-1);
				}
			}
		}
	}

	return id_most_violated;
}

vector<KEdge> heuristicKruskalPrBased(Graph* G, const vector<Edge> &start, const vector<Edge> &target) {

    int nv = Instance::getInstance()->num_vertex;

    qUnion qunion(nv);

    vector<int> contg(nv, 0);

    vector<vector<int> > sol_start(nv);
    vector<vector<int> > sol_target(nv);
    vector<vector<int> > sol_heur(nv);

    vector< KEdge > E;
    vector< KEdge > mst;

    sol_heur.resize(nv);
    for(int i=0; i < nv; i++) {
        sol_start[i].resize(nv, 0);
        sol_target[i].resize(nv, 0);
        sol_heur[i].resize(nv, 0);
    }

    for(int i=0; i < nv - 1; i++) {
        sol_start[start[i].endpoint1][start[i].endpoint2] = sol_start[start[i].endpoint2][start[i].endpoint1] = 1;
        sol_target[target[i].endpoint1][target[i].endpoint2] = sol_target[target[i].endpoint2][target[i].endpoint1] = 1;
    }

    for(int i=0; i < nv; i++) {
        for(int j=i+1; j < nv; j++) {
            KEdge e(i, j, Instance::getInstance()->cost_matrix[i][j]);
            if(sol_start[i][j] == sol_target[i][j] &&  sol_start[i][j] == 1) {
                sol_heur[i][j] = sol_heur[j][i] = 1;
                mst.push_back(e);
                qunion.push(e.p, e.q);
                ++contg[e.p];
                ++contg[e.q];
            }
            else {
                sol_heur[i][j] = sol_heur[j][i] = 0;
                E.push_back(e);
            }
        }
    }

    int size = E.size();
    sort(E.begin(), E.end());

    for(int i=0; i < size; i++) {
        if(contg[E[i].p] <= Instance::getInstance()->degree_constraint[E[i].p]-1 && contg[E[i].q] <= Instance::getInstance()->degree_constraint[E[i].q] - 1) {
            if(qunion.verify(E[i].p, E[i].q)) {
                continue;
            }

            if((int)mst.size() == nv - 2) {
                mst.push_back(E[i]);
                qunion.push(E[i].p, E[i].q);
                ++contg[E[i].p];
                ++contg[E[i].q];
                break;
            }

            vector<int> deltat(nv, 0);
            vector<int> dt(nv, 0);
            for(int j=0; j < nv; j++) {
                dt[qunion.id[j]] += Instance::getInstance()->degree_constraint[j];
            }
            for(int j=0; j<(int) mst.size(); j++) {
                deltat[qunion.id[mst[j].p]]++;
                deltat[qunion.id[mst[j].q]]++;
            }
            if(dt[qunion.id[E[i].p]] - deltat[qunion.id[E[i].p]] > 1 || dt[qunion.id[E[i].q]] - deltat[qunion.id[E[i].q]] > 1) {
                mst.push_back(E[i]);
                qunion.push(E[i].p, E[i].q);
                ++contg[E[i].p];
                ++contg[E[i].q];
            }
        }
    }

    return mst;
}
