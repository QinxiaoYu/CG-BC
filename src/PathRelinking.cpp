/*
 * PathRelinking.cpp
 *
 *  Created on: 12/08/2013
 *      Author: luishbicalho
 */

#include "PathRelinking.h"

TreeGloverModPR::TreeGloverModPR(vector<Edge> &tree, double c) {

	root = 0;
	num_vertex = (int)tree.size()+1;
	tree_cost = c;

	ax.resize(num_vertex, 0);
	px.resize(num_vertex);
	tx.resize(num_vertex, 1);

	adj_matrix.resize(num_vertex);
	for(int i=0; i<num_vertex; i++)
		adj_matrix[i].resize(num_vertex, 0);

	for(unsigned int i = 0; i < tree.size(); i++) {
		adj_matrix[tree[i].endpoint1][tree[i].endpoint2] = true;
		adj_matrix[tree[i].endpoint2][tree[i].endpoint1] = true;
		if(ax[tree[i].endpoint1] < tree[i].cost)
			ax[tree[i].endpoint1] = tree[i].cost;
		if(ax[tree[i].endpoint2] < tree[i].cost)
			ax[tree[i].endpoint2] = tree[i].cost;
	}

	initializeTree(adj_matrix);

}

void TreeGloverModPR::dfs(stack<int> &S, vector<int> &discovered, vector<vector<bool> > &adj_list) {

	while(S.size() > 0) {
		int a = S.top();
		for(int b = 0; b < num_vertex; b++) {
			if(adj_list[a][b] == true) {
				if(discovered[b] == 0) {
					S.push(b);
					discovered[b] = 1;
					px[b] = a;
					dfs(S, discovered, adj_list);
					tx[a] += tx[b];
				}
			}
		}
		S.pop();
		return;
	}
}

void TreeGloverModPR::initializeTree(vector<vector<bool> > &adj_list) {

	stack<int> S;
	vector<int> discovered(num_vertex, 0);
	root = 0;
	px[root] = -1;
	discovered[root] = 1;
	S.push(root);
	dfs(S, discovered, adj_list);

}

void TreeGloverModPR::print_values(vector<int> p, vector<int> t) const {

	cout << "p: ";
	for(int i=0; i<(int)p.size(); i++) {
		cout << p[i] << " ";
	}
	cout << endl;
	cout << "t: ";
	for(int i=0; i<(int)t.size(); i++) {
		cout << t[i] << " ";
	}
	cout << endl;

}

void TreeGloverModPR::print_values_vertical(vector<int> p, vector<int> t) const {

	cout << "\tp\tt" << endl;
	for(int i=0; i<(int)p.size(); i++) {
		cout << i << ":\t" << p[i] << "\t" << t[i] << endl;
	}
}

Edge TreeGloverModPR::performEdgeExchange(Edge ie, vector<vector<int> > &sol, int c) {

	Edge oe;

	incoming_edge = ie;		// Atribui aresta que vai entrar na árvore

	u = ie.endpoint1;
	v = ie.endpoint2;

	// Encontra o nó predecessor em comum das extremidades da aresta que vai
	// entrar na solução e atribui à outgoing_edge a aresta de maior custo que sairá da solução
	oe = findBasisEquivalentPath(sol, c);

	setEdges(ie, oe);
	//	cout << ie.endpoint1 << " " << ie.endpoint2 << " - " << oe.endpoint1 << " " << oe.endpoint2 << endl;
	//	cin.get();
	// ---------------------------------------------------------------//

	double cost = Instance::getInstance()->cost_matrix[u][v] - Instance::getInstance()->cost_matrix[p][q];		// Custo reduzido da aresta à entrar na solução

	tree_cost += cost;

	updatingOperations();		// Atualiza os valores de 'p', 't' e 'a' com a remoção da aresta outgoing_edge

	decidingT2();		// Decide qual das duas árvores deve ser a árvore T1 e a árvore T2

	reRootingT2();		// Atualiza a raiz da árvore T2

	attachT2T1();		// Une as árvores T1 e T2

	return oe;

}

void TreeGloverModPR::setEdges(Edge ie, Edge oe) {

	incoming_edge = ie;
	outgoing_edge = oe;

	if(tx[incoming_edge.endpoint1] > tx[incoming_edge.endpoint2]) {
		u = incoming_edge.endpoint1;
		v = incoming_edge.endpoint2;
	}
	else {
		u = incoming_edge.endpoint2;
		v = incoming_edge.endpoint1;
	}

	if(tx[outgoing_edge.endpoint1] > tx[outgoing_edge.endpoint2]) {
		p = outgoing_edge.endpoint1;
		q = outgoing_edge.endpoint2;
	}
	else {
		p = outgoing_edge.endpoint2;
		q = outgoing_edge.endpoint1;
	}

}

Edge TreeGloverModPR::findBasisEquivalentPath(vector<vector<int> > &sol, int op) {
	//	cout << "op = " << op << endl;
	int x_u = u;
	int x_v = v;
	bool step3, step4, step5, step6;
	Edge most_expensive_edge;
	Edge oe;
	bool first_time;

	first_time = true;		// flag necessária pra pegar a aresta incidente à endpoint1 ou endpoint2 (depende de op)
	// no ciclo formado com a inclusão da aresta e

	step3 = step4 = step5 = step6 = true;

	if(tx[x_u] < tx[x_v]) {
		step3 = false;
		step4 = false;
		step5 = true;
		step6 = false;
	}

	while(true) {

		if(step3) {
			if(tx[x_u] > tx[x_v]) {
				step3 = false;
				step4 = false;
				step5 = false;
				step6 = true;
			}
			else {
				step3 = false;
				step4 = true;
				step5 = false;
				step6 = false;
			}
		}

		if(step4) {
			if(x_u != x_v) {
				step3 = false;
				step4 = false;
				step5 = true;
				step6 = false;
			}
			else {
				z = x_u;
				if(op == 0)
					return most_expensive_edge;
				else
					return oe;
			}
		}

		if(step5) {
			if(op != 0) {
				//				cout << "{" << x_u << ", " << px[x_u] << "} -> " << sol[x_u][px[x_u]] << endl;
				if(px[x_u] != -1 && op == 1 && sol[x_u][px[x_u]] != 1 && first_time == true) {			// x_u não é o nó raiz, endpoint1 está saturado e é a primeira vez que o algoritmo passa por step5?
					oe.setEdge(x_u, px[x_u], Instance::getInstance()->cost_matrix[x_u][px[x_u]]);		// a aresta que deve sair da árvore é {x_u, p[x_u]} pois é a aresta que conecta endpoint1 ao ciclo
					first_time = false;
				}
			}
			else {
				//				cout << "{" << x_u << ", " << px[x_u] << "} -> " << sol[x_u][px[x_u]] << endl;
				if(px[x_u] != -1 && sol[px[x_u]][x_u] != 1)
					if(Instance::getInstance()->cost_matrix[px[x_u]][x_u] > most_expensive_edge.cost)
						most_expensive_edge.setEdge(px[x_u], x_u, Instance::getInstance()->cost_matrix[px[x_u]][x_u]);
			}

			int x_aux = x_u;	// guarda o vértice posterior à x no ciclo formado pela inclusão de e
			int x = px[x_u];
			cout << "u: " << u << " " << px[x_u] << endl;
			while(tx[x] < tx[x_v]) {
				if(px[x] != -1) {
					cout << "u: " << px[x] << " " << x << endl;
					if(Instance::getInstance()->cost_matrix[px[x]][x] > most_expensive_edge.cost)
						most_expensive_edge.setEdge(px[x], x, Instance::getInstance()->cost_matrix[px[x]][x]);
				}

				x_aux = x;
				x = px[x];
			}

			if(op == 2 && first_time == true && sol[x_aux][x] != 1 && x == v) {			// O ponto de equivalencia z é a extremidade endpoint2 e este é o nó saturado?
				oe.setEdge(x_aux, x, Instance::getInstance()->cost_matrix[x_aux][x]);	// a aresta que deve sair da árvore é {x_aux, x} = {endpoint2, x_aux}
				first_time = false;
			}

			x_u = x;

			step3 = true;
			step4 = false;
			step5 = false;
			step6 = false;
		}

		if(step6) {
			if(op != 0) {
				//				cout << "{" << x_v << ", " << px[x_v] << "} -> " << sol[x_v][px[x_v]] << endl;
				if(px[x_v] != -1 && op == 2 && sol[x_v][px[x_v]] != 1 && first_time == true) {			// x_v não é o nó raiz, endpoint2 está saturado e é a primeira vez que o algoritmo passa por step6?
					oe.setEdge(x_v, px[x_v], Instance::getInstance()->cost_matrix[x_v][px[x_v]]);	// a aresta que deve sair da árvore é {x_v, p[x_v]} pois é a aresta que conecta endpoint2 ao ciclo
					first_time = false;
				}
			}
			else {
				//				cout << "{" << x_v << ", " << px[x_v] << "} -> " << sol[x_v][px[x_v]] << endl;
				if(px[x_v] != -1 && sol[px[x_v]][x_v] != 1)
					if(Instance::getInstance()->cost_matrix[px[x_v]][x_v] > most_expensive_edge.cost)
						most_expensive_edge.setEdge(px[x_v], x_v, Instance::getInstance()->cost_matrix[px[x_v]][x_v]);
			}

			int x_aux = x_v;		// guarda o vértice posterior à x no ciclo formado pela inclusão de e
			int x = px[x_v];
			cout << "v: " << v << " " << px[x_v] << endl;
			while(tx[x] < tx[x_u]) {
				if(px[x] != -1) {
					cout << "v: " << px[x] << " " << x << endl;
					if(Instance::getInstance()->cost_matrix[px[x]][x] > most_expensive_edge.cost)
						most_expensive_edge.setEdge(px[x], x, Instance::getInstance()->cost_matrix[px[x]][x]);
				}

				x_aux = x;
				x = px[x];
			}

			if(op == 1 && first_time == true && sol[x_aux][x] != 1 && x == u) {			// O ponto de equivalencia z é a extremidade endpoint1 e este é o nó saturado?
				oe.setEdge(x_aux, x, Instance::getInstance()->cost_matrix[x_aux][x]);		// a aresta que deve sair da árvore é {x_aux, x} = {endpoint1, x_aux}
				first_time = false;
			}

			x_v = x;

			if(tx[x_u] == tx[x_v]) {
				step3 = false;
				step4 = true;
				step5 = false;
				step6 = false;
			}
			else {
				step3 = false;
				step4 = false;
				step5 = true;
				step6 = false;
			}
		}
	}

	Edge e;
	return e;
}

void TreeGloverModPR::updatingOperations() {

	int x_p = p;
	int x_q = q;


	vector<int> p_ast = px;
	vector<int> t_ast = tx;

	/*---------- Update for T - T(q) ----------*/

	//update t(x)
	int x = x_p;
	while(x != root) {
		t_ast[x] = tx[x] - tx[x_q];
		x = px[x];
	}
	t_ast[root] = tx[root] - tx[x_q];

	/*--------------------------------------------*/

	/*------------- Update for T(q) --------------*/
	//update p(x)
	p_ast[x_q] = -1;
	/*--------------------------------------------*/

	px = p_ast;
	tx = t_ast;

	// Atualiza matriz de adjacência
	adj_matrix[p][q] = false;
	adj_matrix[q][p] = false;
	adj_matrix[u][v] = true;
	adj_matrix[v][u] = true;

	// Atualiza os valores de 'a' devido a remoção e a inclusão das arestas
	for(int i=0; i<num_vertex; i++) {
		ax[u] = 0;
		ax[v] = 0;
		ax[p] = 0;
		ax[q] = 0;
		if(adj_matrix[p][i] == true && Instance::getInstance()->cost_matrix[p][i] > ax[p]) {
			ax[p] = Instance::getInstance()->cost_matrix[p][i];
		}
		if(adj_matrix[q][i] == true && Instance::getInstance()->cost_matrix[q][i] > ax[q]) {
			ax[q] = Instance::getInstance()->cost_matrix[q][i];
		}
		if(adj_matrix[u][i] == true && Instance::getInstance()->cost_matrix[u][i] > ax[u]) {
			ax[u] = Instance::getInstance()->cost_matrix[u][i];
		}
		if(adj_matrix[v][i] == true && Instance::getInstance()->cost_matrix[v][i] > ax[v]) {
			ax[v] = Instance::getInstance()->cost_matrix[v][i];
		}
	}

	//	print_values_vertical(p, t);
	//	cout << endl << endl;
	//	cin.get();
}

void TreeGloverModPR::decidingT2() {

	if(tx[root] > tx[q]) {
		x1 = root;
		x2 = q;
	}
	else {
		x1 = q;
		x2 = root;
		root = q;			// observar aqui
	}

	int y = u;
	while(true) {
		if(y == x1) {
			y1 = u;
			y2 = v;
			break;
		}
		else if(y == x2) {
			y2 = u;
			y1 = v;
			break;
		}
		y = px[y];
	}

}

void TreeGloverModPR::reRootingT2() {

	vector<int> p_ast = px;
	vector<int> t_ast = tx;
	int x, y;
	bool nothing_to_do = false;

	// p*(x)
	if(x2 == y2)
		nothing_to_do = true;
	else {
		y = y2;
	}

	if(!nothing_to_do) {
		while(true) {
			p_ast[px[y]] = y;
			if(px[y] == x2) {
				break;
			}
			y = px[y];
		}
		p_ast[y2] = -1;
	}
	//------

	// t*(x)
	t_ast[y2] = tx[x2];

	x = y2;
	while(x != x2) {
		t_ast[px[x]] = t_ast[y2] - tx[x];
		x = px[x];
	}

	px = p_ast;
	tx = t_ast;

	//	print_values_vertical(p, t);
	//	cin.get();

}

void TreeGloverModPR::attachT2T1() {

	vector<int> p_ast = px;
	vector<int> t_ast = tx;

	p_ast[y2] = y1;

	int x = y1;
	while(true) {
		t_ast[x] = tx[x] + tx[y2];
		if(x == x1)
			break;
		x = px[x];
	}

	px = p_ast;
	tx = t_ast;

	//	print_values_vertical(p, t);
	//	cin.get();
}

vector<Edge> TreeGloverModPR::getTree() {

	vector<Edge> tree;
	for(int i=0; i<num_vertex; i++) {
		if(px[i] != -1) {
			int edp1, edp2;
			if(i < px[i]) {
				edp1 = i;
				edp2 = px[i];
			}
			else {
				edp1 = px[i];
				edp2 = i;
			}
			Edge edge(px[i], i, Instance::getInstance()->cost_matrix[edp1][edp2]);
			tree.push_back(edge);
		}
	}

	return tree;

}

double  TreeGloverModPR::getTreeCost() {
	return tree_cost;
}

PathRelinking::PathRelinking(vector<Edge> &s, vector<Edge> &t) {

	double cost_start = 0;
	start = s;
	target = t;
	sol_start.resize(Instance::getInstance()->num_vertex);
	sol_target.resize(Instance::getInstance()->num_vertex);
	degree.resize(Instance::getInstance()->num_vertex, 0);
	sol.resize(Instance::getInstance()->num_vertex);
	for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
		sol_start[i].resize(Instance::getInstance()->num_vertex, 0);
		sol_target[i].resize(Instance::getInstance()->num_vertex, 0);
		sol[i].resize(Instance::getInstance()->num_vertex, 0);
	}

	for(int i=0; i < Instance::getInstance()->num_vertex - 1; i++) {
		sol_start[start[i].endpoint1][start[i].endpoint2] = sol_start[start[i].endpoint2][start[i].endpoint1] = 1;
		sol_target[target[i].endpoint1][target[i].endpoint2] = sol_target[target[i].endpoint2][target[i].endpoint1] = 1;
		cost_start += start[i].cost;
		++degree[start[i].endpoint1];
		++degree[start[i].endpoint2];
	}

	for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
		for(int j=0; j < Instance::getInstance()->num_vertex; j++) {
			if(sol_start[i][j] == sol_target[i][j]) {
				sol[i][j] = sol_target[i][j];
			}
			else
				sol[i][j] = 0;
		}
	}

	tree = new TreeGloverModPR(start, cost_start);

}

PathRelinking::~PathRelinking() {
	delete tree;
}

vector<Edge> PathRelinking::runPR() {

	queue<Edge> Q;

	sort(target.begin(), target.end());

	vector<int> degrees(Instance::getInstance()->num_vertex, 0);
	for(int i=0; i<Instance::getInstance()->num_vertex; i++)
		cout << i << " -> " << Instance::getInstance()->degree_constraint[i] << endl;

	for(unsigned int i=0; i < start.size(); i++) {
		// Calcula os graus dos vértices da solução de início do PR
		++degrees[start[i].endpoint1];
		++degrees[start[i].endpoint2];

		// Coleta as arestas que estão na solução alvo mas não estão na solução de início
		if(sol[target[i].endpoint1][target[i].endpoint2] != 1)
			Q.push(target[i]);
	}

	vector<Edge> best_tree = start;
	double best_cost = INF;

	int cont = 0;
	while(Q.size() > 0) {
		cont++;

		Edge edge = Q.front();
		Q.pop();
		int edp1 = edge.endpoint1;
		int edp2 = edge.endpoint2;

		int op = 3;

		if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2])
			op = 0;
		else if(degrees[edp1] <= Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2])
			op = 1;
		else if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] <= Instance::getInstance()->degree_constraint[edp2])
			op = 2;
//		else {
//			cout << "{" << edp1 << ", " << edp2 << "}" << endl;
//			cin.get();
//			Q.push(edge);
//			tree->performEdgeExchange(edge, sol, op);
//			continue;
//		}
		sol[edp1][edp2] = sol[edp2][edp1] = 1;
		Edge oe = tree->performEdgeExchange(edge, sol, op);
		++degrees[edp1];
        ++degrees[edp2];
        --degrees[oe.endpoint1];
        --degrees[oe.endpoint2];
		if(tree->getTreeCost() < best_cost) {
		    bool valid_tree = true;
		    for(int i=0; i < Instance::getInstance()->num_vertex; i++)
		        if(degrees[i] > Instance::getInstance()->degree_constraint[i]) {
		            valid_tree = false;
		            continue;
		        }
		    if(valid_tree == true) {
                best_cost = tree->getTreeCost();
                best_tree = tree->getTree();
		    }
		}
	}

	return best_tree;
}
