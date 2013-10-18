/*
 * VariableFixing.cpp
 *
 *  Created on: 19/09/2013
 *      Author: luishbicalho
 */

#include "VariableFixing.h"

VariableFixing::VariableFixing(Graph *G, vector<vector<double> > &lagr_cost_matrix, double c) {

    root = 0;
    num_vertex = Instance::getInstance()->num_vertex;



    px.resize(num_vertex);
    tx.resize(num_vertex, 1);

    adj_list.resize(num_vertex);
    for(int i=0; i<num_vertex; i++)
        adj_list[i].resize(num_vertex, 0);

    vector<KEdge> lagr_mst = kruskal(G);
    lb = 0;
    for(unsigned int i=0; i < lagr_mst.size(); i++) {
        adj_list[lagr_mst[i].p][lagr_mst[i].q] = true;
        adj_list[lagr_mst[i].q][lagr_mst[i].p] = true;
        lb += lagr_mst[i].w;
    }
    lb += c;

    initializeTree(adj_list);

    cost_matrix = lagr_cost_matrix;

}

void VariableFixing::dfs(stack<int> &S, vector<int> &discovered, vector<vector<bool> > &adj_list) {

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

void VariableFixing::initializeTree(vector<vector<bool> > &adj_list) {

    stack<int> S;
    vector<int> discovered(num_vertex, 0);
    root = 0;
    px[root] = -1;
    discovered[root] = 1;
    S.push(root);
    dfs(S, discovered, adj_list);

}

Edge VariableFixing::findLargestCostEdge(int u, int v) {

    int x_u = u;
    int x_v = v;
    bool step3, step4, step5, step6;
    Edge most_expensive_edge;

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
                return most_expensive_edge;
            }
        }

        if(step5) {
            if(px[x_u] != -1)
                if(cost_matrix[px[x_u]][x_u] > most_expensive_edge.cost)
                    most_expensive_edge.setEdge(px[x_u], x_u, cost_matrix[px[x_u]][x_u]);

            int x = px[x_u];
            while(tx[x] < tx[x_v]) {
                if(px[x] != -1)
                    if(cost_matrix[px[x]][x] > most_expensive_edge.cost)
                        most_expensive_edge.setEdge(px[x], x, cost_matrix[px[x]][x]);
                x = px[x];
            }
            x_u = x;

            step3 = true;
            step4 = false;
            step5 = false;
            step6 = false;
        }

        if(step6) {
            if(px[x_v] != -1)
                if(cost_matrix[px[x_v]][x_v] > most_expensive_edge.cost)
                    most_expensive_edge.setEdge(px[x_v], x_v, cost_matrix[px[x_v]][x_v]);

            int x = px[x_v];
            while(tx[x] < tx[x_u]) {
                if(px[x] != -1)
                    if(cost_matrix[px[x]][x] > most_expensive_edge.cost)
                        most_expensive_edge.setEdge(px[x], x, cost_matrix[px[x]][x]);
                x = px[x];
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

int VariableFixing::runVariableFixing(double ub, vector<vector<int> > &edges_in_pool, vector<Edge> &all_edges, vector<Edge*> &columns_) {

    int cont = 0;

    for(unsigned int i = 0; i < all_edges.size(); i++) {
        if(adj_list[all_edges[i].endpoint1][all_edges[i].endpoint2] == false && all_edges[i].fixed == false) {
            Edge largest_cost_edge = findLargestCostEdge(all_edges[i].endpoint1, all_edges[i].endpoint2);
            double cr = cost_matrix[all_edges[i].endpoint1][all_edges[i].endpoint2] - largest_cost_edge.cost;

            if(Utils::cmp(cr, ub - lb) > 0) {
                all_edges[i].fixed = true;
                if(edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2] != -1) {
                    columns_[edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2]]->fixed = true;
                }
                cont++;
            }
        }
    }

    return cont;
}

