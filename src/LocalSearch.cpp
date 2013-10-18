/*
 * LocalSearch.cpp
 *
 *  Created on: Oct 10, 2013
 *      Author: luishbicalho
 */

#include "LocalSearch.h"

LocalSearch::LocalSearch(const vector<Edge> &s, double c) {

    sol = s;
    obj_func_value = c;
    tree_sol = new TreeGlover(sol, obj_func_value);
    degrees.resize(Instance::getInstance()->num_vertex, 0);

    for(unsigned int i=0; i < s.size(); i++) {
        ++degrees[s[i].endpoint1];
        ++degrees[s[i].endpoint2];
    }

}

LocalSearch::~LocalSearch() {
    if(tree_sol)
        delete tree_sol;
}

void LocalSearch::runLs(vector<Edge> &edges) {

    for(uint i=0; i < edges.size(); i++) {
        int edp1 = edges[i].endpoint1;
        int edp2 = edges[i].endpoint2;

        // Ambas as extremidades estão folgadas?
        if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2]) {
            Edge outgoing = tree_sol->performEdgeExchange(edges[i]);
            if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {     // A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
                // Calcula os graus decorrente da troca de arestas
                --degrees[outgoing.endpoint1];
                --degrees[outgoing.endpoint2];
                ++degrees[edges[i].endpoint1];
                ++degrees[edges[i].endpoint2];
                //------------------------------------------------
                i = 0;                  // Inspeciona as arestas do começo (first improvement approach)
            }
        }
        // A extremidade 1 está apertada (x(\delta(edp1)) == d_{edp1}) e a extremidade 2 está folgada?
        else if(degrees[edp1] <= Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2]) {
            Edge outgoing = tree_sol->performEdgeExchange(edges[i], 1);
            if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {     // A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
                // Calcula os graus decorrente da troca de arestas
                --degrees[outgoing.endpoint1];
                --degrees[outgoing.endpoint2];
                ++degrees[edges[i].endpoint1];
                ++degrees[edges[i].endpoint2];
                //------------------------------------------------
                i = 0;                  // Inspeciona as arestas do começo (first improvement approach)
            }
        }
        // A extremidade 1 está folgada e a extremidade 2 está apertada (x(\delta(edp2)) == d_{edp2})?
        else if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] <= Instance::getInstance()->degree_constraint[edp2]) {
            Edge outgoing = tree_sol->performEdgeExchange(edges[i], 2);
            if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {     // A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
                // Calcula os graus decorrente da troca de arestas
                --degrees[outgoing.endpoint1];
                --degrees[outgoing.endpoint2];
                ++degrees[edges[i].endpoint1];
                ++degrees[edges[i].endpoint2];
                //------------------------------------------------
                i = 0;                  // Inspeciona as arestas do começo (first improvement approach)
            }
        }
    }

    sol = tree_sol->getTree();
    obj_func_value = tree_sol->getTreeCost();

}


