/*
 * MinCut.cpp
 *
 *  Created on: May 21, 2013
 *      Author: luishbicalho
 */


#include "MinCut.h"


MinCut::MinCut() {

	num_vertex = Instance::getInstance()->num_vertex+2;
	source = Instance::getInstance()->num_vertex;
	sink = Instance::getInstance()->num_vertex+1;

};

MinCut::~MinCut() {

}

int MinCut::solveMinCutWolsey(IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector<vector<int> > &min_cuts) {


    int id_most_violated = -1;
    double most_violated_value = 0;

    Dinic mf_base(num_vertex);

    vector<double> b(Instance::getInstance()->num_vertex, 0);
    for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
        for(int j=i+1; j<Instance::getInstance()->num_vertex; j++) {
            if(edges_in_pool[i][j] != -1) {
                int pos = edges_in_pool[i][j];
                if(xsol[pos] > 0.00001) {       // A aresta i possui valor de relax. linear maior que zero?
                    mf_base.addEdge(i, j, (xsol[pos]/2.0)*100000.0);        // Adiciona arco de ida
                    mf_base.addEdge(j, i, (xsol[pos]/2.0)*100000.0);        // Adiciona arco de volta
                    b[i] += xsol[pos];          // Incrementa a contribuição do arco i ao grau do vértice extremidade 1
                    b[j] += xsol[pos];          // Incrementa a contribuição do arco i ao grau do vértice extremidade 2
                }
            }
        }
    }

    // Adiciona as arestas referentes ao vértice fonte e ao vértice sorvedouro
    // Vértice fonte = n
    // Vértice sorvedouro = n+1

    for(int v=0; v < Instance::getInstance()->num_vertex; v++) {
        if(b[v]/2.0 - 1.0 > 0.00001) {
            mf_base.addEdge(source, v, (b[v]/2.0 - 1.0)*100000.0);
        }
        else {
            mf_base.addEdge(source, v, 0);
        }

        if(1.0 - b[v]/2.0 > 0.00001) {
            mf_base.addEdge(v, sink, (1.0 - b[v]/2.0)*100000.0);
        }
        else {
            mf_base.addEdge(v, sink, 0);
        }
    }

    for(int k=0; k<Instance::getInstance()->num_vertex-2; k++) {

        mf_base.g_adj_list[source][k].cap = INF;
        Dinic mf = mf_base;

        mf.maxFlow(source, sink);

        vector<int> cut(Instance::getInstance()->num_vertex, 0);        // 1 se o vértice pertence ao corte, 0 caso contrário.
        for(int i=0; i<Instance::getInstance()->num_vertex; i++)
            if(mf.level[i] != -1)
                cut[i] = 1;

        double sum = 0;
        int cardS = 0;
        for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
            if(cut[i] == 1) {   // O vertice i pertence a S?
                cardS++;        // Adiciona a unidade referente ao elemento i pertencente à S
                for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
                    if(cut[j] == 1 && edges_in_pool[i][j] != -1) {      // Ambos os vertices, i e j, pertencem à S e A aresta está presente no pool de colunas?
                        sum += xsol[edges_in_pool[i][j]];       // Soma o valor de relaxacao linear referente à aresta {i, j}
                    }
                }
            }
        }

        if(sum - double(cardS-1) > 0.000001 && cardS > 1 && cardS < Instance::getInstance()->num_vertex) {      // O corte atual está violado e S atende as restrições de cardinalidade?
            min_cuts.push_back(cut);
            if(sum - double(cardS-1) > most_violated_value) {
                id_most_violated = (int)min_cuts.size() - 1;
                most_violated_value = sum - double(cardS-1);
            }
        }

        if(b[k]/2.0 - 1.0 > 0.00001 ) {
            mf_base.g_adj_list[source][k].cap = (int)((b[k]/2.0 - 1.0)*100000.0);
        }
        else {
            mf_base.g_adj_list[source][k].cap = 0;
        }
        unsigned int last_pos = mf_base.g_adj_list[k].size() - 1;
        mf_base.g_adj_list[k][last_pos].cap = 100000000;

        cut.clear();
    }

    return id_most_violated;

}

int MinCut::solveMinCutWolseyOrth(IloNumArray &xsol, vector<vector<int> > &edges_in_pool, vector<int> &vec, vector<vector<int> > &min_cuts) {


    int id_most_violated = -1;
    double most_violated_value = 0;

    Dinic mf_base(num_vertex);

    vector<double> b(Instance::getInstance()->num_vertex, 0);
    for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
        for(int j=i+1; j<Instance::getInstance()->num_vertex; j++) {
            if(edges_in_pool[i][j] != -1 && (vec[i] != 1 || vec[j] != 1)) {   // Se a aresta {i,j} está no pool mas não está no corte vec...
                int pos = edges_in_pool[i][j];
                if(xsol[pos] > 0.00001) {       // A aresta i possui valor de relax. linear maior que zero?
                    mf_base.addEdge(i, j, (xsol[pos]/2.0)*100000.0);        // Adiciona arco de ida
                    mf_base.addEdge(j, i, (xsol[pos]/2.0)*100000.0);        // Adiciona arco de volta
                    b[i] += xsol[pos];          // Incrementa a contribuição do arco i ao grau do vértice extremidade 1
                    b[j] += xsol[pos];          // Incrementa a contribuição do arco i ao grau do vértice extremidade 2
                }
            }
        }
    }

    // Adiciona as arestas referentes ao vértice fonte e ao vértice sorvedouro
    // Vértice fonte = n
    // Vértice sorvedouro = n+1

    for(int v=0; v < Instance::getInstance()->num_vertex; v++) {
        if(b[v]/2.0 - 1.0 > 0.00001) {
            mf_base.addEdge(source, v, (b[v]/2.0 - 1.0)*100000.0);
        }
        else {
            mf_base.addEdge(source, v, 0);
        }

        if(1.0 - b[v]/2.0 > 0.00001) {
            mf_base.addEdge(v, sink, (1.0 - b[v]/2.0)*100000.0);
        }
        else {
            mf_base.addEdge(v, sink, 0);
        }
    }

    for(int k=0; k<Instance::getInstance()->num_vertex-2; k++) {

        mf_base.g_adj_list[source][k].cap = INF;
        Dinic mf = mf_base;

        mf.maxFlow(source, sink);

        vector<int> cut(Instance::getInstance()->num_vertex, 0);        // 1 se o vértice pertence ao corte, 0 caso contrário.
        for(int i=0; i<Instance::getInstance()->num_vertex; i++)
            if(mf.level[i] != -1)
                cut[i] = 1;

        double sum = 0;
        int cardS = 0;
        for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
            if(cut[i] == 1) {   // O vertice i pertence a S?
                cardS++;        // Adiciona a unidade referente ao elemento i pertencente à S
                for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
                    if(cut[j] == 1 && edges_in_pool[i][j] != -1) {      // Ambos os vertices, i e j, pertencem à S e A aresta está presente no pool de colunas?
                        sum += xsol[edges_in_pool[i][j]];       // Soma o valor de relaxacao linear referente à aresta {i, j}
                    }
                }
            }
        }

        if(sum - double(cardS-1) > 0.000001 && cardS > 1 && cardS < Instance::getInstance()->num_vertex) {      // O corte atual está violado e S atende as restrições de cardinalidade?
            min_cuts.push_back(cut);
            if(sum - double(cardS-1) > most_violated_value) {
                id_most_violated = (int)min_cuts.size() - 1;
                most_violated_value = sum - double(cardS-1);
            }
        }

        if(b[k]/2.0 - 1.0 > 0.00001 ) {
            mf_base.g_adj_list[source][k].cap = (int)((b[k]/2.0 - 1.0)*100000.0);
        }
        else {
            mf_base.g_adj_list[source][k].cap = 0;
        }
        unsigned int last_pos = mf_base.g_adj_list[k].size() - 1;
        mf_base.g_adj_list[k][last_pos].cap = 100000000;

        cut.clear();
    }

    return id_most_violated;
}

