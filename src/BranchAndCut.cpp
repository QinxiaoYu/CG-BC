/*
 * BranchAndCut.cpp
 *
 *  Created on: 01/10/2013
 *      Author: luishbicalho
 */

#include "BranchAndCut.h"

ILOLAZYCONSTRAINTCALLBACK3(bcMinCutCallback,
        const IloNumVarArray&,  bc_x_,
        vector<Edge>&, columns_,
        vector<vector<int> >&, edges_in_pool) {


    try {
        IloNumArray xsol(getEnv(), columns_.size());            // Valores de relaxação linear das variáveis x_
        getValues(xsol, bc_x_);
        vector<vector<int> > subsets;

        // A separação exada das SEC's
        MinCut mc;
        mc.solveMinCutWolsey(xsol, edges_in_pool, subsets);       // Subconjuntos S de cortes violados

        for(int k=0; k<(int)subsets.size(); k++) {
            IloExpr expr(getEnv());
            int sum_set = 0;        // Calcula o número de elementos no conjunto S da sec k
            for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
                if(subsets[k][i] == 1)
                    sum_set++;
                else
                    continue;
                for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
                    if(subsets[k][i] == 1 && subsets[k][j] == 1 && edges_in_pool[i][j] != -1)     //Ambas as extremidades i e j pertencem à S e a aresta {i,j} está no modelo?
                        expr += bc_x_[edges_in_pool[i][j]];       // Adiciona a variável de decisão referente a aresta {i,j} na nova restrição
                }
            }

            add(expr <= sum_set - 1);

            expr.end();
        }

        xsol.end();
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }
}

ILOUSERCUTCALLBACK3(bcBlossomCallback,
        const IloNumVarArray&,  bc_x_,
        vector<Edge>&, columns_,
        vector<vector<int> >&, edges_in_pool) {


    try {
        IloNumArray xsol(getEnv(), bc_x_.getSize());
        getValues(xsol, bc_x_);

        vector<Edge> support_graph;         // Grafo suporte da solução. O custo da aresta é o valor de relaxação linear desta no mestre
        vector<vector<int> > adj_list(Instance::getInstance()->num_vertex);         // Lista de adjacência do grafo resultante da retirado das arestas com x_e = 1 do grafo suporte

        vector<double> s_line(Instance::getInstance()->num_vertex, 0.0);        // folga do vértice i (b|i| - x(delta(i)))
        vector<double> b_line(Instance::getInstance()->num_vertex);             // grau do vértice i  (b|i| - )

        for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
            b_line[i] = Instance::getInstance()->degree_constraint[i];
            s_line[i] = Instance::getInstance()->degree_constraint[i];
        }

                    // Valores de relaxação linear das variáveis x_
        for(unsigned int i = 0; i < columns_.size(); i++) {
            if(Utils::cmp(xsol[i], 0) > 0) {                // A aresta referente à coluna i pertence à solução?
                Edge edge(columns_[i].endpoint1, columns_[i].endpoint2, xsol[i]);

                support_graph.push_back(edge);          // Adiciona a aresta ao grafo suporte

                s_line[edge.endpoint1] -= xsol[i];
                s_line[edge.endpoint2] -= xsol[i];


                if(Utils::cmp(xsol[i], 1) != 0) {       //O valor de relaxação linear é diferente de 1?
                    adj_list[edge.endpoint1].push_back(edge.endpoint2);
                    adj_list[edge.endpoint2].push_back(edge.endpoint1);
                }
                else {
                    b_line[edge.endpoint1] -= 1;
                    b_line[edge.endpoint2] -= 1;
                }
            }
        }

        vector<vector<int> > H;     // Conjuntos Handle - Cada conjunto handle é um vector no num_vertex posições onde H[i][j] = 1 indica que o vértice j
                                    // pertence ao handle i, e H[i][j] = 0 caso contrário
        vector<vector<Edge> > T;    // Conjuntos Tooth

        Utils::findHandles(adj_list, H);

        T.resize(H.size());
        vector<int> violated_blossom;       // Conjunto de indices de BI's violadas em H

        for(unsigned int i=0; i < H.size(); i++) {
            double sum_slack = 0;
            double sum_degree = 0;
            for(int j=0; j < Instance::getInstance()->num_vertex; j++) {
                if(H[i][j] == 1) {          // O vértice j pertence ao conjunto Handle i?
                    sum_slack += s_line[j];
                    sum_degree += b_line[j];
                }
            }

            if(Utils::cmp(sum_slack, 1) < 0 && ((int)sum_degree % 2) != 0) {            // A blossom está violada?
                violated_blossom.push_back(i);
                for(unsigned int j = 0; j < support_graph.size(); j++) {
                    bool test1 = (H[i][support_graph[j].endpoint1] == 1 && H[i][support_graph[j].endpoint2] != 1);  // endpoint1 pertence à H mas endpoint2 não?
                    bool test2 = (H[i][support_graph[j].endpoint1] != 1 && H[i][support_graph[j].endpoint2] == 1);  // endpoint1 não pertence à H mas endpoint2 pertence?
                    if(test1 || test2) {    // Se somente uma das extremidades pertence à H...
                        T[i].push_back(support_graph[j]);
                    }
                }
            }
        }

        for(uint p = 0; p < violated_blossom.size(); p++) {
            int k = violated_blossom[p];
            IloExpr expr(getEnv());
            int sum_degrees;
            sum_degrees = 0;        // Soma dos graus dos vértices presentes em H

            for(unsigned int i=0; i < H[k].size(); i++) {
                if(H[k][i] == 1) {     // O vértice i pertence à H?
                    sum_degrees += Instance::getInstance()->degree_constraint[i];
                    for(unsigned int j=i+1; j < H[k].size(); j++) {
                        if(H[k][j] == 1)           // O vértice j pertence à H?
                            if(edges_in_pool[i][j] != -1)       // A aresta {i,j} está no modelo?
                                expr += bc_x_[edges_in_pool[i][j]];        // Inclui a aresta {i,j} na restrição que representa a blossom
                    }
                }
            }

            for(unsigned int i=0; i < T[k].size(); i++) {
                expr += bc_x_[edges_in_pool[T[k][i].endpoint1][T[k][i].endpoint2]];          // Adiciona a aresta presente em T (já está no modelo) à restrição
            }

            if(sum_degrees + T[k].size() % 2 != 0)
                add(expr <= floor(0.5*(sum_degrees + T[k].size())));

            expr.end();
        }
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }
}

ILOHEURISTICCALLBACK3(kruskalxCallBack,
        IloNumVarArray, bc_x_,
        vector<Edge>& , columns,
        vector<vector<int> >&, edges_in_pool) {

    try {
        IloNumArray xsol(getEnv());
        getValues(xsol, bc_x_);

        IloNumArray sol_heur(getEnv(), bc_x_.getSize());
        for(uint i =0; i < bc_x_.getSize(); i++)
            sol_heur[i] = 0.0;
        IloNum obj = 0.0;

        Graph *G_bar = new Graph(Instance::getInstance()->num_edges);
        for(uint i=0; i < columns.size(); ++i) {
            if(Utils::cmp(xsol[ edges_in_pool[columns[i].endpoint1][columns[i].endpoint2] ], 0) > 0)
                G_bar->push_edge(columns[i].endpoint1, columns[i].endpoint2,
                    columns[i].cost*(1 - xsol[edges_in_pool[columns[i].endpoint1][columns[i].endpoint2]]));   // O custo é substituído pelo custo complementar através do valor de RL
            else
                G_bar->push_edge(columns[i].endpoint1, columns[i].endpoint2, columns[i].cost);
        }

        vector<KEdge> st_heur = kruskalDegreeConstrained(G_bar);
        vector<Edge> sol;

        if((int)st_heur.size() == Instance::getInstance()->num_vertex - 1) {
            for(uint i=0; i < st_heur.size(); ++i) {
                int edp1 = st_heur[i].p;
                int edp2 = st_heur[i].q;
                obj += Instance::getInstance()->cost_matrix[edp1][edp2];
                Edge e(edp1, edp2, Instance::getInstance()->cost_matrix[edp1][edp2]);
                sol.push_back(e);
            }

            LocalSearch ls(sol, obj);
            ls.runLs(columns);

            sol = ls.tree_sol->getTree();
			obj = ls.tree_sol->getTreeCost();

			for(uint i=0; i < sol.size(); i++) {
            	sol_heur[edges_in_pool[sol[i].endpoint1][sol[i].endpoint2]] = 1.0;
            }

            setSolution(bc_x_, sol_heur, obj);
        }

        xsol.end();
        sol_heur.end();
        delete G_bar;

    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }
}

BranchAndCut::BranchAndCut(ColumnGeneration &cg) : bc_model_(bc_env_) {


    int num_vertex = Instance::getInstance()->num_vertex;
    bc_cplex_ = new IloCplex(bc_model_);
    edges_in_pool.resize(num_vertex);

    for(int i=0; i < num_vertex; i++) {
        edges_in_pool[i].resize(num_vertex, -1);
    }

    for(uint i=0; i < cg.all_edges.size(); i++) {
        if(cg.all_edges[i].fixed == false) {
            columns_.push_back(cg.all_edges[i]);
            edges_in_pool[cg.all_edges[i].endpoint1][cg.all_edges[i].endpoint2] = columns_.size()-1;
            edges_in_pool[cg.all_edges[i].endpoint2][cg.all_edges[i].endpoint1] = columns_.size()-1;
        }
    }

    createVariables();

    createObjectiveFunction();

    createConstraints(cg);

    bc_cplex_->setParam(IloCplex::CutUp, cg.best_primal.obj_fun_value);

}

void BranchAndCut::createVariables() {

    try{
        bc_x_ = IloNumVarArray(bc_env_, columns_.size(), 0.0, 1.0, ILOBOOL);

        // renomeia variaveis x
        char buffer[30];
        for (uint i = 0; i < columns_.size(); ++i) {
            sprintf(buffer, "x_%04d", i + 1);
            bc_x_[i].setName(buffer);
        }
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }

}

void BranchAndCut::createObjectiveFunction() {

    try{
        bc_objective_ = IloAdd(bc_model_, IloMinimize(bc_env_));

        IloExpr obj(bc_env_);
        for (uint i = 0; i < columns_.size(); ++i) {
                obj += columns_[i].cost * bc_x_[i];
        }
        bc_objective_.setExpr(obj);
        obj.end();
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }
}

void BranchAndCut::createConstraints(ColumnGeneration &cg) {

    try {
        int num_vertex = Instance::getInstance()->num_vertex;
        constraints_ = IloRangeArray(bc_env_, num_vertex + 1);

        // Criando restrições constraint_1_: degree constraints
        for (int i = 0; i < num_vertex; ++i) {
           IloExpr expr(bc_env_);

           for (int j = 0; j < (int)columns_.size(); ++j) {
               if(columns_[j].endpoint1 == i || columns_[j].endpoint2 == i)  //alguma extremidade da aresta analisada é o vértice i correspondente à restrição?
                   expr += bc_x_[j];
           }

           constraints_[i] = expr <= Instance::getInstance()->degree_constraint[i];

           bc_model_.add(constraints_[i]);
           expr.end();
        }
        /*-------------------------------------------------*/

        // Criando restrição constraint_2_: MST constraint
        IloExpr expr(bc_env_);
        for (uint i = 0; i < columns_.size(); ++i) {        // Somatorio das variáveis do mestre restrito
               expr += bc_x_[i];
        }

        constraints_[num_vertex] = expr == num_vertex - 1;       // Somatorio das variáveis no mestre restrito deve ser igual à |V| - 1

        bc_model_.add(constraints_[num_vertex]);
        expr.end();
        /*-------------------------------------------------*/

        // SEC's
        for(uint k = 0; k < cg.secs.size(); k++) {
           IloExpr expr(bc_env_);
           int sum_set = 0;        // Calcula o número de elementos no conjunto S da sec k

           for(int i=0; i < num_vertex; i++) {
               if(cg.secs[k]->S[i] == 1) {
                   sum_set++;
               }
               else
                   continue;
               for(int j=i+1; j < num_vertex; j++) {
                   if(cg.secs[k]->S[i] == 1 && cg.secs[k]->S[j] == 1 && edges_in_pool[i][j] != -1) {     //Ambas as extremidades i e j pertencem à S e a aresta {i,j} está no modelo?
                       expr += bc_x_[edges_in_pool[i][j]];       // Adiciona a variável de decisão referente a aresta {i,j} na nova restrição
                   }

               }
           }

           constraints_.add(expr <= sum_set - 1);
           bc_model_.add(constraints_[constraints_.getSize()-1]);

           expr.end();
        }
        /*-------------------------------------------------*/

        // Blossoms
        for(uint k = 0; k < cg.blossoms.size(); k++) {
           IloExpr expr(bc_env_);
           int sum_degrees;
           sum_degrees = 0;        // Soma dos graus dos vértices presentes em H
           for(uint i=0; i < cg.blossoms[k]->H.size(); i++) {
               if(cg.blossoms[k]->H[i] == 1) {     // O vértice i pertence à H?
                   sum_degrees += Instance::getInstance()->degree_constraint[i];
                   for(uint j=i+1; j < cg.blossoms[k]->H.size(); j++) {
                       if(cg.blossoms[k]->H[j] == 1)           // O vértice j pertence à H?
                           if(edges_in_pool[i][j] != -1)       // A aresta {i,j} está no modelo?
                               expr += bc_x_[edges_in_pool[i][j]];        // Inclui a aresta {i,j} na restrição que representa a blossom
                   }
               }
           }

           int cardT = 0;
           for(unsigned int i=0; i < cg.blossoms[k]->T.size(); i++) {
               if(edges_in_pool[cg.blossoms[k]->T[i].endpoint1][cg.blossoms[k]->T[i].endpoint2] != -1) {
                   expr += bc_x_[edges_in_pool[cg.blossoms[k]->T[i].endpoint1][cg.blossoms[k]->T[i].endpoint2]];          // Adiciona a aresta presente em T à restrição
                   cardT++;
               }
           }

           constraints_.add(expr <= floor(0.5*(sum_degrees + cardT)));
           bc_model_.add(constraints_[constraints_.getSize()-1]);
           expr.end();
        }
        /*-------------------------------------------------*/
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }

}

void BranchAndCut::getCgBasis(ColumnGeneration &cg) {

    try{
        IloCplex::BasisStatusArray basis_bc(bc_env_, bc_x_.getSize());

        // Obtem base da Geracao de colunas
        IloCplex::BasisStatusArray basis_var(cg.env_, cg.x_.getSize());
        IloCplex::BasisStatusArray basis_cons(cg.env_, Instance::getInstance()->num_vertex);
        IloCplex::BasisStatusArray basis_secs(cg.env_, cg.num_secs);
        IloCplex::BasisStatusArray basis_blossoms(cg.env_, cg.num_blossoms);

        cg.cplex_->getBasisStatuses(basis_var, cg.x_);
        cg.cplex_->getBasisStatuses(basis_cons, cg.degree_constraints_);
        basis_cons.add(cg.cplex_->getBasisStatus(cg.st_constraint_));
        cg.cplex_->getBasisStatuses(basis_secs, cg.secs_);
        basis_cons.add(basis_secs);
        cg.cplex_->getBasisStatuses(basis_blossoms, cg.blossoms_);
        basis_cons.add(basis_blossoms);
        /*-----------------------------------------*/

        for(uint i = 0; i < columns_.size(); i++) {
            int edp1 = columns_[i].endpoint1;
            int edp2 = columns_[i].endpoint2;
            if(cg.edges_in_pool[edp1][edp2] != -1) {
                if(cg.columns_[cg.edges_in_pool[edp1][edp2]]->fixed == false) {
                    basis_bc[i] = basis_var[cg.edges_in_pool[edp1][edp2]];
                }
            }
        }

        bc_cplex_->setBasisStatuses(basis_bc, bc_x_, basis_cons, constraints_);
    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
}

}

void BranchAndCut::setCplexParameters() {

    bc_cplex_->setParam(IloCplex::PreInd, 0);
    bc_cplex_->setParam(IloCplex::AggInd, 0);
    bc_cplex_->setParam(IloCplex::HeurFreq, -1);

    bc_cplex_->setParam(IloCplex::FracCuts, -1);
    bc_cplex_->setParam(IloCplex::FlowCovers, -1);
    bc_cplex_->setParam(IloCplex::GUBCovers, -1);
    bc_cplex_->setParam(IloCplex::Covers, -1);

    //versao 11
    bc_cplex_->setParam(IloCplex::ZeroHalfCuts, -1);
    bc_cplex_->setParam(IloCplex::ImplBd, -1);
    bc_cplex_->setParam(IloCplex::Cliques, -1);
    bc_cplex_->setParam(IloCplex::DisjCuts, -1);
    bc_cplex_->setParam(IloCplex::FlowPaths, -1);
    bc_cplex_->setParam(IloCplex::MIRCuts, -1);
    bc_cplex_->setParam(IloCplex::Threads, 1);

    bc_cplex_->setParam(IloCplex::EpGap, 1e-6);
    bc_cplex_->setDeleteMode(IloCplex::FixBasis);
}

void BranchAndCut::solveBc(double time_limit) {

    try {
			setCplexParameters();

			bc_cplex_->extract(bc_model_);

			bc_cplex_->add(bcMinCutCallback(bc_env_, bc_x_, columns_, edges_in_pool));

			bc_cplex_->add(bcBlossomCallback(bc_env_, bc_x_, columns_, edges_in_pool));

			bc_cplex_->add(kruskalxCallBack(bc_env_, bc_x_, columns_, edges_in_pool));

			bc_cplex_->setParam(IloCplex::TiLim, time_limit);

			bc_cplex_->solve();

			cout << "===> " << bc_cplex_->getObjValue() << " / " << bc_cplex_->getStatus() << " " << bc_cplex_->getNnodes() << endl;

    } catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    } catch (...) {
        cerr << "Error" << endl;
    }

}
