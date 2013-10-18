/*
 * GeracaoColunas.cpp
 *
 *  Created on: 08/12/2012
 *      Author: luis
 */

#include "ColumnGeneration.h"

ColumnGeneration::ColumnGeneration() : model_(env_) {

	// Inicializa matriz de adjacencia das arestes no pool de colunas
	INITIAL_POOL_SIZE = Instance::getInstance()->num_vertex + 1;
	edges_in_pool.resize(Instance::getInstance()->num_vertex);
	for(int i=0; i<Instance::getInstance()->num_vertex; i++)
		edges_in_pool[i].resize(Instance::getInstance()->num_vertex, -1); 
	/*-------------*/

	time_cg = num_secs = num_blossoms = fixed = cont_var_fixed = cont_secs_rm = cont_debug = 0;
	best_dual = 0;

	try {
		cplex_ = new IloCplex(model_);
	} catch (IloException& e) {
		cerr << "ERROR: " << e.getMessage() << endl;
	} catch (...) {
		cerr << "Error" << endl;
	}

	for(int i = 0; i < Instance::getInstance()->num_vertex; ++i)
		for(int j = i + 1; j<Instance::getInstance()->num_vertex; ++j) {
			if(Instance::getInstance()->degree_constraint[i] != 1 || Instance::getInstance()->degree_constraint[j] != 1) {		// A aresta pode estar contida e uma solução para o DCMST?
				Edge edge(i, j, Instance::getInstance()->cost_matrix[i][j]);
				all_edges.push_back(edge);
			}
			else {
			    fixed++;
			}
		}

	G = new Graph(Instance::getInstance()->num_edges);
	for(int i = 0; i < (int)all_edges.size(); ++i) {
		G->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, Instance::getInstance()->cost_matrix[all_edges[i].endpoint1][all_edges[i].endpoint2]);
	}

	tblossom = tsec = tpl = tpric = tprimal = tplsec = 0;
	debug.open("sintese.txt", std::ios_base::app);

}

ColumnGeneration::~ColumnGeneration() {

	delete cplex_;
	delete G;
	debug.close();

}

void ColumnGeneration::createVariables() {

    x_ = IloNumVarArray(env_, (int)columns_.size(), 0, 1);

    char buffer[30];

    // renomeia variaveis x
    for (int i = 0; i < (int)columns_.size(); ++i) {
        sprintf(buffer, "x_%04d", i + 1);
        x_[i].setName(buffer);
    }


}

void ColumnGeneration::createObjectiveFunction() {

	objective_ = IloAdd(model_, IloMinimize(env_));

	IloExpr obj(env_);

	for (int i = 0; i < (int) columns_.size(); ++i) {
		obj += columns_[i]->cost * x_[i];
	}

	objective_.setExpr(obj);
	obj.end();

}

void ColumnGeneration::createConstraints() {

	int num_vertex = Instance::getInstance()->num_vertex;
	int i, j;

	// Criando restrições constraint_1_: degree constraints
	degree_constraints_ = IloRangeArray(env_, num_vertex);

	for (i = 0; i < num_vertex; ++i) {
		IloExpr expr(env_);

		for (j = 0; j < (int)columns_.size(); ++j) {
			if(columns_[j]->endpoint1 == i || columns_[j]->endpoint2 == i)	//alguma extremidade da aresta analisada é o vértice i correspondente à restrição?
				expr += x_[j];
		}

		degree_constraints_[i] = expr <= Instance::getInstance()->degree_constraint[i];

		char buffer[20];
		sprintf(buffer, "R1_%03d", i);
		degree_constraints_[i].setName(buffer);

		model_.add(degree_constraints_[i]);
		expr.end();
	}
	/*-------------------------------------------------*/

	// Criando restrição constraint_2_: MST constraint
	IloExpr expr(env_);
	for (i = 0; i < (int)columns_.size(); ++i) {		// Somatorio das variáveis do mestre restrito
		expr += x_[i];
	}
	st_constraint_ = expr == num_vertex - 1;				// Somatorio das variáveis no mestre restrito deve ser igual à |V| - 1

	st_constraint_.setName("R2____");
	model_.add(st_constraint_);
	expr.end();
	/*-------------------------------------------------*/

	secs_ = IloRangeArray(env_);

	blossoms_ = IloRangeArray(env_);
}

void ColumnGeneration::createModel() {

	createVariables();
	createObjectiveFunction();
	createConstraints();

}

void ColumnGeneration::generateInitialColumns() {


//	for(int i = 0; i < (int)all_edges.size(); ++i) {
//		G->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, Instance::getInstance()->cost_matrix[all_edges[i].endpoint1][all_edges[i].endpoint2]);
//		Edge *edge = new Edge(all_edges[i].endpoint1, all_edges[i].endpoint2, all_edges[i].cost);
//		columns_.push_back(edge);
//		edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2] = edges_in_pool[all_edges[i].endpoint2][all_edges[i].endpoint1] = i;
//	}
//	return;

    // Gera a AGM e AG com grau restrito heurístico e inclui as colunas desta no pool de colunas inicial
	vector<KEdge> AGDegreeConstrained = kruskalDegreeConstrained(G);
	vector<KEdge> AGM = kruskal(G);

	for(int i = 0; i < (int)AGDegreeConstrained.size(); ++i) {
		Edge *e = new Edge();
		e->endpoint1 = AGDegreeConstrained[i].p;
		e->endpoint2 = AGDegreeConstrained[i].q;
		e->cost = Instance::getInstance()->cost_matrix[e->endpoint1][e->endpoint2];
		columns_.push_back(e);
		edges_in_pool[e->endpoint1][e->endpoint2] = edges_in_pool[e->endpoint2][e->endpoint1] = (int)columns_.size() - 1;
	}

	for(int i = 0; i < (int)AGM.size(); ++i) {
		Edge *e = new Edge();
		e->endpoint1 = AGM[i].p;
		e->endpoint2 = AGM[i].q;
		e->cost = Instance::getInstance()->cost_matrix[e->endpoint1][e->endpoint2];
		if(edges_in_pool[e->endpoint1][e->endpoint2] == -1) {
			columns_.push_back(e);
			edges_in_pool[e->endpoint1][e->endpoint2] = edges_in_pool[e->endpoint2][e->endpoint1] = (int)columns_.size() - 1;
		}
	}

	/*-------------------------------------------------------------*/

	// Inclui mais colunas...
	sort(all_edges.begin(), all_edges.end());
	vector<int> counter(Instance::getInstance()->num_vertex, 0);

	for(int i=0; i<(int)all_edges.size(); i++) {
		if(counter[all_edges[i].endpoint1] < 3 || counter[all_edges[i].endpoint2] < 3) {
			Edge *e = new Edge();
			e->endpoint1 = all_edges[i].endpoint1;
			e->endpoint2 = all_edges[i].endpoint2;
			e->cost = all_edges[i].cost;
			if(edges_in_pool[e->endpoint1][e->endpoint2] == -1 && edges_in_pool[e->endpoint2][e->endpoint1] == -1) {
				++counter[e->endpoint1];
				++counter[e->endpoint2];
				columns_.push_back(e);
				edges_in_pool[e->endpoint1][e->endpoint2] = edges_in_pool[e->endpoint2][e->endpoint1] = (int)columns_.size() - 1;
			}
		}
	}
	/*---------------------------------------------------------------*/

}

void ColumnGeneration::getDualPrices(IloNumArray & u, IloNum & m, IloNumArray &p, IloNumArray &b) {

	try {
		cplex_->getDuals(u, degree_constraints_);	// Obtem o valor das variáveis duais correspondentes às restrições de grau
		cplex_->getDuals(p, secs_);                // Obtem o valor das variáveis duais correspondentes às SEC's
		cplex_->getDuals(b, blossoms_);	// Obtem o valor das variáveis duais correspondentes às BI's
		m = cplex_->getDual(st_constraint_);	// Obtem o valor da variável dual correspondente a restrição de Arvore Geradora

	} catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}

}

double ColumnGeneration::pricing(const IloNumArray & u, const IloNum & m, const IloNumArray &p, IloNumArray &b, Columns & new_columns) {

	int i;

	Edge best_edge(0, 0, 999999999);		// Melhor aresta a ser retornada pelo pricing. Custo inicial = infinito

	int cont = 0;
	for(i=0; i < (int)all_edges.size(); i++) {
		if(edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2] == -1) {		// Para cada aresta não presente no pool de colunas:
			double sum_p = 0;				// Soma das variáveis duais referente as retriçoes SEC
			double sum_b = 0;				// Soma das variáveis duais referente as BI's

			for(int k=0; k<(int)secs.size(); k++) {
				if(secs[k]->S[all_edges[i].endpoint1] == 1 && secs[k]->S[all_edges[i].endpoint2] == 1)		// Ambas as extremidades pertencem ao conj. S da SEC k?
					sum_p += p[k];
			}

			for(int k=0; k<(int)blossoms.size(); k++) {
				if(blossoms[k]->H[all_edges[i].endpoint1] == 1 && blossoms[k]->H[all_edges[i].endpoint2] == 1)		// Ambas as extremidades pertencem ao handle H da blossom k?
					sum_b += b[k];
			}

			all_edges[i].cost_reduced = all_edges[i].cost - u[all_edges[i].endpoint1] - u[all_edges[i].endpoint2] - m - sum_p - sum_b;			// Calculo do custo reduzido

			if(Utils::cmp(all_edges[i].cost_reduced, 0, 1e-4) < 0)
				cont++;
		}
		else {		// A aresta {i, j} já está incluída no mestre
			all_edges[i].cost_reduced = 0;
		}
	}

	sort(all_edges.begin(), all_edges.end());

	for(i = 0; i < ceil(cont*FACTOR_PRIC); i++) {
		Edge *edge = new Edge(all_edges[i].endpoint1, all_edges[i].endpoint2, all_edges[i].cost);
		new_columns.push_back(edge);
	}
	return all_edges[0].cost_reduced;
}

void ColumnGeneration::addNewColumns(const Columns & new_columns, int &cont_cols) {

	int num_columns = static_cast<int> (new_columns.size());

	for (int i = 0; i < num_columns; ++i) {
		cont_cols++;

		IloNumVar new_variable = IloNumVar(env_, 0, 1);
		Edge *new_col = new Edge(new_columns[i]->endpoint1, new_columns[i]->endpoint2, new_columns[i]->cost);
		// Insere a nova coluna no pool de colunas
		columns_.push_back(new_col);

		// indice da nova coluna para indexar nas restricoes
		int new_column_index = static_cast<int> (columns_.size() - 1);

		// Insere a nova coluna no pool de colunas
		edges_in_pool[new_columns[i]->endpoint1][new_columns[i]->endpoint2] = new_column_index;
		edges_in_pool[new_columns[i]->endpoint2][new_columns[i]->endpoint1] = new_column_index;

		// Muda o label da variavel
		//testes
		char buffer[30];
		sprintf(buffer, "x_%04d", new_column_index + 1);
		new_variable.setName(buffer);

		// adiciona nova coluna no pool de variavel do cplex
		x_.add(new_variable);

		// adiciona variavel no modelo do cplex
		model_.add(new_variable);

		// Altera o valor do coeficiente da nova variável na função objetivo
		objective_.setLinearCoef(x_[new_column_index], columns_[new_column_index]->cost);

		// seta coeficiente das restrições de grau (restrições 1)
		for (int j = 0; j < Instance::getInstance()->num_vertex; ++j) {
			if(columns_[new_column_index]->endpoint1 == j || columns_[new_column_index]->endpoint2 == j)
				degree_constraints_[j].setLinearCoef(x_[new_column_index], 1);
		}

		// seta coeficiente da restrição de arv. geradora (restrição 2)
		st_constraint_.setLinearCoef(x_[new_column_index], 1);

		// seta coeficiente das restrições secs
		IloNumArray aux(env_, secs.size());
		for(int i=0; i<(int)secs.size(); i++) {
			// Testa se nova coluna possui ambas extremidades no conjunto S avaliado
			if(secs[i]->S[columns_[new_column_index]->endpoint1] == 1 && secs[i]->S[columns_[new_column_index]->endpoint2] == 1)
			    secs_[i].setLinearCoef(x_[new_column_index], 1);
		}
		/*--------------------------------------*/

		// seta coeficiente das restrições de blossom
		for(int i=0; i<(int)blossoms.size(); i++) {
			// Testa se nova coluna possui ambas extremidades no conjunto H avaliado
			if(blossoms[i]->H[columns_[new_column_index]->endpoint1] == 1 && blossoms[i]->H[columns_[new_column_index]->endpoint2] == 1) {
				blossoms_[i].setLinearCoef(x_[new_column_index], 1);
			}
		}
	}

}

void ColumnGeneration::addNewSECs(vector<vector<int> > subset) {

	for(int k=0; k<(int)subset.size(); k++) {
		IloExpr expr(env_);
		int sum_set = 0;		// Calcula o número de elementos no conjunto S da sec k
		for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
			if(subset[k][i] == 1)
				sum_set++;
			else
				continue;
			for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
				if(subset[k][i] == 1 && subset[k][j] == 1 && edges_in_pool[i][j] != -1)		//Ambas as extremidades i e j pertencem à S e a aresta {i,j} está no modelo?
					expr += x_[edges_in_pool[i][j]];		// Adiciona a variável de decisão referente a aresta {i,j} na nova restrição
			}
		}

		Sec* s = new Sec(subset[k]);		// Cria uma nova sec
		secs.push_back(s);										// Adiciona a nova sec no pool de secs
		IloRange sec = expr <= (sum_set - 1);
		char buffer[20];
        sprintf(buffer, "R3_%04d", num_secs+1);
		sec.setName(buffer);
		secs_.add(sec);

		model_.add(secs_[num_secs++]);						// Inclui a nova sec no modelo
		expr.end();
	}

}

void ColumnGeneration::findNewSECs() {

	do {
		cont_debug++;

		IloNumArray xsol(env_, columns_.size());			// Valores de relaxação linear das variáveis x_
		cplex_->getValues(xsol, x_);

		vector<vector<int> > subsets, subsets_heur;
		int best_cut_id_heur = -1;
		int best_cut_id = heuristicKruskalSEC(G, xsol, edges_in_pool, subsets);

		if(subsets.size() == 0) {	// Se nenhum corte foi encontrado através da heurística baseada em kruskal...
			// A separação exada das SEC's é realizada
			MinCut mc;
			best_cut_id = mc.solveMinCutWolsey(xsol, edges_in_pool, subsets);		// Subconjuntos S de cortes violados

			if(best_cut_id >= 0) {		// Se foi encontrado cortes...
				// Mais uma separação é utilizada a fim de se obter cortes ortogonais
				// As arestas presente no corte mais violado através do método anterior são removidas
				MinCut mc_heur;
				best_cut_id_heur = mc_heur.solveMinCutWolseyOrth(xsol, edges_in_pool, subsets[best_cut_id], subsets_heur);
			}

		}

		vector<vector<int> > best_cuts;
		if(subsets.size() > 0) {
			best_cuts.push_back(subsets[best_cut_id]);
			double euc_norm_mv = Utils::euclidianNorm(subsets[best_cut_id], cplex_, xsol, edges_in_pool);	// Norma Euclidiana do corte mais violado
			for(int k=0; k<(int)subsets.size(); k++) {
				if(k != best_cut_id) {
					double euc_norm_cut = Utils::euclidianNorm(subsets[k], cplex_, xsol, edges_in_pool);						// Norma Euclidiana do corte atual
					double inner_product = 0;		// Valor do produto escalar
					for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
						for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
							if(subsets[k][i] == 1 && subsets[k][j] == 1 && subsets[best_cut_id][i] == 1 && subsets[best_cut_id][j] == 1) {	// A aresta está presente em ambas SECs?
								if(edges_in_pool[i][j] != -1)
									inner_product += xsol[edges_in_pool[i][j]]/euc_norm_mv * xsol[edges_in_pool[i][j]]/euc_norm_cut;
							}
						}
					}
					if(inner_product <= FACTOR_SEC) {
						best_cuts.push_back(subsets[k]);
					}
				}
			}

			if(subsets_heur.size() > 0) {
				best_cuts.push_back(subsets_heur[best_cut_id_heur]);
				euc_norm_mv = Utils::euclidianNorm(subsets_heur[best_cut_id_heur], cplex_, xsol, edges_in_pool);	// Norma Euclidiana do corte mais violado
				for(int k=0; k<(int)subsets_heur.size(); k++) {
					if(k != best_cut_id_heur) {
						double euc_norm_cut = Utils::euclidianNorm(subsets_heur[k], cplex_, xsol, edges_in_pool);						// Norma Euclidiana do corte atual
						double inner_product = 0;		// Valor do produto escalar
						for(int i=0; i < Instance::getInstance()->num_vertex; i++) {
							for(int j=i+1; j < Instance::getInstance()->num_vertex; j++) {
								if(subsets_heur[k][i] == 1 && subsets_heur[k][j] == 1 && subsets_heur[best_cut_id_heur][i] == 1 && subsets_heur[best_cut_id_heur][j] == 1) {	// A aresta está presente em ambas SECs?
									if(edges_in_pool[i][j] != -1)
										inner_product += xsol[edges_in_pool[i][j]]/euc_norm_mv * xsol[edges_in_pool[i][j]]/euc_norm_cut;
								}
							}
						}
						if(inner_product <= FACTOR_SEC) {
							best_cuts.push_back(subsets_heur[k]);
						}
					}
				}
			}

			addNewSECs(best_cuts);

			start = Utils::getTime();
				cplex_->solve();
			end = Utils::getTime();
			tplsec += end - start;

			manageRows();
			xsol.end();
		}
		else {
			break;
		}
	} while(true);

}

void ColumnGeneration::blossomHeuristicEdgeWeight1(IloNumArray &xsol, vector<vector<int> > &Hs, vector<vector<Edge> > &Ts) {

	vector<Edge> support_graph;			// Grafo suporte da solução. O custo da aresta é o valor de relaxação linear desta no mestre
	vector<vector<int> > adj_list(Instance::getInstance()->num_vertex);			// Lista de adjacência do grafo resultante da retirado das arestas com x_e = 1 do grafo suporte

	vector<double> s_line(Instance::getInstance()->num_vertex, 0.0);		// folga do vértice i (b|i| - x(delta(i)))
	vector<double> b_line(Instance::getInstance()->num_vertex);				// grau do vértice i  (b|i| - )

	for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
		b_line[i] = Instance::getInstance()->degree_constraint[i];
		s_line[i] = Instance::getInstance()->degree_constraint[i];
	}

				// Valores de relaxação linear das variáveis x_
	for(unsigned int i = 0; i<columns_.size(); i++) {
		if(Utils::cmp(xsol[i], 0) > 0) {				// A aresta referente à coluna i pertence à solução?
			Edge edge(columns_[i]->endpoint1, columns_[i]->endpoint2, xsol[i]);

			support_graph.push_back(edge);			// Adiciona a aresta ao grafo suporte

            s_line[edge.endpoint1] -= xsol[i];
            s_line[edge.endpoint2] -= xsol[i];


			if(Utils::cmp(xsol[i], 1) != 0) {		//O valor de relaxação linear é diferente de 1?
				adj_list[edge.endpoint1].push_back(edge.endpoint2);
				adj_list[edge.endpoint2].push_back(edge.endpoint1);
			}
			else {
				b_line[edge.endpoint1] -= 1;
				b_line[edge.endpoint2] -= 1;
			}
		}
	}

	vector<vector<int> > H;		// Conjuntos Handle - Cada conjunto handle é um vector no num_vertex posições onde H[i][j] = 1 indica que o vértice j
								// pertence ao handle i, e H[i][j] = 0 caso contrário
	vector<vector<Edge> > T;	// Conjuntos Tooth

	Utils::findHandles(adj_list, H);

	T.resize(H.size());
	vector<int> violated_blossom;		// Conjunto de indices de BI's violadas em H

	for(unsigned int i=0; i < H.size(); i++) {
        double sum_slack = 0;
        double sum_degree = 0;
        for(int j=0; j < Instance::getInstance()->num_vertex; j++) {
            if(H[i][j] == 1) {			// O vértice j pertence ao conjunto Handle i?
                sum_slack += s_line[j];
                sum_degree += b_line[j];
            }
        }

        if(Utils::cmp(sum_slack, 1) < 0 && ((int)sum_degree % 2) != 0) {			// A blossom está violada?
            violated_blossom.push_back(i);
            for(unsigned int j = 0; j < support_graph.size(); j++) {
                bool test1 = (H[i][support_graph[j].endpoint1] == 1 && H[i][support_graph[j].endpoint2] != 1);	// endpoint1 pertence à H mas endpoint2 não?
                bool test2 = (H[i][support_graph[j].endpoint1] != 1 && H[i][support_graph[j].endpoint2] == 1);  // endpoint1 não pertence à H mas endpoint2 pertence?
                if(test1 || test2) {	// Se somente uma das extremidades pertence à H...
                    T[i].push_back(support_graph[j]);
                }
            }
        }
	}

    Hs.clear();
    Ts.clear();
    Hs.resize(violated_blossom.size());
    Ts.resize(violated_blossom.size());
    for(unsigned int i=0; i<violated_blossom.size(); i++) {
        int index = violated_blossom[i];
        Hs[i] = H[index];
        Ts[i] = T[index];
    }

}

void ColumnGeneration::blossomHeuristicBridge(IloNumArray &xsol, vector<vector<int> > &Hs, vector<vector<Edge> > &Ts) {

	vector<Edge> support_graph;
	vector<vector<int> > bridges (Instance::getInstance()->num_vertex);
	vector<double> s_line(Instance::getInstance()->num_vertex);		        // folga do vértice i (b|i| - x(delta(i)))
	vector<double> b_line(Instance::getInstance()->num_vertex);				// grau do vértice i  (b|i| - |e \in \delta(i) : x_e = 1|)

	for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
		b_line[i] = Instance::getInstance()->degree_constraint[i];
		s_line[i] = Instance::getInstance()->degree_constraint[i];
	}

	for(unsigned int i = 0; i<columns_.size(); i++) {
		if(Utils::cmp(xsol[i], 0) > 0) {				// A aresta referente à coluna i pertence à solução?
			Edge edge(columns_[i]->endpoint1, columns_[i]->endpoint2, xsol[i]);
			support_graph.push_back(edge);			// Adiciona a aresta ao grafo suporte

			s_line[edge.endpoint1] -= xsol[i];
			s_line[edge.endpoint2] -= xsol[i];

			if(Utils::cmp(xsol[i], 1) == 0) {		// O valor de relaxação linear é igual de 1?
				b_line[edge.endpoint1] -= 1;
				b_line[edge.endpoint2] -= 1;
			}
		}
	}

	for(int i=0; i<Instance::getInstance()->num_vertex; i++)
		bridges[i].resize(Instance::getInstance()->num_vertex, 0);

	Utils::findBridges(support_graph, bridges);

	vector<Edge> new_graph;
	vector<vector<int> > adj_list(Instance::getInstance()->num_vertex);
	for(unsigned int i=0; i<support_graph.size(); i++) {
		if(bridges[support_graph[i].endpoint1][support_graph[i].endpoint2] != 1) {	// Se a aresta i da árvore não for uma ponte...
			new_graph.push_back(support_graph[i]);
			adj_list[support_graph[i].endpoint1].push_back(support_graph[i].endpoint2);
			adj_list[support_graph[i].endpoint2].push_back(support_graph[i].endpoint1);
		}
	}

	vector<vector<int> > H;		// Conjuntos Handle
	vector<vector<Edge> > T;	// Conjuntos Tooth
	Utils::findHandles(adj_list, H);

	T.resize(H.size());
	vector<int> violated_blossom;		// Conjunto de indices de BI's violadas em H

	for(unsigned int i=0; i < H.size(); i++) {
//	    int cont = 0;
//        for(int j=0; j < Instance::getInstance()->num_vertex; j++)
//            if(H[i][j] == 1)
//                cont++;
	    /*if(cont < Instance::getInstance()->num_vertex)*/ {

            double sum_slack = 0;
            double sum_degree = 0;
            for(int j=0; j < Instance::getInstance()->num_vertex; j++) {
                if(H[i][j] == 1) {			// O vértice j pertence ao conjunto Handle i?
                    sum_slack += s_line[j];
                    sum_degree += b_line[j];
                }
            }


            if(Utils::cmp(sum_slack, 1) < 0 && ((int)sum_degree % 2) != 0) {			// A blossom está violada?
//                cout << "[b]" << sum_slack << " " << sum_degree << endl;
                violated_blossom.push_back(i);
                for(unsigned int j = 0; j < support_graph.size(); j++) {
                    bool test1 = (H[i][support_graph[j].endpoint1] == 1 && H[i][support_graph[j].endpoint2] != 1);	// endpoint1 pertence à H mas endpoint2 não?
                    bool test2 = (H[i][support_graph[j].endpoint1] != 1 && H[i][support_graph[j].endpoint2] == 1);  // endpoint1 não pertence à H mas endpoint2 pertence?
                    if(test1 || test2) {	// Se somente uma das extremidades pertence à H...
                        T[i].push_back(support_graph[j]);
                    }
                }
            }
	    }
	}

	Hs.clear();
	Ts.clear();
	Hs.resize(violated_blossom.size());
	Ts.resize(violated_blossom.size());
	for(unsigned int i=0; i<violated_blossom.size(); i++) {
		int index = violated_blossom[i];
//		addNewBlossom(H[index], T[index]);
        Hs[i] = H[index];
        Ts[i] = T[index];
	}
}

void ColumnGeneration::addNewBlossom(vector<int> &H, vector<Edge> &T) {

	IloExpr expr(env_);
	int sum_degrees;
	sum_degrees = 0;		// Soma dos graus dos vértices presentes em H

	for(unsigned int i=0; i < H.size(); i++) {
		if(H[i] == 1) {		// O vértice i pertence à H?
			sum_degrees += Instance::getInstance()->degree_constraint[i];
			for(unsigned int j=i+1; j < H.size(); j++) {
				if(H[j] == 1)			// O vértice j pertence à H?
					if(edges_in_pool[i][j] != -1)		// A aresta {i,j} está no modelo?
						expr += x_[edges_in_pool[i][j]];		// Inclui a aresta {i,j} na restrição que representa a blossom
			}
		}
	}

	for(unsigned int i=0; i < T.size(); i++) {
		expr += x_[edges_in_pool[T[i].endpoint1][T[i].endpoint2]];			// Adiciona a aresta presente em T (já está no modelo) à restrição
	}

	IloRange blossom;
	blossom = expr <= floor(0.5*(sum_degrees + T.size()));
	Blossom* b = new Blossom(H, T);		// Cria uma nova blossom
	blossoms.push_back(b);															// Adiciona a nova blossom ao conjunto de blossoms

	char buffer[20];
	sprintf(buffer, "R4_%04d", num_blossoms+1);
	blossom.setName(buffer);
	blossoms_.add(blossom);
	model_.add(blossoms_[num_blossoms++]);				// Adiciona a nova blossom ao modelo
	expr.end();
}

//-%%%%%%%%% alterado %%%%%%%%%%%%%%%-/
bool ColumnGeneration::findNewBlossoms() {

    vector<vector<int> > H;
    vector<vector<Edge> > T;
    IloNumArray xsol(env_, columns_.size());
    cplex_->getValues(xsol, x_);

    blossomHeuristicEdgeWeight1(xsol, H, T);

//	if(H.size() == 0)
//		blossomHeuristicBridge(xsol, H, T);

    for(unsigned int i = 0; i < H.size(); i++)
        addNewBlossom(H[i], T[i]);

    if(H.size() > 0) {
	    cplex_->solve();
	    return true;
    }

    return false;
}

void ColumnGeneration::manageRows() {

	IloRangeArray auxrange(env_);		// Array de restrições que armazena àquelas que estão folgadas

	for(int i = 0; i < num_secs; i++) {

		double slack = cplex_->getSlack(secs_[i]);
		if(fabs(slack) < ERROR) {		// A restrição i está justa?
			secs[i]->cont = 0;
		}
		else {
			secs[i]->cont++;
		}

		if(secs[i]->cont == NUM_ITER_SEC) {		// Se a restrição i não está justa a NUM_ITER_SEC iterações ela é removida do pool de SEC's
			auxrange.add(secs_[i]);
			secs_[i] = secs_[num_secs-1];
			secs_.setSize(num_secs-1);
			num_secs--;

			delete secs[i];
			secs[i] = secs[secs.size()-1];
			secs.pop_back();
			cont_secs_rm++;
		}
	}

	if(auxrange.getSize() > 0) {		// Existe restrições para serem removidas do modelo?
		model_.remove(auxrange);
		cplex_->solve();
	}
}

Solution ColumnGeneration::kruskalx() {

	IloNumArray xsol(env_, columns_.size());			// Valores de relaxação linear das variáveis x_
    cplex_->getValues(xsol, x_);

	Graph *G_bar = new Graph(Instance::getInstance()->num_edges);
	for(unsigned int i=0; i<all_edges.size(); ++i) {
		if(edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2] != -1) {	// Se a aresta estiver no pool de colunas...
			G_bar->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2,
					all_edges[i].cost*(1 - xsol[edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2]]));	// O custo é substituído pelo custo complementar através do valor de RL
		}
		else
			G_bar->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, all_edges[i].cost);
	}

	vector<KEdge> st_heur = kruskalDegreeConstrained(G_bar);
	Solution sol;
	sol.obj_fun_value = 0.0;
	for(unsigned int i=0; i < st_heur.size(); ++i) {
		int edp1 = st_heur[i].p;
		int edp2 = st_heur[i].q;
		double cost = Instance::getInstance()->cost_matrix[edp1][edp2];
		Edge edge(edp1, edp2, cost);
		sol.spanning_tree.push_back(edge);
		sol.obj_fun_value += cost;
	}

	return sol;
}

Solution ColumnGeneration::lagrangeanPrimalHeuristic() {

	IloNumArray xsol(env_, columns_.size());			// Valores de relaxação linear das variáveis x_
	IloNumArray degree_duals(env_, Instance::getInstance()->num_vertex);
	IloNumArray blossoms_duals(env_, blossoms.size());

	cplex_->getDuals(degree_duals, degree_constraints_);	// Obtem o valor das variáveis duais correspondentes às restrições de grau
	cplex_->getDuals(blossoms_duals, blossoms_);           // Obtem o valor das variáveis duais correspondentes às BI's
    cplex_->getValues(xsol, x_);

	Graph *G_lagrangean = new Graph(Instance::getInstance()->num_edges);
	for(unsigned int i=0; i < all_edges.size(); ++i) {
		double lagrangean_cost = all_edges[i].cost;
        lagrangean_cost -= degree_duals[all_edges[i].endpoint1];
        lagrangean_cost -= degree_duals[all_edges[i].endpoint2];

		for(unsigned int j=0; j < blossoms.size(); j++) {
			if(blossoms[j]->H[all_edges[i].endpoint1] == 1 && blossoms[j]->H[all_edges[i].endpoint2] == 1) {
				lagrangean_cost -= blossoms_duals[j];
			}
			else {
				for(unsigned int k = 0; k < blossoms[j]->T.size(); k++) {
					if(blossoms[j]->T[k].endpoint1 == all_edges[i].endpoint1 && blossoms[j]->T[k].endpoint2 == all_edges[i].endpoint2) {
						lagrangean_cost -= blossoms_duals[j];
					}
					else if(blossoms[j]->T[k].endpoint1 == all_edges[i].endpoint2 && blossoms[j]->T[k].endpoint2 == all_edges[i].endpoint1) {
						lagrangean_cost -= blossoms_duals[j];
					}
				}
			}
		}

		G_lagrangean->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, lagrangean_cost);
	}

	vector<KEdge> st_lagrangean = kruskal(G_lagrangean);
	vector<vector<int> > st_lagr_adj_matrix(Instance::getInstance()->num_vertex);
	for(int i=0; i < Instance::getInstance()->num_vertex; i++)
		st_lagr_adj_matrix[i].resize(Instance::getInstance()->num_vertex, 0);

	for(unsigned int i=0; i < st_lagrangean.size(); i++) {
		st_lagr_adj_matrix[st_lagrangean[i].p][st_lagrangean[i].q] = 1;
		st_lagr_adj_matrix[st_lagrangean[i].q][st_lagrangean[i].p] = 1;
	}

	Graph *G_dg = new Graph(Instance::getInstance()->num_edges);
	for(unsigned int i=0; i < all_edges.size(); ++i) {
		if(st_lagr_adj_matrix[all_edges[i].endpoint1][all_edges[i].endpoint1] == 1) {
			G_dg->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, 0);
		}
		else {
			G_dg->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, all_edges[i].cost);
		}
	}

	vector<KEdge> st_heur = kruskalDegreeConstrained(G_dg);
	Solution sol;
	sol.obj_fun_value = 0.0;
	for(unsigned int i=0; i < st_heur.size(); ++i) {
		int edp1 = st_heur[i].p;
		int edp2 = st_heur[i].q;
		double cost = Instance::getInstance()->cost_matrix[edp1][edp2];
		Edge edge(edp1, edp2, cost);
		sol.spanning_tree.push_back(edge);
		sol.obj_fun_value += cost;
	}

	sol = localSearchEdgeExchange(sol);
	if(sol.obj_fun_value < best_primal.obj_fun_value)
		best_primal = sol;

	return sol;
}

Solution ColumnGeneration::localSearchEdgeExchange(Solution s) {

	TreeGlover tree_sol(s.spanning_tree, s.obj_fun_value);
	vector<int> degrees(Instance::getInstance()->num_vertex, 0);
	for(unsigned int i=0; i < s.spanning_tree.size(); i++) {
		++degrees[s.spanning_tree[i].endpoint1];
		++degrees[s.spanning_tree[i].endpoint2];
	}

	for(unsigned int i=0; i < all_edges.size(); i++) {
		int edp1 = all_edges[i].endpoint1;
		int edp2 = all_edges[i].endpoint2;

		// Ambas as extremidades estão folgadas?
		if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2]) {
			Edge outgoing = tree_sol.performEdgeExchange(all_edges[i]);
			if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {		// A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
				// Calcula os graus decorrente da troca de arestas
				--degrees[outgoing.endpoint1];
				--degrees[outgoing.endpoint2];
				++degrees[all_edges[i].endpoint1];
				++degrees[all_edges[i].endpoint2];
				//------------------------------------------------
				i = 0;					// Inspeciona as arestas do começo (first improvement approach)
			}
		}
		// A extremidade 1 está apertada (x(\delta(edp1)) == d_{edp1}) e a extremidade 2 está folgada?
		else if(degrees[edp1] <= Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] < Instance::getInstance()->degree_constraint[edp2]) {
			Edge outgoing = tree_sol.performEdgeExchange(all_edges[i], 1);
			if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {		// A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
				// Calcula os graus decorrente da troca de arestas
				--degrees[outgoing.endpoint1];
				--degrees[outgoing.endpoint2];
				++degrees[all_edges[i].endpoint1];
				++degrees[all_edges[i].endpoint2];
				//------------------------------------------------
				i = 0;					// Inspeciona as arestas do começo (first improvement approach)
			}
		}
		// A extremidade 1 está folgada e a extremidade 2 está apertada (x(\delta(edp2)) == d_{edp2})?
		else if(degrees[edp1] < Instance::getInstance()->degree_constraint[edp1] && degrees[edp2] <= Instance::getInstance()->degree_constraint[edp2]) {
			Edge outgoing = tree_sol.performEdgeExchange(all_edges[i], 2);
			if(!(outgoing.endpoint1 == 0 && outgoing.endpoint2 == 0)) {		// A BL foi realizada, ou seja, a aresta retornada é diferente de {0, 0}?
				// Calcula os graus decorrente da troca de arestas
				--degrees[outgoing.endpoint1];
				--degrees[outgoing.endpoint2];
				++degrees[all_edges[i].endpoint1];
				++degrees[all_edges[i].endpoint2];
				//------------------------------------------------
				i = 0;					// Inspeciona as arestas do começo (first improvement approach)
			}
		}
	}

	Solution sol(tree_sol.getTree(), tree_sol.getTreeCost());
	return sol;

}

// Problema: quando a inserção de uma aresta da solução target possui as duas extremidades saturadas na solução atual!
Solution ColumnGeneration::runPathRelinking() {

	int sol1, sol2;
	sol1 = rand()%elite_solutions.size();
	do {
		sol2 = rand()%elite_solutions.size();
	} while (sol1 == sol2);

	sort(elite_solutions[sol1].spanning_tree.begin(), elite_solutions[sol1].spanning_tree.end());
	sort(elite_solutions[sol2].spanning_tree.begin(), elite_solutions[sol2].spanning_tree.end());

	cout << "Solução 1: " << elite_solutions[sol1].obj_fun_value << endl;
	for(unsigned int i=0; i < elite_solutions[sol1].spanning_tree.size(); i++)
		cout << "{" << elite_solutions[sol1].spanning_tree[i].endpoint1 << ", " << elite_solutions[sol1].spanning_tree[i].endpoint2 << "}" << endl;
	cout << "Solução 2: " << elite_solutions[sol2].obj_fun_value << endl;
	for(unsigned int i=0; i < elite_solutions[sol2].spanning_tree.size(); i++)
			cout << "{" << elite_solutions[sol2].spanning_tree[i].endpoint1 << ", " << elite_solutions[sol2].spanning_tree[i].endpoint2 << "}" << endl;

	cin.get();

 	PathRelinking PR(elite_solutions[sol1].spanning_tree, elite_solutions[sol2].spanning_tree);

	vector<Edge> sol_pr = PR.runPR();
	Solution sol(sol_pr, 0);
	return sol;

}

void ColumnGeneration::primalSolution() {

	Solution sol_kruskalx = kruskalx();
	if(sol_kruskalx.obj_fun_value < best_primal.obj_fun_value) {
		best_primal = sol_kruskalx;
		elite_solutions.push_back(best_primal);
	}
	Solution sol_local_search = localSearchEdgeExchange(sol_kruskalx);
	if(sol_local_search.obj_fun_value < best_primal.obj_fun_value) {
		best_primal = sol_local_search;
		elite_solutions.push_back(best_primal);
	}

	if(elite_solutions.size() > 1) {
	    int sol1, sol2;
        sol1 = rand()%elite_solutions.size();
        do {
            sol2 = rand()%elite_solutions.size();
        } while (sol1 == sol2);
        vector<KEdge> st = heuristicKruskalPrBased(G, elite_solutions[sol1].spanning_tree, elite_solutions[sol2].spanning_tree/*best_primal.spanning_tree*/);
        vector<Edge> st_aux(st.size());
        double cost_st_aux = 0;
        for(unsigned int i = 0; i < st.size(); i++) {
            st_aux[i].setEdge(st[i].p, st[i].q, st[i].w);
            cost_st_aux += st[i].w;

        }

        Solution prbased(st_aux, cost_st_aux);
        sol_local_search = localSearchEdgeExchange(prbased);

        if(sol_local_search.obj_fun_value < best_primal.obj_fun_value) {
            best_primal = sol_local_search;
            elite_solutions.push_back(best_primal);
        }
	}
}

int ColumnGeneration::variableFixingLPR() {

    int cont_edges = 0;

    IloNum st_duals;
    IloNumArray degree_duals(env_, Instance::getInstance()->num_vertex);
    IloNumArray secs_duals(env_, secs.size());
    IloNumArray blossoms_duals(env_, blossoms.size());

    getDualPrices(degree_duals, st_duals, secs_duals, blossoms_duals);

    for(int i = 0; i < (int)all_edges.size(); i++) {
        if(all_edges[i].fixed == false) {
            double sum_p = 0;               // Soma das variáveis duais referente as retriçoes SEC
            double sum_b = 0;               // Soma das variáveis duais referente as BI's

            for(int k=0; k<(int)secs.size(); k++) {
                if(secs[k]->S[all_edges[i].endpoint1] == 1 && secs[k]->S[all_edges[i].endpoint2] == 1)      // Ambas as extremidades pertencem ao conj. S da SEC k?
                    sum_p += secs_duals[k];
            }

            for(int k=0; k<(int)blossoms.size(); k++) {
                if(blossoms[k]->H[all_edges[i].endpoint1] == 1 && blossoms[k]->H[all_edges[i].endpoint2] == 1)      // Ambas as extremidades pertencem ao handle H da blossom k?
                    sum_b += blossoms_duals[k];
            }

            all_edges[i].cost_reduced = all_edges[i].cost - degree_duals[all_edges[i].endpoint1] - degree_duals[all_edges[i].endpoint2] - st_duals - sum_p - sum_b;       // Calculo do custo reduzido
            if(Utils::cmp(all_edges[i].cost_reduced, best_primal.obj_fun_value - best_dual) > 0) {
                cont_edges++;
                all_edges[i].fixed = true;
                if(edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2] != -1) {
                    columns_[edges_in_pool[all_edges[i].endpoint1][all_edges[i].endpoint2]]->fixed = true;
                }
            }
        }
    }

    return cont_edges;
}

int ColumnGeneration::variableFixingLagr() {

    double lagr_const = 0;
    vector<vector<double> > lagr_cost_matrix(Instance::getInstance()->num_vertex);
    for(int i=0; i < Instance::getInstance()->num_vertex; i++)
        lagr_cost_matrix[i].resize(Instance::getInstance()->num_vertex, 0.0);

    IloNumArray degree_duals(env_, Instance::getInstance()->num_vertex);
    IloNumArray blossoms_duals(env_, blossoms.size());

    cplex_->getDuals(degree_duals, degree_constraints_); // Obtem o valor das variáveis duais correspondentes às restrições de grau
    for(int i=0; i<Instance::getInstance()->num_vertex; i++) {
        lagr_const += degree_duals[i] * degree_constraints_[i].getUB();
    }

    cplex_->getDuals(blossoms_duals, blossoms_);       // Obtem o valor das variáveis duais correspondentes às BI's

    for(unsigned int i=0; i<blossoms.size(); i++) {
        lagr_const += blossoms_duals[i] * blossoms_[i].getUB();
    }

    Graph *G_lagrangean = new Graph(Instance::getInstance()->num_edges);
    for(unsigned int i=0; i < all_edges.size(); ++i) {
        double lagrangean_cost = all_edges[i].cost;
        lagrangean_cost -= degree_duals[all_edges[i].endpoint1];
        lagrangean_cost -= degree_duals[all_edges[i].endpoint2];

        for(unsigned int j=0; j < blossoms.size(); j++) {
            if(blossoms[j]->H[all_edges[i].endpoint1] == 1 && blossoms[j]->H[all_edges[i].endpoint2] == 1) {
                lagrangean_cost -= blossoms_duals[j];
            }
            else {
                for(unsigned int k = 0; k < blossoms[j]->T.size(); k++) {
                    if(blossoms[j]->T[k].endpoint1 == all_edges[i].endpoint1 && blossoms[j]->T[k].endpoint2 == all_edges[i].endpoint2) {
                        lagrangean_cost -= blossoms_duals[j];
                    }
                    else if(blossoms[j]->T[k].endpoint1 == all_edges[i].endpoint2 && blossoms[j]->T[k].endpoint2 == all_edges[i].endpoint1) {
                        lagrangean_cost -= blossoms_duals[j];
                    }
                }
            }
        }

        G_lagrangean->push_edge(all_edges[i].endpoint1, all_edges[i].endpoint2, lagrangean_cost);
        lagr_cost_matrix[all_edges[i].endpoint1][all_edges[i].endpoint2] = lagr_cost_matrix[all_edges[i].endpoint2][all_edges[i].endpoint1] = lagrangean_cost;
    }

    VariableFixing vf(G_lagrangean, lagr_cost_matrix, lagr_const);

    int num_edges = vf.runVariableFixing(best_primal.obj_fun_value, edges_in_pool, all_edges, columns_);

    return num_edges;
}

vector<int> ColumnGeneration::bfs(int source, char checked[], vector<vector<int> > &lista_adj) {

	vector<int> s;

	queue<int> Q;
	checked[source] = 'b';
	Q.push(source);
	s.push_back(source);
	while(Q.empty() == false) {
		int u = Q.front();
		Q.pop();
		for(int i=0; i< (int)lista_adj[u].size(); i++) {
			int v = lista_adj[u][i];
			if(checked[v] == 'w') {
				checked[v] = 'b';
				s.push_back(v);
				Q.push(v);
			}
		}
	}

	sort(s.begin(), s.end());

	return s;

}

vector<int> ColumnGeneration::checa_solucao_bfs (vector<Edge> v) {

	vector<int> melhor;
	double valor_melhor = 0;

	vector<vector<int> > lista_adj(Instance::getInstance()->num_vertex);
	char* checked = new char[Instance::getInstance()->num_vertex];
	for(int i = 0; i<Instance::getInstance()->num_vertex; i++)
		checked[i] = 'w';

	for(int i=0; i<(int)v.size(); i++) {
		lista_adj[v[i].endpoint1].push_back(v[i].endpoint2);
		lista_adj[v[i].endpoint2].push_back(v[i].endpoint1);
	}

	vector<int> s;
	int cont = 0;
	for(int i = 0; i<Instance::getInstance()->num_vertex; i++) {
		if(checked[i] == 'w') {
			cont++;
			double soma_s = 0;
			s = bfs(i, checked, lista_adj);
			cout << "Subset " << cont << ": ";
			for(int j=0; j< (int)s.size(); j++) {
				cout << s[j] << " ";
				for(int k =j+1; k<(int)s.size(); k++) {
					if(edges_in_pool[s[j]][s[k]]!=-1)
						soma_s += (double) cplex_->getValue(x_[edges_in_pool[s[j]][s[k]]]);
				}
			}
			cout << "\t sum_e(s) = " << soma_s << "\t|S| = " << s.size();
			if(soma_s > s.size() - 1 && (soma_s - s.size() + 1) > valor_melhor) {
				valor_melhor = soma_s - s.size();
				melhor = s;
			}
			cout << endl;
		}

	}

	return melhor;

}

void ColumnGeneration::setCplexParameters() {

	cplex_->setParam(IloCplex::PreInd, 0);
	cplex_->setParam(IloCplex::AggInd, 0);
	cplex_->setParam(IloCplex::HeurFreq, -1);

	cplex_->setParam(IloCplex::FracCuts, -1);
	cplex_->setParam(IloCplex::FlowCovers, -1);
	cplex_->setParam(IloCplex::GUBCovers, -1);
	cplex_->setParam(IloCplex::Covers, -1);

	//versao 11
	cplex_->setParam(IloCplex::ZeroHalfCuts, -1);
	cplex_->setParam(IloCplex::ImplBd, -1);
	cplex_->setParam(IloCplex::Cliques, -1);
	cplex_->setParam(IloCplex::DisjCuts, -1);
	cplex_->setParam(IloCplex::FlowPaths, -1);
	cplex_->setParam(IloCplex::MIRCuts, -1);
	cplex_->setParam(IloCplex::Threads, 1);

	cplex_->setDeleteMode(IloCplex::FixBasis);
//	master_->setDeleteMode(IloCplex::LeaveBasis);

	//metodo da barreira
//	master_->setParam(IloCplex::RootAlg, IloCplex::Barrier);


}

void ColumnGeneration::runCG(int time_limit) {

    // Armazena os valores das variaveis duais
	IloNumArray degree_duals(env_, Instance::getInstance()->num_vertex);
	IloNum st_duals;
	bool find_violated_blossoms;

	// Novas colunas a serem inseridas
	Columns new_columns;

	int cont_deg = 0;
	double value = 9999999;
	double cost_reduced;

	try {
		// Gera um numero pre-definido de colunas para o pool inicial
		generateInitialColumns();

		setlocale(LC_ALL, "ptb");
		tot_start = Utils::getTime();

		createModel();
		setCplexParameters();

		cplex_->extract(model_);

		int cont_cols = 0;		// Contador de colunas geradas por pricing


		debug << "------------------------------------------------------------------------------------------------------" << endl;
		debug << "@@ ---- " << Instance::getInstance()->instancename << " ---- @@" << endl;

		start = Utils::getTime();
			cplex_->solve();
		end = Utils::getTime();
		tpl += end - start;

        do {
        	start = Utils::getTime();
            	findNewSECs();
            end = Utils::getTime();
            tsec += end - start;

            start = Utils::getTime();
            	find_violated_blossoms = findNewBlossoms();
            end = Utils::getTime();
			tblossom += end - start;
        } while(find_violated_blossoms == true);
        debug << __LINE__ << "Separação!" << endl;

        start = Utils::getTime();
        	primalSolution();
        	lagrangeanPrimalHeuristic();
        end = Utils::getTime();
		tprimal += end - start;
		debug << __LINE__ << "Primal!" << endl;

        debug << Instance::getInstance()->instancename << endl;
		debug << "LB: " << best_dual << " UB: " << best_primal.obj_fun_value << "%GAP: " << (best_primal.obj_fun_value - best_dual)/best_primal.obj_fun_value
		      << " #cols: " << columns_.size() << " #secs: " << secs.size() << " #blossom: " << blossoms.size();

		do {

			// Armazena os valores das variaveis duais
			IloNumArray secs_duals(env_, secs.size());
			IloNumArray blossoms_duals(env_, blossoms.size());

			getDualPrices(degree_duals, st_duals, secs_duals, blossoms_duals);

			// Calculo  do Pricing
			cout << "Pricing..." << endl;
			start = Utils::getTime();
				cost_reduced = pricing(degree_duals, st_duals, secs_duals, blossoms_duals, new_columns);
			end = Utils::getTime();
			tpric += end - start;
			debug << __LINE__ << "Pricing!" << endl;

			double inf_lim = cplex_->getObjValue();
			inf_lim = cplex_->getObjValue() - cost_reduced * (Instance::getInstance()->num_vertex-1);
			if(best_dual < inf_lim)
				best_dual = inf_lim;

			if (new_columns.size() < 1) { // == 0
				best_dual = cplex_->getObjValue();
				if(best_primal.obj_fun_value - best_dual > 1) {
				    cont_var_fixed += variableFixingLPR();
                }
				break;
			}

			// Adiciona novas colunas ao modelo
			addNewColumns(new_columns, cont_cols);
			new_columns.clear();

			start = Utils::getTime();
			cplex_->solve();
			end = Utils::getTime();
			tpl += end - start;

			if(Utils::cmp(cplex_->getObjValue(), value, 1e-4) < 0)		// Há degeneração?
				value = Utils::cmp(cplex_->getObjValue());
			else
				cont_deg++;

			manageRows();

			cout << "Separacao..." << endl;
	        do {
	        	start = Utils::getTime();
	            	findNewSECs();
	            end = Utils::getTime();
	            tsec += end - start;

	            start = Utils::getTime();
	            	find_violated_blossoms = findNewBlossoms();
	            end = Utils::getTime();
				tblossom += end - start;
	        } while(find_violated_blossoms == true);
	        debug << __LINE__ << "Separação (loop)!" << endl;

	        cout << "Primal..." << endl;
	        start = Utils::getTime();
	        	primalSolution();
	        	lagrangeanPrimalHeuristic();
	        end = Utils::getTime();
			tprimal += end - start;
			debug << __LINE__ << "Primal (loop)!" << endl;

			if(best_dual < best_primal.obj_fun_value && (best_primal.obj_fun_value - best_dual)/best_primal.obj_fun_value < 0.10)
				cont_var_fixed += variableFixingLagr();

			double aux = Utils::getTime();
			if(Utils::cmp(aux - tot_start, (double)time_limit, 1e-2) >= 0) {
				best_dual = cplex_->getObjValue();
				best_dual -= cost_reduced * (Instance::getInstance()->num_vertex-1);
				if(best_primal.obj_fun_value - best_dual > 1) {
				    cont_var_fixed += variableFixingLPR();
                }
				break;
			}
		} while (true);

		tot_end = Utils::getTime();
		time_cg = tot_end - tot_start;

		debug << "LB: " << best_dual << " UB: " << best_primal.obj_fun_value << "%GAP: " << (best_primal.obj_fun_value - best_dual)/best_primal.obj_fun_value
		      << " #cols: " << columns_.size() << " #secs: " << secs.size() << " #blossom: " << blossoms.size();
		debug << "Tempo de separação das SECs: " << tsec << endl;
		debug << "Tempo de separação das blossoms: " << tblossom << endl;
		debug << "Tempo do pricing: " << tpric << endl;
		debug << "Tempo da heurística primais: " << tprimal << endl;
		debug << "Tempo de PL do mestre: " << tpl << endl;
		debug << "Tempo de PL das SECs: " << tplsec << endl;
		debug << "------------------------------------------------------------------------------------------------------" << endl;
	} catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	} catch (...) {
		cerr << "Error" << endl;
	}

}
