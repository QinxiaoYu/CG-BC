/*
 * Instance.h
 *
 *  Created on: 09/12/2012
 *      Author: luis
 */

#ifndef INSTANCE_H_
#define INSTANCE_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

/*
 *  --- Classe Instance ---
 *  Singleton que contém os dados da instância a ser submetida ao método
 * */
class Instance {

public:
	static Instance* instance;
	int num_vertex;         // número de vértices do grafo
	int num_edges;          // número de arestas do grafo
	vector<int> degree_constraint;          // cada posição [i] indica o número máximo de arestas que podem incidir (grau máximo) ao vértice i
	vector<vector<double> > cost_matrix;    // cada posição [i][j] indica o custo da aresta {i, j}
	string instancename;                    // nome do arquivo de entrada instância
	string outputname;                      // nome do arquivo de saída (resultados) da instância

	static Instance* getInstance();         // retorna o ponteiro referente ao objeto da instância
	void finalizeInstance();                // encerra a instância atual
	Instance();

	/*
	 *  Construtor
	 *  Input:
	 *          name: nome do arquivo com os dados da instância
	 *        output: nome do arquivo que será gravado os resultados para esta instância
	 *     tipo_inst: tipo da instância, formato pelo qual os dados estão no arquivo de entrada
	 *                   tipo_inst = 1 -> ANDINST
	 *                   tipo_inst = 2 e 3 -> DE e R123
	 * */
	Instance(string name, string output, int tipo_inst);
	~Instance();
	void factory(string name, string output, int tipo_inst);        // constói uma nova instância
	void printValues();                                             // imprime os atributos da instância


};

#endif /* INSTANCE_H_ */
