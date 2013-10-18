/*
 * Instance.cpp
 *
 *  Created on: 03/05/2013
 *      Author: luishbicalho
 */

#include "Instance.h"

Instance* Instance::instance = new Instance();

Instance* Instance::getInstance() {
	return instance;
}

Instance::Instance () {
	num_vertex = 0;
	num_edges = 0;
}

Instance::~Instance() {

};

Instance::Instance (string name, string output, int tipo_inst) {

	if(name == "inst10_1.txt" || name == "inst10_2.txt" || name == "inst5_1.txt") {
		instancename = name;
		outputname = "teste.txt";

		num_vertex = 0;
		num_edges = 0;
		ifstream dados(name.c_str());
		if(!dados) {
			cout << "Erro na leitura do arquivo." << endl;
			exit(1);
		}

		dados >> num_vertex;
		num_edges = (num_vertex * (num_vertex -1))/2;

		degree_constraint.resize(num_vertex);
		cost_matrix.resize(num_vertex);
		for(int i=0; i<num_vertex; i++)
			cost_matrix[i].resize(num_vertex, 0);

		int u;
		for(int i=0; i < num_vertex; i++) {
			dados >> u;
			degree_constraint[i] = u;
		}


		for(int i=0; i<num_vertex; i++)
			for(int j=0; j<num_vertex; j++)
				dados >> cost_matrix[i][j];

		for(int i=0; i < num_vertex; i++) {
			cout << degree_constraint[i] << endl;
		}


		for(int i=0; i<num_vertex; i++) {
			for(int j=0; j<num_vertex; j++)
				cout << cost_matrix[i][j] << " ";
			cout << endl;
		}
		return;
	}

	instancename = name;
	outputname = output;

	int trash;
	num_vertex = 0;
	num_edges = 0;

	ifstream dados(name.c_str());
	if(!dados) {
		cout << "Erro na leitura do arquivo." << endl;
		exit(1);
	}

	int u, v, count;
	double cost;
	switch(tipo_inst) {
	case 1:		//ANDINST

		dados >> num_vertex;
		dados >> num_edges;

		degree_constraint.resize(num_vertex);
		cost_matrix.resize(num_vertex);
		for(int i=0; i<num_vertex; i++)
			cost_matrix[i].resize(num_vertex);

		for(int i=0; i < num_vertex; i++)
			cost_matrix[i][i] = 0;

		count = 0;
		while(count != num_edges) {
			dados >> u;
			dados >> v;
			dados >> cost;
			cost_matrix[u - 1][v - 1] = cost;
			cost_matrix[v - 1][u - 1] = cost;
			++count;
		}

		for(int i=0; i < num_vertex; i++) {
			dados >> trash;
			dados >> u;
			degree_constraint[i] = u;
		}
		break;
	case 2:		//DE
	case 3:		//R123
		dados >> trash;
		dados >> num_vertex;
		dados >> num_edges;

		degree_constraint.resize(num_vertex);
		cost_matrix.resize(num_vertex);
		for(int i=0; i<num_vertex; i++)
			cost_matrix[i].resize(num_vertex);

		for(int i=0; i < num_vertex; i++)
			cost_matrix[i][i] = 0;

		for(int i=0; i < num_vertex; i++) {
			dados >> u;
			degree_constraint[i] = u;
		}

		count = 0;
		while(count != num_edges) {
			dados >> u;
			dados >> v;
			dados >> cost;
			cost_matrix[u - 1][v - 1] = cost;
			cost_matrix[v - 1][u - 1] = cost;
			++count;
		}
		break;
	}
	dados.close();

}

void Instance::factory(string name, string output, int tipo_inst) {

	delete instance;
	instance = new Instance(name, output, tipo_inst);

}

void Instance::finalizeInstance() {

	delete instance;

}

void Instance::printValues() {
	cout << "|V| = " << num_vertex << endl;
	cout << "|E| = " << num_edges << endl;
	cin.get();
	cout << "Graus: " << endl;
	for(int i=0; i<num_vertex; i++)
		cout << "v = "  << i+1 << ": d_v = " << degree_constraint[i] << endl;
	cin.get();
	for(int i=0; i<num_vertex; i++)
		for(int j=i+1; j<num_vertex; j++) {
			cout << i+1 << " " << j+1 << " " << cost_matrix[i][j] << endl;
		}
	cin.get();
	cout << endl;
}
