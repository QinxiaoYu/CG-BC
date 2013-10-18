/*
 * main.cpp
 *
 *  Created on: 08/12/2012
 *      Author: luis
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Instance.h"
#include "ColumnGeneration.h"
#include "BranchAndCut.h"

#include <ilcplex/ilocplex.h>
#define time_limit 14400

using namespace std;

int main(int argc, char* argv[]) {

	string inst_path;
	string outputname;
	int inst_type = 0;

	if(argc > 1) {
		string aux(argv[1]);
		inst_path = aux;
		if(argc > 2) {
			string aux2(argv[2]);
			outputname = aux2;
		}
		else
			outputname = "saida.txt";
		if(argc > 3) {
			inst_type = atoi(argv[3]);
		}
	}
	else {
		//		inst_path = "/home/luishbicalho/Dropbox/Projeto/GC/allinstances/ANDINST/tb1ct100_1.txt"; inst_type = 1;
		//		inst_path = "/compartilhamentos/homes/luishbicalho/Dropbox/Projeto/GC/allinstances/ANDINST/tb1ct100_1.txt"; inst_type = 1;
		inst_path = "/home/luishbicalho/Dropbox/Projeto/GC/allinstances/DE/dcmst100_1.dat"; inst_type = 2;
		//		inst_path = "/compartilhamentos/homes/luishbicalho/Dropbox/Projeto/GC/allinstances/DE/dcmst100_1.dat"; inst_type = 2;
		//		inst_path = "inst10_2.txt";
		//		inst_path = "inst5_1.txt";
		outputname = "saida.txt";
	}

	srand(50);

	cout << inst_path << "..." << inst_type << endl;

	Instance::getInstance()->factory(inst_path, outputname, inst_type);

	ofstream results(Instance::getInstance()->outputname.c_str(), std::ios_base::app);

	ColumnGeneration cg;

	cg.runCG(time_limit);
	BranchAndCut bc(cg);

	double time_bc, start_bc, end_bc;
	time_bc = start_bc = end_bc = 0;

	double tl_bc = time_limit - cg.time_cg;
	if(cg.best_primal.obj_fun_value - cg.best_dual > 1 && tl_bc > 0) {
        start_bc = Utils::getTime();
        bc.solveBc(tl_bc);
        end_bc = Utils::getTime();
        time_bc = end_bc - start_bc;
        results << inst_path << " " << std::fixed << cg.best_dual << " " << cg.best_primal.obj_fun_value << " " << setprecision(3) << 100*(cg.best_primal.obj_fun_value - cg.best_dual)/cg.best_primal.obj_fun_value
                << " " << cg.columns_.size() << " " << (double)100*cg.cont_var_fixed/Instance::getInstance()->num_edges << " " << cg.time_cg << " " << cg.cplex_->getStatus() << " "
                << bc.bc_cplex_->getNnodes() << " " << time_bc << " " << bc.bc_cplex_->getObjValue() << " "  << time_bc + cg.time_cg << " " << bc.bc_cplex_->getStatus() << endl;
	}
	else {
	    results << inst_path << " " << std::fixed << cg.best_dual << " " << cg.best_primal.obj_fun_value << " " << setprecision(3) << 100*(cg.best_primal.obj_fun_value - cg.best_dual)/cg.best_primal.obj_fun_value
                << " " << cg.columns_.size() << " " << (double)100*cg.cont_var_fixed/Instance::getInstance()->num_edges << " " << cg.time_cg << " " << cg.cplex_->getStatus() << " "
                << "0" << " " << "0" << " " << cg.best_primal.obj_fun_value << " "  << time_bc + cg.time_cg  << endl;
	}

	Instance::getInstance()->finalizeInstance();

	return 0;
}
