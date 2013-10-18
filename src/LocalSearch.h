/*
 * LocalSearch.h
 *
 *  Created on: Oct 10, 2013
 *      Author: luishbicalho
 */

#ifndef LOCALSEARCH_H_
#define LOCALSEARCH_H_

#include <vector>
#include "Edge.h"
#include "TreeGlover.h"

using namespace std;

class LocalSearch {
    public:
        vector<Edge> sol;
        double obj_func_value;
        vector<int> degrees;
        TreeGlover* tree_sol;

        LocalSearch(const vector<Edge> &s, double c);
        ~LocalSearch();
        void runLs(vector<Edge> &edges);
};


#endif /* LOCALSEARCH_H_ */
