/*
 * TypedefsAndDefines.h
 *
 *  Created on: Sep 12, 2013
 *      Author: luishbicalho
 */

#ifndef TYPEDEFSANDDEFINES_H_
#define TYPEDEFSANDDEFINES_H_

#include <vector>
#include "Edge.h"

using std::vector;

#define FACTOR_SEC 0.2
#define FACTOR_PRIC 0.10
#define ERROR 0.00001
#define NUM_ITER_SEC 30
#define NUM_ITER_COL 10

#define INF 0x3F3F3F3F // Signed int
#define EPS 1e-10
#define PI 3.14159265358979323846

#define FORN(i, a, b) for(typeof(b) i = (a); i < (b); i++)
#define ALL(x) x.begin(), x.end()
#define FOREACH(i, c) for(typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)

typedef vector<Edge*> Columns;


#endif /* TYPEDEFSANDDEFINES_H_ */
