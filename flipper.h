#pragma once
#include <vector>
#include <galois++/element.h>
#include "matrix.h"

class Flipper {

private:
	int numEqns, numVars, numNodes;
	// x is the current variables, b is the desired solution as in Ax=b, xiPrime is the "ideal" solution
	std::vector<Galois::Element> x, b, xiPrime;
	Galois::Field* gf;
	// the parity check matrix
	Matrix *matrix;
	// t is the number of satisfied equations for each 
	std::vector<int> t;
	// beta is the reliability ratio, defined to be ti/di where di is the degree
	std::vector<double> beta;

	// variables involved in path finding
	// stored as al[u] = {v,var, deg1node} where var is the intermediate degree 2 variable node, deg1node is the degree1node if applicable, -1 otherwise.
	std::vector<std::vector<std::tuple<int, int, int>>> al; 
	// distance from one eqn node to another eqn node, for APSP
	std::vector<std::vector<int>> d;
	// next is the next node in the path
	std::vector<std::vector<std::tuple<int, int, int>>> next;
	double maxBeta;

public:
	Flipper(Matrix* m) {
		matrix = m;

		numEqns = m->getNumRows();
		numVars = m->getNumCols();
		numNodes = numEqns + numVars;

		gf = m->getField();
		x.assign(numVars, Galois::Element(gf,0));
		xiPrime.assign(numVars, Galois::Element(gf, 0));
		t.resize(numVars);
		beta.resize(numVars);
		al.resize(numEqns);
		b.assign(numEqns, Galois::Element(gf, 0));
		d.assign(numEqns, std::vector<int>(numEqns, 1e9));
		next.assign(numEqns, std::vector<std::tuple<int,int,int>>(numEqns, { -1,-1,-1 }));
	}

	void setupAL();

	void floydWarshall();
	std::vector<std::tuple<int, int, int>> getPath(int u, int v);

	void assignX(std::vector<int> xs);
	void assignB(std::vector<int> bs);

	void update(int i);
	void getXiPrimeTiBi(int i);
	int getMaxBetaIndex();
	void flip(int i);

	// solving the equations using extended bit flipping algorithm
	void solve_extended_bit_flipping_consecutively();
	void solve_extended_bit_flipping_with_Omega();

	// solving the equations using multiple threshold decoding
	bool solve_single_threshold_decoding(int threshold);
	bool solve_multiple_threshold_decoding(std::vector<int> thresholds);

	friend std::ostream& operator<<(std::ostream& output, Flipper& flipper);

};