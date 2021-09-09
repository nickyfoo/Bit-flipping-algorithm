#include <iostream>
#include <fstream>
#include <galois++/field.h>
#include <galois++/element.h>
#include "matrix.h"
#include "flipper.h"
#include <random>

// For taking input
/*
int main(int argc, char *argv[]) {
	std::ifstream inputfile;
	inputfile.open("./sample_F2.txt", std::ios::in);
	if (inputfile.is_open()) {
		int fieldOrder; inputfile >> fieldOrder;
		Galois::Field gf(fieldOrder);
		std::cout << gf << std::endl;
		int row, col; inputfile >> row >> col;
		Matrix m(row, col, &gf);
		int n; inputfile >> n;
		for (int line = 0; line < n; line++) {
			int i, j, x;
			inputfile >> i >> j >> x;
			m.insert(i, j, Galois::Element(&gf, x));
		}
		Flipper f(&m);
		std::vector<int> x(col);
		for (int i = 0; i < col; i++) {
			inputfile >> x[i];
		}
		f.assignX(x);
		std::vector<int> b(row);
		for (int i = 0; i < row; i++) {
			inputfile >> b[i];
		}
		f.assignB(b); 
		f.FloydWarshall();
		m.findSingleVar();
		std::cout << m;
		f.solve_extended_bit_flipping_consecutively();
		////f.solve_consequently();
		//std::vector<int> thresholds;
		//for (int i = 0; i < row; i++) thresholds.push_back(i);
		//if (f.solve_multiple_threshold_decoding(thresholds)) {
		//	std::cout << "Successfully decoded!\n";
		//}
		//else {
		//	std::cout << "Failed to decode!\n";
		//}
		inputfile.close();
	}
}
*/

int main() {
	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row = 5, col = 10;
	std::uniform_int_distribution<> posgen(0, col - 1);
	std::uniform_int_distribution<> xgen(1, gf.q - 1);
	int DijkstraNumZero = 0, DijkstraNumSatisfied = 0, DijkstraNumFailed = 0;
	int SingleThresholdNumZero = 0, SingleThresholdNumSatisfied = 0, SingleThresholdNumFailed = 0;
	int MultThresholdNumZero = 0, MultThresholdNumSatisfied = 0, MultThresholdNumFailed = 0;
	std::vector<int> thresholds = { 0,1,2 };
	for (int i = 0; i < 1000; i++) {
		Matrix m(row, col, &gf);
		m.random();
		//std::cout << m;
		Flipper f(&m);
		std::vector<int> x(col,0);
		x[posgen(gen)] = xgen(gen);
		x[posgen(gen)] = xgen(gen);
		x[posgen(gen)] = xgen(gen);
		// Run Dijkstra
		f.assignX(x);
		f.Dijkstra_solve_extended_bit_flipping_consecutively();
		bool isZero = f.isZeroVector(), isSatisfied = f.isMatrixSatisfied();
		if (isZero) DijkstraNumZero++;
		if (isSatisfied) DijkstraNumSatisfied++;
		if (!isZero && !isSatisfied) DijkstraNumFailed++;
		std::cout << "Dijkstra: " << f.isZeroVector() << ' ' << f.isMatrixSatisfied() << ' ';

		// Run Single threshold
		f.assignX(x);
		f.solve_single_threshold_decoding(0);
		isZero = f.isZeroVector(), isSatisfied = f.isMatrixSatisfied();
		if (isZero) SingleThresholdNumZero++;
		if (isSatisfied) SingleThresholdNumSatisfied++;
		if (!isZero && !isSatisfied) SingleThresholdNumFailed++;
		std::cout << "SingleThreshold: " << f.isZeroVector() << ' ' << f.isMatrixSatisfied() << ' ';


		// Run Single threshold
		f.assignX(x);
		f.solve_multiple_threshold_decoding(thresholds);
		isZero = f.isZeroVector(), isSatisfied = f.isMatrixSatisfied();
		if (isZero) MultThresholdNumZero++;
		if (isSatisfied) MultThresholdNumSatisfied++;
		if (!isZero && !isSatisfied) MultThresholdNumFailed++;
		std::cout << "MultThreshold: " << f.isZeroVector() << ' ' << f.isMatrixSatisfied() << '\n';

	}
	std::cout << DijkstraNumZero << ' ' << DijkstraNumSatisfied << ' ' << DijkstraNumFailed << '\n';
	std::cout << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << '\n';
	std::cout << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
	/*
	Flipper f(&m);

	std::vector<int> x(col);
	for (int i = 0; i < col; i++) {
		inputfile >> x[i];
	}
	f.assignX(x);
	std::vector<int> b(row);
	for (int i = 0; i < row; i++) {
		inputfile >> b[i];
	}
	f.assignB(b);
	f.FloydWarshall();
	m.findSingleVar();
	std::cout << m;
	f.solve_extended_bit_flipping_consecutively();
	*/

}