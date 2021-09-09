#include <iostream>
#include <fstream>
#include <galois++/field.h>
#include <galois++/element.h>
#include "matrix.h"
#include "flipper.h"

// For taking input
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
