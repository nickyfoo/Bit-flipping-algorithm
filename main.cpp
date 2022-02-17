#include <iostream>
#include <fstream>
#include <galois++/field.h>
#include <galois++/element.h>
#include "matrix.h"
#include "flipper.h"
#include <random>
#include <fstream>
#include <chrono>
#include <unordered_set>
#include <string>

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


// Single run trying to decode
/*
int main() {
	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row = 100, col = 200;
	std::uniform_int_distribution<> posgen(0, col - 1);
	std::uniform_int_distribution<> xgen(1, gf.q - 1);
	int DijkstraNumZero = 0, DijkstraNumFailed = 0;
	int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
	int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
	std::vector<int> thresholds = { 0,1,2 };
	for (int i = 0; i < 1000; i++) {
		Matrix m(row, col, &gf);
		m.random();
		//std::cout << m;
		Flipper f(&m);
		std::vector<int> x(col,0);
		x[posgen(gen)] = xgen(gen);
		x[posgen(gen)] = xgen(gen);
		//x[posgen(gen)] = xgen(gen);
		// Run Dijkstra
		f.assignX(x);
		f.Dijkstra_solve_extended_bit_flipping_consecutively();
		bool isZero = f.isZeroVector();
		if (isZero) DijkstraNumZero++;
		else DijkstraNumFailed++;
		std::cout << "Dijkstra: " << f.isZeroVector() << ' ';

		// Run Single threshold
		f.assignX(x);
		f.solve_single_threshold_decoding(0);
		isZero = f.isZeroVector();
		if (isZero) SingleThresholdNumZero++;
		else SingleThresholdNumFailed++;
		std::cout << "SingleThreshold: " << f.isZeroVector() << ' ';


		// Run Single threshold
		f.assignX(x);
		f.solve_multiple_threshold_decoding(thresholds);
		isZero = f.isZeroVector();
		if (isZero) MultThresholdNumZero++;
		else MultThresholdNumFailed++;
		std::cout << "MultThreshold: " << f.isZeroVector() << '\n';

	}
	std::cout << DijkstraNumZero << ' ' << DijkstraNumFailed << '\n';
	std::cout << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << '\n';
	std::cout << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
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
	
}
*/

// Running large number of decodings to get statistics
/*
std::ostream& printStats(std::ostream& output, std::vector<std::pair<double, int>>& v) {
	int n = v.size();
	double sum = 0;
	for (auto& [x,i] : v) sum += x;
	output << v[0].first << "," << v[0].second << ' ' << sum / n << ' ' << v[n - 1].first << "," << v[n - 1].second << '\n';
	return output;
}

int main() {
	std::ofstream logfile;
	logfile.open("log.txt");

	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row = 100, col = 200;
	std::uniform_int_distribution<> posgen(0, col - 1);
	std::uniform_int_distribution<> xgen(1, gf.q - 1);
	std::vector<int> thresholds = { 0,1,2 };

	//statistics:
	double DijkstraSum = 0, SingleThresholdSum = 0, MultThresholdSum = 0;
	std::vector<std::pair<double,int>> DijkstraResult, SingleThresholdResult, MultThresholdResult;
	int numMatrices = 1000;
	for (int i = 0; i < numMatrices; i++) {
		int DijkstraNumZero = 0, DijkstraNumFailed = 0, DijkstraNumPaths = 0;
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		Matrix m(row, col, &gf);
		m.random();
		logfile << "Index " << i << '\n';
		m.printEntries(logfile);

		Flipper f(&m);
		for (int j = 0; j < 100; j++) {

			std::vector<int> x(col, 0);
			std::vector<int> pos, val;
			logfile << "Input vector: ";
			for (int i = 0; i < 2; i++) {
				pos.push_back(posgen(gen));
				val.push_back(xgen(gen));
				logfile << pos[i] << ": " << val[i] << ' ';
				x[pos[i]] = val[i];
			}
			logfile << '\n';

			// Run Dijkstra
			f.assignX(x);
			f.Dijkstra_solve_extended_bit_flipping_consecutively();
			bool isZero = f.isZeroVector();
			if (isZero) DijkstraNumZero++;
			else DijkstraNumFailed++;
			DijkstraNumPaths += f.Dijkstra_getNumPaths();
			logfile << "Dijkstra: " << f.isZeroVector() << ' ' << f.Dijkstra_getNumPaths() << ' ';

			// Run Single threshold
			f.assignX(x);
			f.solve_single_threshold_decoding(0);
			isZero = f.isZeroVector();
			if (isZero) SingleThresholdNumZero++;
			else SingleThresholdNumFailed++;
			logfile << "SingleThreshold: " << f.isZeroVector() << ' ';


			// Run Single threshold
			f.assignX(x);
			f.solve_multiple_threshold_decoding(thresholds);
			isZero = f.isZeroVector();
			if (isZero) MultThresholdNumZero++;
			else MultThresholdNumFailed++;
			logfile << "MultThreshold: " << f.isZeroVector() << '\n';

			logfile << DijkstraNumZero << ' ' << DijkstraNumFailed << ' ' << DijkstraNumPaths << '\n';
			logfile << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << '\n';
			logfile << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';

		}
		std::cout << "Index: " << i << '\n';
		std::cout << DijkstraNumZero << ' ' << DijkstraNumFailed << ' ' << DijkstraNumPaths << '\n';
		std::cout << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << '\n';
		std::cout << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		DijkstraResult.push_back({ DijkstraNumZero, i });
		SingleThresholdResult.push_back({ SingleThresholdNumZero, i });
		MultThresholdResult.push_back({ MultThresholdNumZero, i });
		DijkstraSum += DijkstraNumZero;
		SingleThresholdSum += SingleThresholdNumZero;
		MultThresholdSum += MultThresholdNumZero;
	}
	std::sort(DijkstraResult.begin(), DijkstraResult.end());
	std::sort(SingleThresholdResult.begin(), SingleThresholdResult.end());
	std::sort(MultThresholdResult.begin(), MultThresholdResult.end());
	printStats(logfile, DijkstraResult);
	printStats(logfile, SingleThresholdResult);
	printStats(logfile, MultThresholdResult);
	printStats(std::cout, DijkstraResult);
	printStats(std::cout, SingleThresholdResult);
	printStats(std::cout, MultThresholdResult);
}
*/

// Generating matrices with Gallager's method
/*
int main() {
	std::ofstream logfile;
	logfile.open("density 0.02 errors 4.txt");

	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row_limit = 1000, col_limit = 1000;
	int row_inc = 100, col_inc = 100;
	int nummatrices = 10;
	int numtrials = 100; 
	int numErrors = 4;
	double densityThreshold = 0.02;
	std::vector<int> thresholds;
	for (int i = 0; i <= numErrors; i++) {
		thresholds.push_back(i);
	}
	logfile << row_limit << ' ' << col_limit << ' ' << row_inc << ' ' << col_inc << ' ' << nummatrices << ' ' << numtrials << '\n';
	for (int row = 100; row <= row_limit; row += row_inc) {
		for (int col = 100; col <= col_limit; col += col_inc) {

			// Increasing density of matrices
			for (int rowWt = 1; rowWt <= col; rowWt++) {
				if (double(rowWt) / col != densityThreshold) continue;
				if (col % rowWt == 0) {
					int numRowsPerLayer = col / rowWt;
					if (row % numRowsPerLayer == 0) {

						std::uniform_int_distribution<> posgen(0, col - 1);
						std::uniform_int_distribution<> xgen(1, gf.q - 1);
						int colWt = row / numRowsPerLayer;
						std::cout << "row: " << row << " col: " << col << " rowWt: " << rowWt << " colWt: " << colWt << " density: " << double(rowWt)/col << '\n';
						logfile << row << ' ' << col << ' ' << rowWt << ' ' << colWt << ' ' << double(rowWt) / col << ' ';

						bool isZero, isSatisfied;
						int SingleThresholdNumZero = 0, SingleThresholdNumSatisfied = 0, SingleThresholdNumFailed = 0;
						int MultThresholdNumZero = 0, MultThresholdNumSatisfied = 0, MultThresholdNumFailed = 0;

						for (int k = 0; k < nummatrices; k++) {

							Matrix m(0, 0, &gf);
							for (int i = 0; i < colWt; i++) {
								m.addLayer(rowWt, col);
							}
							//std::cout << m << '\n';
							Flipper f(&m);

							auto starttime = std::chrono::steady_clock::now();
							for (int j = 0; j < numtrials; j++) {
								std::vector<int> x(col, 0);
								std::vector<int> pos, val;
								for (int i = 0; i < numErrors; i++) {
									pos.push_back(posgen(gen));
									val.push_back(xgen(gen));
									x[pos[i]] = val[i];
								}
								// Run Single threshold
								f.assignX(x);
								f.solve_single_threshold_decoding(0);
								isZero = f.isZeroVector();
								isSatisfied = f.isMatrixSatisfied();
								if (isSatisfied) SingleThresholdNumSatisfied++;
								if (isZero) SingleThresholdNumZero++;
								else SingleThresholdNumFailed++;

								// Run Mult threshold
								f.assignX(x);
								f.solve_multiple_threshold_decoding(thresholds);
								isZero = f.isZeroVector();
								isSatisfied = f.isMatrixSatisfied();
								if (isSatisfied) MultThresholdNumSatisfied++;
								if (isZero) MultThresholdNumZero++;
								else MultThresholdNumFailed++;
							}

							auto endtime = std::chrono::steady_clock::now();
							std::chrono::duration<double> timeelapsed = endtime - starttime;
							std::cout << k << ' ' << timeelapsed.count() << '\n';

						}

						std::cout << "Single: " << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << '\n';
						std::cout << "Mult: " << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
						logfile << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << ' ';
						logfile << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
					}

				}
			}
		}
	}
	//std::uniform_int_distribution<> posgen(0, col - 1);
	//std::uniform_int_distribution<> xgen(1, gf.q - 1);
	//std::vector<int> thresholds = { 0,1,2 };
	//Matrix m(0, 0, &gf);
	//m.addLayer(4, 16);
	//m.addLayer(4, 16);
	//std::cout << m << '\n';
	//Flipper f(&m);
}
*/

/*
// Generating matrices with up to colWt column weight
int main() {
	std::ofstream logfile;
	logfile.open("increasing colwt errors 3.txt");

	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row_limit = 1000, col_limit = 1000;
	int row_inc = 100, col_inc = 100;
	int nummatrices = 1;
	int numtrials = 100;
	int numErrors = 3;
	std::vector<int> thresholds;
	for (int i = 0; i <= numErrors; i++) {
		thresholds.push_back(i);
	}
	logfile << row_limit << ' ' << col_limit << ' ' << row_inc << ' ' << col_inc << ' ' << nummatrices << ' ' << numtrials << '\n';
	for (int row = 100; row <= row_limit; row += row_inc) {
		for (int col = 100; col <= col_limit; col += col_inc) {

			// Increasing colWt of the matrices
			for (int colWt = 1; colWt < row/2; colWt++) {
				std::uniform_int_distribution<> posgen(0, col - 1);
				std::uniform_int_distribution<> xgen(1, gf.q - 1);
				std::cout << "row: " << row << " col: " << col << " colWt: " << colWt << '\n';
				logfile << row << ' ' << col << ' ' << colWt << ' ';

				bool isZero, isSatisfied;
				int SingleThresholdNumZero = 0, SingleThresholdNumSatisfied = 0, SingleThresholdNumFailed = 0;
				int MultThresholdNumZero = 0, MultThresholdNumSatisfied = 0, MultThresholdNumFailed = 0;

				for (int k = 0; k < nummatrices; k++) {
					Matrix m(row, col, &gf);
					m.random_colwt(colWt);
					Flipper f(&m);

					auto starttime = std::chrono::steady_clock::now();
					for (int j = 0; j < numtrials; j++) {
						std::vector<int> x(col, 0);
						std::vector<int> pos, val;
						for (int i = 0; i < numErrors; i++) {
							pos.push_back(posgen(gen));
							val.push_back(xgen(gen));
							x[pos[i]] = val[i];
						}
						// Run Single threshold
						f.assignX(x);
						f.solve_single_threshold_decoding(0);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) SingleThresholdNumSatisfied++;
						if (isZero) SingleThresholdNumZero++;
						else SingleThresholdNumFailed++;

						// Run Mult threshold
						f.assignX(x);
						f.solve_multiple_threshold_decoding(thresholds);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) MultThresholdNumSatisfied++;
						if (isZero) MultThresholdNumZero++;
						else MultThresholdNumFailed++;
					}

					auto endtime = std::chrono::steady_clock::now();
					std::chrono::duration<double> timeelapsed = endtime - starttime;
					std::cout << k << ' ' << timeelapsed.count() << '\n';
				}

				std::cout << "Single: " << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << '\n';
				std::cout << "Mult: " << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				logfile << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << ' ';
				logfile << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				if (SingleThresholdNumZero == numtrials) break;
			}
		}
	}
}
*/

/*
// Generating matrices with fixed size, fixed colwt.
int main() {
	std::ofstream logfile;
	std::string filename = "test 100x100 errors 10.txt";
	std::cout << filename << '\n';
	logfile.open(filename);
	

	int fieldOrder = 5;
	Galois::Field gf(fieldOrder);
	std::cout << gf << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row_limit = 100, col_limit = 100;
	int row_inc = 100, col_inc = 100;
	int nummatrices = 10;
	int numtrials = 100;
	int numErrors = 10;
	std::vector<int> thresholds;
	for (int i = 0; i <= numErrors; i++) {
		thresholds.push_back(i);
	}
	logfile << row_limit << ' ' << col_limit << ' ' << row_inc << ' ' << col_inc << ' ' << nummatrices << ' ' << numtrials << '\n';
	for (int row = row_limit; row == row_limit; row += row_inc) {
		for (int col = col_limit; col == col_limit; col += col_inc) {

			// Increasing colWt of the matrices
			for (int colWt = 1; colWt < row / 2; colWt++) {
				std::uniform_int_distribution<> posgen(0, col - 1);
				std::uniform_int_distribution<> xgen(1, gf.q - 1);
				std::cout << "row: " << row << " col: " << col << " colWt: " << colWt << '\n';
				logfile << row << ' ' << col << ' ' << colWt << ' ';

				bool isZero, isSatisfied;
				int SingleThresholdNumZero = 0, SingleThresholdNumSatisfied = 0, SingleThresholdNumFailed = 0;
				int MultThresholdNumZero = 0, MultThresholdNumSatisfied = 0, MultThresholdNumFailed = 0;

				for (int k = 0; k < nummatrices; k++) {
					Matrix m(row, col, &gf);
					m.random_colwt(colWt);
					Flipper f(&m);

					auto starttime = std::chrono::steady_clock::now();
					for (int j = 0; j < numtrials; j++) {
						std::vector<int> x(col, 0);
						std::vector<int> pos, val;
						std::unordered_set<int> pos_filled;
						while (pos_filled.size() < numErrors) {
							int newpos = posgen(gen);
							if (pos_filled.find(newpos) == pos_filled.end()) {
								pos_filled.insert(newpos);
								pos.push_back(posgen(gen));
								val.push_back(xgen(gen));
							}
						}
						for (int i = 0; i < numErrors; i++) {
							x[pos[i]] = val[i];
						}
						// Run Single threshold
						f.assignX(x);
						f.solve_single_threshold_decoding(0);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) SingleThresholdNumSatisfied++;
						if (isZero) SingleThresholdNumZero++;
						else SingleThresholdNumFailed++;

						// Run Mult threshold
						f.assignX(x);
						f.solve_multiple_threshold_decoding(thresholds);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) MultThresholdNumSatisfied++;
						if (isZero) MultThresholdNumZero++;
						else MultThresholdNumFailed++;
					}

					auto endtime = std::chrono::steady_clock::now();
					std::chrono::duration<double> timeelapsed = endtime - starttime;
					std::cout << k << ' ' << timeelapsed.count() << '\n';
				}

				std::cout << "Single: " << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << '\n';
				std::cout << "Mult: " << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				logfile << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << ' ';
				logfile << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				// Stop if fully decoded, to save on computations.
				//if (SingleThresholdNumZero == numtrials*nummatrices) break;
			}
		}
	}
}
*/

/*
// Generating matrices by finite affine plane, X = F_q X F_q, D = lines.
int main() {
	std::ofstream logfile;
	std::string filename = "test 100x100 errors 10.txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	int q = 5;
	int fieldOrder = q*q*q;
	Galois::Field F_q(q);
	std::cout << F_q << std::endl;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	int row_limit = 100, col_limit = 100;
	int row_inc = 100, col_inc = 100;
	int nummatrices = 10;
	int numtrials = 100;
	int numErrors = 10;
	std::vector<int> thresholds;
	for (int i = 0; i <= numErrors; i++) {
		thresholds.push_back(i);
	}

	// Pick numrows <= q*q
	int numrows = 10;

	Matrix affine(q*q, q * q, &F_q, &F_q);
	affine.affine_plane(q*q);
	std::cout << affine << '\n';

	Matrix diff_set(q * q, q * q, &F_q, &F_q);
	diff_set.difference_set(q * q);
	std::cout << diff_set << '\n';
	/*
	Galois::Field gf(fieldOrder);
	Matrix proj_geom(q * q + q + 1, q * q + q + 1, &F_q, &gf);
	proj_geom.projective_geometry(q * q + q + 1);
	std::cout << proj_geom;
	*/
	
	/*
	logfile << row_limit << ' ' << col_limit << ' ' << row_inc << ' ' << col_inc << ' ' << nummatrices << ' ' << numtrials << '\n';
	for (int row = row_limit; row == row_limit; row += row_inc) {
		for (int col = col_limit; col == col_limit; col += col_inc) {

			// Increasing colWt of the matrices
			for (int colWt = 1; colWt < row / 2; colWt++) {
				std::uniform_int_distribution<> posgen(0, col - 1);
				std::uniform_int_distribution<> xgen(1, gf.q - 1);
				std::cout << "row: " << row << " col: " << col << " colWt: " << colWt << '\n';
				logfile << row << ' ' << col << ' ' << colWt << ' ';

				bool isZero, isSatisfied;
				int SingleThresholdNumZero = 0, SingleThresholdNumSatisfied = 0, SingleThresholdNumFailed = 0;
				int MultThresholdNumZero = 0, MultThresholdNumSatisfied = 0, MultThresholdNumFailed = 0;

				for (int k = 0; k < nummatrices; k++) {
					Matrix m(row, col, &gf);
					m.random_colwt(colWt);
					Flipper f(&m);

					auto starttime = std::chrono::steady_clock::now();
					for (int j = 0; j < numtrials; j++) {
						std::vector<int> x(col, 0);
						std::vector<int> pos, val;
						std::unordered_set<int> pos_filled;
						while (pos_filled.size() < numErrors) {
							int newpos = posgen(gen);
							if (pos_filled.find(newpos) == pos_filled.end()) {
								pos_filled.insert(newpos);
								pos.push_back(posgen(gen));
								val.push_back(xgen(gen));
							}
						}
						for (int i = 0; i < numErrors; i++) {
							x[pos[i]] = val[i];
						}
						// Run Single threshold
						f.assignX(x);
						f.solve_single_threshold_decoding(0);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) SingleThresholdNumSatisfied++;
						if (isZero) SingleThresholdNumZero++;
						else SingleThresholdNumFailed++;

						// Run Mult threshold
						f.assignX(x);
						f.solve_multiple_threshold_decoding(thresholds);
						isZero = f.isZeroVector();
						isSatisfied = f.isMatrixSatisfied();
						if (isSatisfied) MultThresholdNumSatisfied++;
						if (isZero) MultThresholdNumZero++;
						else MultThresholdNumFailed++;
					}

					auto endtime = std::chrono::steady_clock::now();
					std::chrono::duration<double> timeelapsed = endtime - starttime;
					std::cout << k << ' ' << timeelapsed.count() << '\n';
				}

				std::cout << "Single: " << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << '\n';
				std::cout << "Mult: " << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				logfile << SingleThresholdNumZero << ' ' << SingleThresholdNumSatisfied << ' ' << SingleThresholdNumFailed << ' ';
				logfile << MultThresholdNumZero << ' ' << MultThresholdNumSatisfied << ' ' << MultThresholdNumFailed << '\n';
				// Stop if fully decoded, to save on computations.
				//if (SingleThresholdNumZero == numtrials*nummatrices) break;
			}
		}
	}

}
*/

std::tuple<int,int,int,int> run_trials(int numTrials, int numErrors, Matrix *m) {
	int rows = m->getNumRows();
	int cols = m->getNumCols();

	std::vector<int> thresholds;
	for (int i = 0; i <= numErrors; i++) {
		thresholds.push_back(i);
	}
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::uniform_int_distribution<> posgen(0, cols - 1);
	std::uniform_int_distribution<> xgen(1, m->getField()->q - 1);

	bool isZero;
	int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
	int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;

	Flipper f(m);

	auto starttime = std::chrono::steady_clock::now();
	for (int j = 0; j < numTrials; j++) {
		std::vector<int> x(cols, 0);
		std::vector<int> pos, val;
		std::unordered_set<int> pos_filled;
		while (pos_filled.size() < numErrors) {
			int newpos = posgen(gen);
			if (pos_filled.find(newpos) == pos_filled.end()) {
				pos_filled.insert(newpos);
				pos.push_back(posgen(gen));
				val.push_back(xgen(gen));
			}
		}
		for (int i = 0; i < numErrors; i++) {
			x[pos[i]] = val[i];
		}

		// Run Single threshold
		f.assignX(x);
		f.solve_single_threshold_decoding(0);
		isZero = f.isZeroVector();

		if (isZero) SingleThresholdNumZero++;
		else SingleThresholdNumFailed++;

		// Run Mult threshold
		f.assignX(x);
		f.solve_multiple_threshold_decoding(thresholds);
		isZero = f.isZeroVector();

		if (isZero) MultThresholdNumZero++;
		else MultThresholdNumFailed++;
	}

	auto endtime = std::chrono::steady_clock::now();
	std::chrono::duration<double> timeelapsed = endtime - starttime;
	std::cout << timeelapsed.count() << '\n';

	return { SingleThresholdNumZero, SingleThresholdNumFailed, MultThresholdNumZero, MultThresholdNumFailed };
}

// Make the matrix here
void test_projective_geom(int lowErrors, int highErrors, int numMatrices, int numTrials, int numRows, int q) {
	int fieldOrder = q * q * q;
	int numCols = q * q + q + 1;

	std::ofstream logfile;
	const std::string filename = "projective geom " + std::to_string(numRows) + "x" + std::to_string(numCols) + " " 
		+ std::to_string(lowErrors) + " " + std::to_string(highErrors) + ".txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	Galois::Field F_q(q);
	Galois::Field gf(fieldOrder);
	std::cout << "Fields set up\n";

	for (int numErrors = lowErrors; numErrors <= highErrors; numErrors++) {
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		for (int i = 0; i < numMatrices; i++) {
			Matrix m(numRows, q * q + q + 1, &F_q, &gf);
			m.projective_geometry(numRows);
			auto[singleSuccess, singleFail, multSuccess, multFail] = run_trials(numTrials, numErrors, &m);
			SingleThresholdNumZero += singleSuccess;
			SingleThresholdNumFailed += singleFail;
			MultThresholdNumZero += multSuccess;
			MultThresholdNumFailed += multFail;
			std::cout << i << ' ' << singleSuccess << ' ' << singleFail << ' ' << multSuccess << ' ' << multFail << '\n';
		}
		std::cout << numErrors << " Errors: " << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile << numErrors << ' ' << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile.flush();
		if (MultThresholdNumZero == 0) break;
	}
	logfile.close();
}

// Make the matrix here
// has rowWt q
void test_affine_plane(int lowErrors, int highErrors, int numMatrices, int numTrials, int numRows, int q) {
	int numCols = q * q;

	std::ofstream logfile;
	const std::string filename = "affine plane " + std::to_string(numRows) + "x" + std::to_string(numCols) + " "
		+ std::to_string(lowErrors) + " " + std::to_string(highErrors) + ".txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	Galois::Field F_q(q);
	std::cout << "Fields set up\n";

	for (int numErrors = lowErrors; numErrors <= highErrors; numErrors++) {
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		for (int i = 0; i < numMatrices; i++) {
			Matrix m(numRows, q * q, &F_q, &F_q);
			m.affine_plane(numRows);
			auto [singleSuccess, singleFail, multSuccess, multFail] = run_trials(numTrials, numErrors, &m);
			SingleThresholdNumZero += singleSuccess;
			SingleThresholdNumFailed += singleFail;
			MultThresholdNumZero += multSuccess;
			MultThresholdNumFailed += multFail;
			std::cout << i << ' ' << singleSuccess << ' ' << singleFail << ' ' << multSuccess << ' ' << multFail << '\n';
		}
		std::cout << numErrors << " Errors: " << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile << numErrors << ' ' << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile.flush();
		if (MultThresholdNumZero == 0) break;
	}
	logfile.close();
}


// Make the matrix here
void test_difference_set(int lowErrors, int highErrors, int numMatrices, int numTrials, int numRows, int q) {
	int numCols = q * q;

	std::ofstream logfile;
	const std::string filename = "difference set " + std::to_string(numRows) + "x" + std::to_string(numCols) + " "
		+ std::to_string(lowErrors) + " " + std::to_string(highErrors) + ".txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	Galois::Field F_q(q);
	std::cout << "Fields set up\n";

	for (int numErrors = lowErrors; numErrors <= highErrors; numErrors++) {
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		for (int i = 0; i < numMatrices; i++) {
			Matrix m(numRows, q * q, &F_q, &F_q);
			m.difference_set(numRows);
			auto [singleSuccess, singleFail, multSuccess, multFail] = run_trials(numTrials, numErrors, &m);
			SingleThresholdNumZero += singleSuccess;
			SingleThresholdNumFailed += singleFail;
			MultThresholdNumZero += multSuccess;
			MultThresholdNumFailed += multFail;
			std::cout << i << ' ' << singleSuccess << ' ' << singleFail << ' ' << multSuccess << ' ' << multFail << '\n';
		}
		std::cout << numErrors << " Errors: " << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile << numErrors << ' ' << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile.flush();
		if (MultThresholdNumZero == 0) break;
	}
	logfile.close();
}

// Make the matrix here
void test_random_col_wt(int lowErrors, int highErrors, int numMatrices, int numTrials, int numRows, int q, int colWt) {
	int numCols = q * q;

	std::ofstream logfile;
	const std::string filename = "col wt " + std::to_string(colWt) + " " + std::to_string(numRows) + "x" + std::to_string(numCols) + " "
		+ std::to_string(lowErrors) + " " + std::to_string(highErrors) + ".txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	Galois::Field F_q(q);
	std::cout << "Fields set up\n";

	for (int numErrors = lowErrors; numErrors <= highErrors; numErrors++) {
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		for (int i = 0; i < numMatrices; i++) {
			Matrix m(numRows, q * q, &F_q, &F_q);
			m.random_colwt(colWt);
			auto [singleSuccess, singleFail, multSuccess, multFail] = run_trials(numTrials, numErrors, &m);
			SingleThresholdNumZero += singleSuccess;
			SingleThresholdNumFailed += singleFail;
			MultThresholdNumZero += multSuccess;
			MultThresholdNumFailed += multFail;
			std::cout << i << ' ' << singleSuccess << ' ' << singleFail << ' ' << multSuccess << ' ' << multFail << '\n';
		}
		std::cout << numErrors << " Errors: " << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile << numErrors << ' ' << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile.flush();
		if (MultThresholdNumZero == 0) break;
	}
	logfile.close();
}

// Make the matrix here
void test_fixed_row_wt(int lowErrors, int highErrors, int numMatrices, int numTrials, int numRows, int q, int rowWt) {
	int numCols = q * q;

	std::ofstream logfile;
	const std::string filename = "fixed row wt " + std::to_string(rowWt) + " " + std::to_string(numRows) + "x" + std::to_string(numCols) + " "
		+ std::to_string(lowErrors) + " " + std::to_string(highErrors) + ".txt";
	std::cout << filename << '\n';
	logfile.open(filename);

	Galois::Field F_q(q);
	std::cout << "Fields set up\n";

	for (int numErrors = lowErrors; numErrors <= highErrors; numErrors++) {
		int SingleThresholdNumZero = 0, SingleThresholdNumFailed = 0;
		int MultThresholdNumZero = 0, MultThresholdNumFailed = 0;
		for (int i = 0; i < numMatrices; i++) {
			Matrix m(numRows, q * q, &F_q, &F_q);
			m.fixed_row_wt(numRows, q * q, rowWt);
			auto [singleSuccess, singleFail, multSuccess, multFail] = run_trials(numTrials, numErrors, &m);
			SingleThresholdNumZero += singleSuccess;
			SingleThresholdNumFailed += singleFail;
			MultThresholdNumZero += multSuccess;
			MultThresholdNumFailed += multFail;
			std::cout << i << ' ' << singleSuccess << ' ' << singleFail << ' ' << multSuccess << ' ' << multFail << '\n';
		}
		std::cout << numErrors << " Errors: " << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile << numErrors << ' ' << SingleThresholdNumZero << ' ' << SingleThresholdNumFailed << ' ' << MultThresholdNumZero << ' ' << MultThresholdNumFailed << '\n';
		logfile.flush();
		if (MultThresholdNumZero == 0) break;
	}
	logfile.close();
}

int main() {
	int lowErrors = 1;
	int highErrors = 100;
	int totalErrors = 10;
	int numMatrices = 10;
	int numTrials = 100;
	int q = 25;
	// Choose numRows <= q*q
	int numRows = 300;

	// Testing for row wt of projective geometry
	//test_fixed_row_wt(lowErrors, highErrors, numMatrices, numTrials, numRows, q, q + 1);

	//test_projective_geom(lowErrors, highErrors, numMatrices, numTrials, numRows, q);
	test_affine_plane(lowErrors, highErrors, numMatrices, numTrials, numRows, q);
	test_difference_set(54, 100, numMatrices, numTrials, numRows, q);
	//test_random_col_wt(lowErrors, highErrors, numMatrices, numTrials, numRows, q, 1);
}