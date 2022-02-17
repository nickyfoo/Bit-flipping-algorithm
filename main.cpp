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