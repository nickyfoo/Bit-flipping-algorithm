#include <galois++/element.h>
#include "matrix.h"
#include <map>
#include <random>

bool Matrix::insert(int i, int j, Galois::Element x) {
	if (0 <= i && i < m && 0 <= j && j < n) {
		if (filledEntries.find({ i,j }) != filledEntries.end()) {
			return false;
		}
		equationWiseEntries[i].push_back({ j,x });
		variableWiseEntries[j].push_back({ i,x });
		filledEntries.insert({ i,j });
		return true;
	}
	return false;
}

void Matrix::random() {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> colWeight(1, 2);
	std::uniform_int_distribution<> randRow(0, m-1);
	std::uniform_int_distribution<> entry(1, (gf->q)-1);
	for (int j = 0; j < n; j++) {
		int columnWeight = colWeight(gen);
		for (int k = 0; k < columnWeight; k++) {
			int i = randRow(gen);
			int entryVal = entry(gen);
			insert(i, j, Galois::Element(gf, entryVal));
		}
	}
}


/**
	\brief solves cx=d for parity check equation eqn, and variable var. 

	\param eqn   Parity check equation index (row of the parity check matrix)
	\param var	 variable node that we are solving for
	\param x	 vector of all x_i, input variables
*/
Galois::Element Matrix::solveEqn(int eqn, int var, std::vector<Galois::Element> x, Galois::Element d){
	Galois::Element ans = d;
	Galois::Element coeff(gf, 0);
	for (auto& [j, entry] : equationWiseEntries[eqn]) {
		if (j != var) {
			ans = ans - (entry * x[j]);
		}
		else {
			coeff = entry;
		}
	}
	ans = ans / coeff;
	return ans;
}

bool Matrix::isEqnSatisfied(int eqn, std::vector<Galois::Element> x, Galois::Element d) {
	Galois::Element ans = d;
	for (auto& [j, entry] : equationWiseEntries[eqn]) {
		ans = ans - (entry * x[j]);
	}
	//todo: check if this equality is overloaded.
	return ans == Galois::Element(gf, 0);
}

int Matrix::getDegree(int i) {
	return variableWiseEntries[i].size();
}


/**
	\brief Populates isAdjToSingleVar with Eqn nodes that are adj to a single variable node for path finding.
*/
void Matrix::findSingleVar() {
	//Iterate over variables
	for (int i = 0; i < n; i++) {
		if (getDegree(i) == 1) {
			for (auto& [j, entry] : variableWiseEntries[i]) {
				//Equation node j is adj to variable node i with degree 1
				indexOfAdjSingleVar[j] = i;
			}
		}
	}
}

int Matrix::getIndexOfAdjSingleVar(int i) {
	return indexOfAdjSingleVar[i];
}

int Matrix::getNumRows() {
	return m;
}
int Matrix::getNumCols() {
	return n;
}

Galois::Field* Matrix::getField() {
	return gf;
}


std::vector<std::pair<int, Galois::Element>>* Matrix::getAdjacentEquations(int var) {
	return &variableWiseEntries[var];
}

std::vector<std::pair<int, Galois::Element>>* Matrix::getAdjacentVariables(int eqn) {
	return &equationWiseEntries[eqn];
}

std::ostream& operator<<(std::ostream& output, Matrix&matrix) {
	for (int i = 0; i < matrix.m; i++) {
		int rowindex = 0;
		if (matrix.equationWiseEntries[i].size() == 0) {
			for (int j = 0; j < matrix.n; j++) {
				output << "0\t";
			}
		}
		else {
			auto [jj, entry] = matrix.equationWiseEntries[i][rowindex];
			for (int j = 0; j < matrix.n; j++) {
				if (j == jj) {
					output << entry.value() << '\t';
					rowindex++;
					if (rowindex == matrix.equationWiseEntries[i].size()) {
						rowindex--;
					}
					jj = matrix.equationWiseEntries[i][rowindex].first;
					entry = matrix.equationWiseEntries[i][rowindex].second;
				}
				else {
					output << "0\t";
				}
			}
		}
		output << '\n';
	}
	return output;
}

