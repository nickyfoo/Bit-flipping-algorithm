#pragma once

#include <galois++/element.h>
#include <set>

// A sparse m x n matrix
class Matrix {
private:
	int m = 0, n = 0;
	std::vector<std::vector<std::pair<int, Galois::Element>>> equationWiseEntries, variableWiseEntries;
	std::set<std::pair<int, int>> filledEntries;
	std::vector<int> indexOfAdjSingleVar;
	Galois::Field* gf;
public:

	Matrix(int row, int col, Galois::Field* gfp) {
		m = row;
		n = col;
		equationWiseEntries.resize(row);
		variableWiseEntries.resize(col);
		indexOfAdjSingleVar.assign(row, -1);
		gf = gfp;
		
	}

	void random();

	bool insert(int i, int j, Galois::Element x);
	Galois::Element solveEqn(int eqn, int var, std::vector<Galois::Element> x, Galois::Element d);
	bool isEqnSatisfied(int eqn, std::vector<Galois::Element> x, Galois::Element d);
	
	int getDegree(int i);
	//Remember to call findSingleVar after populating entries
	void findSingleVar();
	int getIndexOfAdjSingleVar(int i);

	int getNumRows();
	int getNumCols();
	Galois::Field* getField();
	std::vector<std::pair<int, Galois::Element>>* getAdjacentEquations(int var);
	std::vector<std::pair<int, Galois::Element>>* getAdjacentVariables(int eqn);

	friend std::ostream& operator<<(std::ostream& output, Matrix& matrix);

};