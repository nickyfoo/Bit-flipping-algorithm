#pragma once

#include <galois++/element.h>
#include <set>
#include <map>

// A sparse m x n matrix
class Matrix {
private:
	int m = 0, n = 0;
	std::vector<std::vector<std::pair<int, Galois::Element>>> equationWiseEntries, variableWiseEntries;
	std::set<std::pair<int, int>> filledEntries;
	std::vector<int> indexOfAdjSingleVar;
	Galois::Field* gf;

	Galois::Field* construction_helper;
public:

	Matrix(int row, int col, Galois::Field* gfp, Galois::Field* ch) {
		m = row;
		n = col;
		equationWiseEntries.resize(row);
		variableWiseEntries.resize(col);
		indexOfAdjSingleVar.assign(row, -1);
		gf = gfp;
		construction_helper = ch;
	}




	void random();
	void random_colwt(int colwt);
	void addLayer(int rowWt, int numCols);

	// Fixed row wt
	void fixed_row_wt(int numrows, int numCols, int rowWt);

	// Finite affine plane generation, X = F_q X F_q
	void affine_plane(int numrows);

	// map from smallest elt of line (representative) to index.
	std::map<int, int> line_rep_to_index;
	std::map<int, std::set<int>> lines;
	void init_line_tables();

	int find_primitive();
	// Returns a vector of m^2 + m + 1 points (lines)
	std::vector<int> get_line();
	Galois::Element field_trace(Galois::Element);
	void evaluate_trace();
	std::set<int> trace_range;
	std::set<int> get_coset(std::vector<int> line, int offset);
	void projective_geometry(int numrows);

	void difference_set(int numrows);
	std::set<std::pair<int, int>> get_D();

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
  std::ostream& printEntries(std::ostream& output);

};