#include <galois++/element.h>
#include "matrix.h"
#include <map>
#include <random>
#include <algorithm>
#include <numeric>

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

void Matrix::random_colwt(int columnWeight) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> randRow(0, m - 1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);
	for (int j = 0; j < n; j++) {
		std::set<int> posfilled;
		while (posfilled.size() < columnWeight) {
			int i = randRow(gen);
			while (posfilled.find(i) != posfilled.end()) {
				i = randRow(gen);
			}
			posfilled.insert(i);
			int entryVal = entry(gen);
			insert(i, j, Galois::Element(gf, entryVal));
		}
	}
}

// numCols must be a multiple of rowWt
void Matrix::addLayer(int rowWt, int numCols) {
	int numAddedRows = numCols / rowWt;
	int prevSize = m;
	m += numCols / rowWt;
	n = numCols;
	equationWiseEntries.resize(m);
	variableWiseEntries.resize(n);

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> randPos(0,n-1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);
	
	std::vector<int> colPermutation(n);
	for (int i = 0; i < n; i++) {
		colPermutation[i] = i;
	}
	for (int i = n - 1; i >= 0; i--) {
		std::swap(colPermutation[i], colPermutation[randPos(gen)]);
	}

	for (int i = 0; i < numAddedRows; i++) {
		for (int j = 0; j < rowWt; j++) {
			int entryVal = entry(gen);
			insert(prevSize + i, colPermutation[i * rowWt + j], Galois::Element(gf, entryVal));
		}
		std::sort(equationWiseEntries[prevSize + i].begin(), equationWiseEntries[prevSize + i].end(), [](auto& a, auto& b) {
			return a.first < b.first;
			});
	}

}

// Implement this by choosing numrows out of all q^2 choose q lines.
// Gives a numrows x q^2 matrix
void Matrix::fixed_row_wt(int numrows, int numCols, int rowWt){
	int q = gf->q;
	int n = numCols;
	std::string starting_string = std::string(numCols - rowWt, '0') + std::string(rowWt, '1');
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> randPos(0, n - 1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);

	std::set<std::string> encountered;
	while (encountered.size() < numrows) {

		// gen a random permutation
		std::string nextstring = starting_string;
		for (int j = n - 1; j >= 0; j--) {
			std::swap(nextstring[j], nextstring[randPos(gen)]);
		}


		if (encountered.find(nextstring) == encountered.end()) {
			for (int i = 0; i < n; i++) {
				if (nextstring[i] == '1') {
					int entryVal = entry(gen);
					insert(encountered.size(), i, Galois::Element(gf, entryVal));
				}
			}
			encountered.insert(nextstring);
		}
	}
}

// Implement this by choosing numrows out of all q^2 choose q lines.
// Gives a numrows x q^2 matrix
// RowWt = q
void Matrix::affine_plane(int numrows) {
	int q = gf->q;
	int n = q * q;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> startPos(0, n - 1);
	std::uniform_int_distribution<> dGen(0, q - 1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);

	// Each line is determined by two points.
	std::set<std::pair<int,int>> encountered;
	while(encountered.size()<numrows){

		int start = startPos(gen);
		int i = start / q, j = start % q;
		int di = dGen(gen), dj = dGen(gen);
		if (di == 0 && dj == 0) continue;
		std::set<int> line;

		//(i,j) + k(di, dj) to get <di,dj>
		for (int k = 0; k < q; k++) {
			int nexti = (Galois::Element(gf, i) + Galois::Element(gf, di) * Galois::Element(gf, k)).value();
			int nextj = (Galois::Element(gf, j) + Galois::Element(gf, dj) * Galois::Element(gf, k)).value();
			line.insert(nexti * q + nextj);
		}

		int rep1 = *line.begin();
		int rep2 = *(++line.begin());
		
		if (encountered.find({ rep1, rep2 }) == encountered.end()) {
			for (auto& pos : line) {
				int entryVal = entry(gen);
				insert(encountered.size(), pos, Galois::Element(gf, entryVal));
			}
			
			encountered.insert({ rep1, rep2 });
		}
	}
}

int Matrix::find_primitive() {
	int q = gf->q;
	Galois::Element id(gf, 1);
	for (int i = 1; i < q; i++) {
		Galois::Element current(gf, i);
		int order = 1;
		while (current != id) {
			current = current * Galois::Element(gf, i);
			order += 1;
		}
	}
	return 1;
}


std::vector<int> Matrix::get_line() {
	int q = construction_helper->q;
	std::vector<int> line;
	Galois::Element zero(construction_helper);
	for (int i = 1; i < q; i++) {
		if (field_trace(Galois::Element(construction_helper, i)) == zero) {
			line.push_back(i);
		}
	}
	return line;
}

// Multiplicative coset
std::set<int> Matrix::get_coset(std::vector<int> line, int offset) {
	std::set<int> ans;
	for (auto& x : line) {
		ans.insert((Galois::Element(construction_helper, x) * Galois::Element(construction_helper, offset)).value());
	}
	return ans;
}

// Assume that gf is F_q^3. Tr(x) = x + x^q + x^q^2
Galois::Element Matrix::field_trace(Galois::Element x) {
	int q = cbrt(construction_helper->q);
	Galois::Element ans = x;
	
	Galois::Element inter = x;
	for (int i = 1; i < q; i++) {
		inter *= x;
	}
	ans += inter;

	Galois::Element inter2 = inter;
	for (int i = 1; i < q; i++) {
		inter2 *= inter;
	}
	ans += inter2;
	trace_range.insert(ans.value());
	return ans;
}

void Matrix::init_line_tables() {
	std::set<int> nums;
	int q = construction_helper->q;
	for (int i = 1; i < construction_helper->q; i++) {
		nums.insert(i);
	}

	int currindex = 0;

	while (nums.size()) {
		int smallest = *nums.begin();
		line_rep_to_index[smallest] = currindex;

		std::set<int> line;
		for(auto&i:trace_range){
			if (i == 0) continue;
			line.insert((Galois::Element(construction_helper, smallest) * Galois::Element(construction_helper, i)).value());
			nums.erase((Galois::Element(construction_helper, smallest) * Galois::Element(construction_helper, i)).value());
		}

		for (auto& x : line) {
			lines[x] = line;
		}

		currindex++;


	}
}

// Assuming F has cube order = q^3.
// Gives a numrows x q*2 + q + 1 matrix. each line is coset of kernel of Field trace, i.e. q^2 elts.
// rowWt = q+1lines
void Matrix::projective_geometry(int numrows) {
	std::vector<int> line = get_line();
	init_line_tables();
	int q = construction_helper->q;
	int p = cbrt(q);

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> offsetgen(1, (construction_helper->q) - 1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);

	std::set<int> encountered;
	while (encountered.size() < numrows*(p-1)) {
		// gen a random offset for the coset
		int offset = offsetgen(gen);

		if (encountered.find(offset) == encountered.end()) {
			std::set<int> coset = get_coset(line, offset);
			while (coset.size()) {
				int smallest = *coset.begin();
				int entryVal = entry(gen);
				insert(encountered.size()/(p-1), line_rep_to_index[smallest], Galois::Element(gf, entryVal));
				for (auto&i:trace_range) {
					if(i==0) continue;
					coset.erase((Galois::Element(construction_helper, smallest) * Galois::Element(construction_helper, i)).value());
				}
			}
			for (auto& x : lines[offset]) {
				encountered.insert(x);
			}
		}
	}
}

// points are F_q X F_q, lines are cosets of D = {(a, a^2)|a in F_q}, D+(h,k)
void Matrix::difference_set(int numrows) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> cosetgen(0, (gf->q) - 1);
	std::uniform_int_distribution<> entry(1, (gf->q) - 1);

	std::set<std::pair<int, int>> D = get_D();
	
	std::set<std::pair<int, int>> encountered; 
	while (encountered.size() < numrows) {
		int h = cosetgen(gen);
		int k = cosetgen(gen);
		
		if (encountered.find({ h,k }) == encountered.end()) {
			for (auto& [a, b] : D) {
				int entryVal = entry(gen);
				int i = (Galois::Element(gf, a) + Galois::Element(gf, h)).value();
				int j = (Galois::Element(gf, b) + Galois::Element(gf, k)).value();
				insert(encountered.size(), i * gf->q + j, Galois::Element(gf, entryVal));
			}
			encountered.insert({ h,k });
		}
	}
}

std::set<std::pair<int, int>> Matrix::get_D() {
	std::set<std::pair<int, int>> D;
	for (int i = 0; i < gf->q; i++) {
		D.insert({ i, (Galois::Element(gf,i) * Galois::Element(gf,i)).value() });
	}
	return D;
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

std::ostream& Matrix::printEntries(std::ostream& output){
  for (int i = 0; i < m; i++) {
    for (auto& [j, entry] : equationWiseEntries[i]) {
			output << "(" << i << ", " << j << ") = " << entry << '\n';
		}
	}
	return output;
}

std::ostream& operator<<(std::ostream& output, Matrix&matrix) {
	for (int i = 0; i < matrix.m; i++) {
		int rowindex = 0;
		if (matrix.equationWiseEntries[i].size() == 0) {
			for (int j = 0; j < matrix.n; j++) {
				output << "0 ";
			}
		}
		else {
			std::sort(matrix.equationWiseEntries[i].begin(), matrix.equationWiseEntries[i].end(), [](auto a, auto b) {
				return a.first <= b.first;
			});
			auto [jj, entry] = matrix.equationWiseEntries[i][rowindex];
			for (int j = 0; j < matrix.n; j++) {
				if (j == jj) {
					output << entry.value() << ' ';
					rowindex++;
					if (rowindex == matrix.equationWiseEntries[i].size()) {
						rowindex--;
					}
					jj = matrix.equationWiseEntries[i][rowindex].first;
					entry = matrix.equationWiseEntries[i][rowindex].second;
				}
				else {
					output << "0 ";
				}
			}
		}
		output << '\n';
	}
	return output;
}

