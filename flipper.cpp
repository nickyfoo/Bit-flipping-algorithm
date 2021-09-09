#include "flipper.h"
#include <map>
#include <set>
#include <algorithm>
#include <unordered_set>


void Flipper::setupAL() {
	for (int i = 0; i < numVars; i++) {
		if (matrix->getDegree(i) == 2) {
			std::vector<std::pair<int, Galois::Element>> adjacentEquations = *(matrix->getAdjacentEquations(i));

			for (auto& [u, uEntry] : adjacentEquations) {
				for (auto& [v, vEntry] : adjacentEquations) {
					if (u == v) continue;
					int deg1VarNode = matrix->getIndexOfAdjSingleVar(v);
					al[u].push_back({ v,i, deg1VarNode });
				}
			}
		}
	}
}


void Flipper::FloydWarshall() {
	for (int i = 0; i < numEqns; i++) {
		for (auto& [j, intermediate, deg1Node] : al[i]) {
			APSPDistance[i][j] = 1;
			APSPNext[i][j] = { j, intermediate, deg1Node };
		}
	}

	for (int i = 0; i < numEqns; i++) {
		APSPDistance[i][i] = 0;
		APSPNext[i][i] = { i, -1, matrix->getIndexOfAdjSingleVar(i) };
	}
	for (int k = 0; k < numEqns; k++) {
		for (int i = 0; i < numEqns; i++) {
			for (int j = 0; j < numEqns; j++) {
				if (APSPDistance[i][j] > APSPDistance[i][k] + APSPDistance[k][j]) {
					APSPDistance[i][j] = APSPDistance[i][k] + APSPDistance[k][j];
					APSPNext[i][j] = APSPNext[i][k];
				}
			}
		}
	}
}

void Flipper::Dijkstra(int startNode, std::vector<int>& d, std::vector<std::tuple<int, int, int, int>>& prev) {
	std::set<std::pair<int, int>> pq;
	d.assign(numEqns, 1e9);
	prev.assign(numEqns, { -1,-1,-1,-1 });
	d[startNode] = 0;
	pq.insert({ d[startNode], startNode });
	while (pq.size()) {
		auto [du, u] = *pq.begin();
		pq.erase(pq.begin());
		if (du > d[u]) continue;
		for (auto& [v, var, deg1node] : al[u]) {
			if (du + 1 < d[v]) {
				d[v] = du + 1;
				prev[v] = { u,v,var,deg1node };
				pq.insert({ d[v],v });
			}
		}
	}
}

 std::vector<std::tuple<int, int, int>> Flipper::Dijkstra_getPath(int u, int v, std::vector<std::tuple<int, int, int, int>>& prev) {
	std::vector<std::tuple<int, int, int>> ans;
	int curr = v;
	while (curr != u) {
		auto&[uu, vv, var, deg1node] = prev[curr];
		curr = uu;
		ans.push_back({ vv, var, deg1node });
	}
	ans.push_back({ u,-1,-1 });
	std::reverse(ans.begin(), ans.end());
	return ans;
}

std::vector<std::tuple<int,int,int>> Flipper::FloydWarshall_getPath(int u, int v) {
	std::vector<std::tuple<int, int, int>> ans;
	if (APSPNext[u][v] == std::make_tuple(-1, -1, -1)) return ans;
	ans.push_back(APSPNext[u][u]);
	while (u != v) {
		auto& [nextNode, intermediate, degree1Node] = APSPNext[u][v];
		ans.push_back(APSPNext[u][v]);
		u = nextNode;
	}
	return ans;
}


void Flipper::assignX(std::vector<int> xs) {
	for (int i = 0; i < numVars; i++) {
		x[i] = Galois::Element(gf, xs[i]);
	}
}

void Flipper::assignB(std::vector<int> bs) {
	for (int i = 0; i < numEqns; i++) {
		b[i] = Galois::Element(gf, bs[i]);
	}
}

void Flipper::getXiPrimeTiBi(int i) {
	std::map<int, int> eltToNumSatisfied;
	for (auto& [j, entry] : *(matrix->getAdjacentEquations(i))) {
		Galois::Element ans = matrix->solveEqn(j, i, x, b[j]);
		eltToNumSatisfied[ans.value()]++;
	}

	std::map<int, std::vector<int>> numSatisfiedToElt;
	for (auto& [eltVal, numSatisfied] : eltToNumSatisfied) {
		numSatisfiedToElt[numSatisfied].push_back(eltVal);
	}

	auto& [numSatisfied, eltVals] = *numSatisfiedToElt.rbegin();
	t[i] = numSatisfied;
	beta[i] = (double)numSatisfied / matrix->getDegree(i);
	if (find(eltVals.begin(), eltVals.end(), x[i].value()) != eltVals.end()) {
		xiPrime[i] = Galois::Element(gf, x[i].value());
		beta[i] = -beta[i];
	}
	else {
		int pos = rand() % eltVals.size();
		xiPrime[i] = Galois::Element(gf, eltVals[pos]);
	}
}


int Flipper::getMaxBetaIndex() {
	auto it = std::max_element(beta.begin(), beta.end());
	maxBeta = *it;
	return it - beta.begin();
}

void Flipper::update(int i) {
	getXiPrimeTiBi(i);
	//update for all affected equations
	for (auto& [adjEqn, entry] : *(matrix->getAdjacentEquations(i))) {
		for (auto& [j, entry] : *(matrix->getAdjacentVariables(adjEqn))) {
			if (i != j) {
				getXiPrimeTiBi(j);
			}
		}
	}
}

bool Flipper::FloydWarshall_findFlippingPath(int startNode) {
	int endNode = -1, deg1EqnNode = -1, deg1VarNode = -1; //deg1EqnNode is adj to deg1VarNode
	bool foundPath = false;
	std::unordered_set<int> affectedNodes; // for update;
	for (int i = 0; i < numEqns; i++) {
		//keep track of deg1 nodes we have seen
		if (deg1EqnNode == -1 && APSPDistance[startNode][i] != 1e9) {
			int index = matrix->getIndexOfAdjSingleVar(i);
			if (index != -1) {
				deg1EqnNode = i;
				deg1VarNode = index;
			}
		}
		if (i == startNode) continue;
		// Reachable and not satisfied
		if (APSPDistance[startNode][i] != 1e9) {
			if (!matrix->isEqnSatisfied(i, x, b[i])) {
				endNode = i;
				break;
			}
		}
	}
	//there are 2 unsatisfied nodes
	if (endNode != -1) {
		std::vector<std::tuple<int, int, int>> path = FloydWarshall_getPath(startNode, endNode);
		//start flipping the intermediates
		for (int i = 1; i < path.size(); i++) {
			auto& [u, prevIntermediate, prevDeg1Node] = path[i - 1];
			auto& [v, intermediate, deg1Node] = path[i];
			x[intermediate] = matrix->solveEqn(u, intermediate, x, b[u]);
			affectedNodes.insert(intermediate);
		}
		foundPath = true;
	}
	else if(deg1EqnNode!=-1) {
		std::vector<std::tuple<int, int, int>> path = FloydWarshall_getPath(startNode, deg1EqnNode);
		//start flipping the intermediates
		for (int i = 1; i < path.size(); i++) {
			auto& [u, prevIntermediate, prevDeg1Node] = path[i - 1];
			auto& [v, intermediate, deg1Node] = path[i];
			x[intermediate] = matrix->solveEqn(u, intermediate, x, b[u]);
			affectedNodes.insert(intermediate);
		}
		auto& [v, intermediate, deg1Node] = *path.rbegin();
		x[deg1VarNode] = matrix->solveEqn(v, deg1VarNode, x, b[v]);
		affectedNodes.insert(deg1VarNode);
		foundPath = true;
	}

	for (auto& u : affectedNodes) {
		update(u);
	}
	return foundPath;
}

void Flipper::FloydWarshall_flip(int i) {
	//step 4
	if (maxBeta > 0) {
		x[i] = xiPrime[i];
		update(i);
	}
	//step 5
	else {
		for (auto& [adjEqn, entry] : *(matrix->getAdjacentEquations(i))) {
			if (!matrix->isEqnSatisfied(adjEqn, x, b[adjEqn])) {
				if (FloydWarshall_findFlippingPath(adjEqn)) {
					return;
				}
			}
		}
	}
}


void Flipper::Dijkstra_flip(int i) {
	//step 4
	if (maxBeta > 0) {
		x[i] = xiPrime[i];
		update(i);
	}
	//step 5
	else {
		int startNode;
		int endNode = -1, deg1EqnNode = -1, deg1VarNode = -1; //deg1EqnNode is adj to deg1VarNode
		std::unordered_set<int> affectedNodes; // for update;
		for (auto& [adjEqn, entry] : *(matrix->getAdjacentEquations(i))) {
			if (!matrix->isEqnSatisfied(adjEqn, x, b[adjEqn])) {
				if (Dijkstra_findFlippingPath(adjEqn)) {
					return;
				}
			}
		}
	}

}


bool Flipper::Dijkstra_findFlippingPath(int startNode) {
	std::vector<int> d;
	std::vector<std::tuple<int, int, int, int>> prev;
	Dijkstra(startNode, d, prev);
	int endNode = -1, deg1EqnNode = -1, deg1VarNode = -1; //deg1EqnNode is adj to deg1VarNode
	bool foundPath = false;
	std::unordered_set<int> affectedNodes; // for update;
	for (int i = 0; i < numEqns; i++) {
		//keep track of deg1 nodes we have seen
		if (deg1EqnNode == -1 && d[i] != 1e9) {
			int index = matrix->getIndexOfAdjSingleVar(i);
			if (index != -1) {
				deg1EqnNode = i;
				deg1VarNode = index;
			}
		}
		if (i == startNode) continue;
		// Reachable and not satisfied
		if (d[i] != 1e9) {
			if (!matrix->isEqnSatisfied(i, x, b[i])) {
				endNode = i;
				break;
			}
		}
	}
	//there are 2 unsatisfied nodes
	if (endNode != -1) {
		std::vector<std::tuple<int, int, int>> path = Dijkstra_getPath(startNode, endNode, prev);
		//start flipping the intermediates
		for (int i = 1; i < path.size(); i++) {
			auto& [u, prevIntermediate, prevDeg1Node] = path[i - 1];
			auto& [v, intermediate, deg1Node] = path[i];
			x[intermediate] = matrix->solveEqn(u, intermediate, x, b[u]);
			affectedNodes.insert(intermediate);
		}
		foundPath = true;
	}
	else if (deg1EqnNode != -1) {
		std::vector<std::tuple<int, int, int>> path = Dijkstra_getPath(startNode, deg1EqnNode, prev);
		//start flipping the intermediates
		for (int i = 1; i < path.size(); i++) {
			auto& [u, prevIntermediate, prevDeg1Node] = path[i - 1];
			auto& [v, intermediate, deg1Node] = path[i];
			x[intermediate] = matrix->solveEqn(u, intermediate, x, b[u]);
			affectedNodes.insert(intermediate);
		}
		auto& [v, intermediate, deg1Node] = *path.rbegin();
		x[deg1VarNode] = matrix->solveEqn(v, deg1VarNode, x, b[v]);
		affectedNodes.insert(deg1VarNode);
		foundPath = true;
	}

	for (auto& u : affectedNodes) {
		update(u);
	}
	return foundPath;
}

//Flips one at a time, checking for new max beta value every time.
void Flipper::solve_extended_bit_flipping_consecutively() {
	for (int i = 0; i < numVars; i++) {
		getXiPrimeTiBi(i);
	}
	//floydWarshall();
	int index = getMaxBetaIndex();
	std::cout << *this << '\n';
	while (beta[index] != -1) {
		std::cout << "current beta is: " << beta[index] << ", flipping: " << index << '\n';
		Dijkstra_flip(index);
		std::cout << *this << '\n';
		index = getMaxBetaIndex();
	}
	std::cout << *this << '\n';
}


void Flipper::solve_extended_bit_flipping_with_Omega() {
	for (int i = 0; i < numVars; i++) {
		getXiPrimeTiBi(i);
	}
	FloydWarshall();
	int index = getMaxBetaIndex();
	std::set<int> Omega;
	for (int i = 0; i < numVars; i++) {
		if (beta[i] == maxBeta) Omega.insert(i);
	}
	std::cout << "Omega contains: ";
	for (auto& x : Omega) std::cout << x << ' ';
	std::cout << '\n';
	std::cout << *this << '\n';
	while (beta[index] != -1) {
		while (Omega.size() && maxBeta != -1) {
			std::cout << "current beta is: " << maxBeta<< ", current Omega size is: " << Omega.size() << ", flipping: " << index << '\n';

			std::cout << "Omega contains: ";
			for (auto& x : Omega) std::cout << x << ' ';
			std::cout << '\n';
			FloydWarshall_flip(*Omega.begin());
			std::set<int> updatedOmega;
			for (auto&i : Omega) {
				if (beta[i] == maxBeta) {
					std::cout << i << ' ' << beta[i] << ' ' << maxBeta << '\n';
					updatedOmega.insert(i);
				}
			}
			Omega = updatedOmega;
			std::cout << *this << '\n';
		}
		index = getMaxBetaIndex();
		for (int i = 0; i < numVars; i++) {
			if (beta[i] == maxBeta) Omega.insert(i);
		}
	}
}

// returns failure flag
bool Flipper::solve_single_threshold_decoding(int threshold) {
	std::cout << "Original x is: ";
	for (auto& xi : x) std::cout << xi << ' ';
	std::cout << '\n';
	bool flipped = true;
	while (flipped) {
		flipped = false;
		for (int i = 0; i < numVars; i++) {
			//Calculate messages for xi
			int numZeroMessages = 0;
			std::map<int, int> eltToNumSatisfied;
			for (auto& [j, entry] : *(matrix->getAdjacentEquations(i))) {
				Galois::Element ans = matrix->solveEqn(j, i, x, b[j]);
				if (ans == x[i]) numZeroMessages++;
				else eltToNumSatisfied[ans.value()]++;
			}

			//std::cout << "eltToNumSatisfied: ";
			//for (auto& [eltVal, numSatisfied]:eltToNumSatisfied) {
			//	std::cout << "(" << eltVal << "," << numSatisfied << ") ";
			//}
			//std::cout << '\t';
			//std::cout << "numZeroMessages: " << numZeroMessages << "\n";

			std::map<int,int> numSatisfiedToElt;
			for (auto& [eltVal, numSatisfied] : eltToNumSatisfied) {
				numSatisfiedToElt[numSatisfied]=eltVal;
			}
			if (numSatisfiedToElt.size() == 0) continue;
			auto& [numSatisfied, eltVal] = *numSatisfiedToElt.rbegin();
			if (numSatisfied - numZeroMessages > threshold) {
				std::cout << i << ": numSatisfied = " << numSatisfied << ", numZeroMessages = " << numZeroMessages << ", threshold = " << threshold << ", Flipping " << i << '\n';
				x[i] = Galois::Element(gf, eltVal);
				flipped = true;
				std::cout << "x is: ";
				for (auto& xi : x) std::cout << xi << ' ';
				std::cout << '\n';
			}
		}
	}
	for (int i = 0; i < numEqns; i++) {
		if (!matrix->isEqnSatisfied(i, x, b[i])) return false;
	}
	return true;
}

bool Flipper::solve_multiple_threshold_decoding(std::vector<int> thresholds) {
	int numThresholds = thresholds.size();
	sort(thresholds.begin(), thresholds.end());
	for (int i = numThresholds-1; i>=0; i--){
		int threshold = thresholds[i];
		std::cout << "Multiple threshold decoding with threshold = " << threshold << '\n';
		solve_single_threshold_decoding(threshold);
	}

	for (int i = 0; i < numEqns; i++) {
		if (!matrix->isEqnSatisfied(i, x, b[i])) return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& output, Flipper& flipper) {
	output << "numEqns: " << flipper.numEqns << " numVars: " << flipper.numVars << '\n';
	output << '\n';
	output << "Input variables x_i:\n";
	for (auto&elt : flipper.x) {
		output << elt << '\t';
	}
	output << '\n';

	output << "Solution b_i:\n";
	for (auto& elt : flipper.b) {
		output << elt << '\t';
	}
	output << '\n';

	output << "Ideal values x_i':\n";
	for (auto& elt : flipper.xiPrime) {
		output << elt << '\t';
	}
	output << '\n';

	output << "Num satisfied eqns by ideal val t_i':\n";
	for (auto& elt : flipper.t) {
		output << elt << '\t';
	}
	output << '\n';

	output << "Reliability ratio beta_i':\n";
	for (auto& elt : flipper.beta) {
		output << elt << '\t';
	}
	output << '\n';

	return output;
}