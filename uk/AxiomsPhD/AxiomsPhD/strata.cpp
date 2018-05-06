#include "dt.hpp"
#include "axioms.hpp"

dt::strata::strata(const matrix& inP, const std::vector<matrix>& inXs, 
			const std::vector<std::string>& inTimes, const std::vector<std::string>& inSNames,
			const std::vector<std::string>& inGNames, const std::string& inName) : name(inName) {
	P = inP;
	Xs = inXs;
	times = inTimes;
	sNames = inSNames;
	gNames = inGNames;
}

dt::strata dt::strata::combineStrata(std::vector<intervals> splits, const std::vector<std::string>& sNames, const std::string& name) const {
	dt::strata result;
	dt::matrix X;
	std::vector<double> row(P[0].size());
	int t, T, S, s, N, i, k;
	
	T = this->times.size();
	S = splits[0].size(); // Number of strata in result
	N = this->P[0].size();

	// Check splits for correctness
	// Number of splits
	if (splits.size() != T) {
		throw dt::error("combineStrata: number of split intervals is not equal to that of periods");
	}
	// Consistency of intervals, their limits and number of them
	for (t = 0; t < T; t++) {
		if (!splits[t].checkCosistency()) {
			throw dt::error("combineStrata: inconsistent intervals in " + splits[t].toString());
		}
		if ((splits[t][0].lb <= 0) || splits[t][splits[t].size() - 1].rb > this->Xs.size()) {
			throw dt::error("combineStrata: either split intervals " + splits[t].toString() + " cover zero or negative numbers "
				+ "or or they cover values greater than the number of available strata");
		}
		if (splits[t].size() != S) {
			throw dt::error("combineStrata: different number of untervals in splits data");
		}
	}

	result.gNames = this->gNames;
	result.name = name;
	result.P = this->P;
	result.sNames = sNames;
	result.times = this->times;

	// Fill result.Xs with empty matrices.
	for (s = 0; s < S; s++) {
		result.Xs.push_back(X);
	}
	for (t = 0; t < T; t++) {
		for (s = 0; s < S; s++) {
			if (splits[t][s] != dt::empty_interval) {
				for (i = 0; i < N; i++) {
					row[i] = 0.0;
					for (k = splits[t][s].lb - 1; k < splits[t][s].rb; k++) {
						row[i] += this->Xs[k](t, i);
					}
				}
				result.Xs[s].push_back(row);
			}
		}
	}

	return result;
}

dt::strata dt::strata::combineStrata(std::vector<intervals> splits, const std::string& name) const {
	std::vector<std::string> sNames;
	for (int s = 0; s < splits[0].size(); s++) {
		sNames.push_back("");
	}
	return this->combineStrata(splits, sNames, name);
}

std::vector<dt::tradestat> dt::strata::getClasses(const dt::splitData& splits,
									  const std::vector<std::string>& cNames) const {
	std::vector<dt::tradestat> result(splits[0].inters.size());
	std::vector<double> row(P[0].size());
	int t, T, S, s, N, i, k, p, iPeriod;
	
	T = this->times.size();
	S = splits[0].inters.size(); // Number of strata in result
	N = this->P[0].size();
	
	// Check splits for correctness
	// Consistency of intervals, their limits and number of them
	for (t = 0; t < splits.size(); t++) {
		if (splits[t].inters.size() != S) {
			throw dt::error("getClasses: different number of untervals in splits data");
		}
		for (s = 0; s < S; s++) {
			if (splits[t].inters[s] == dt::empty_interval) {
				continue;
			}
			else if ((splits[t].inters[s].lb <= 0) || splits[t].inters[s].rb > this->Xs.size()) {
				throw dt::error("getClasses: either split interval " + splits[t].inters[s].toString() + " covers zero or negative numbers "
					+ "or or it covers values greater than the number of available strata");
			}
		}
	}

	// Split itself
	for (t = 0; t < splits.size(); t++) {
		for (p = 0; p < this->times.size(); p++) {
			if (this->times[p] == splits[t].period) {
				break;
			}
		}
		if (p == this->times.size()) {
			throw dt::error("getClasses: year " + splits[t].period + " is not presented in strata times data");
		}
		for (s = 0; s < S; s++) {
			if (splits[t].inters[s] != dt::empty_interval) {
				result[s].P.push_back(this->P[p]);
				result[s].times.push_back(this->times[p]);
				for (i = 0; i < N; i++) {
					row[i] = 0.0;
					for (k = splits[t].inters[s].lb - 1; k < splits[t].inters[s].rb; k++) {
						row[i] += this->Xs[k](t, i);
					}
				}
				result[s].X.push_back(row);
			}
		}
	}

	// Copy commodity names, name the resulting trade statistics and compute Paasche matrices
	for (s = 0; s < S; s++) {
		result[s].name = cNames[s];
		result[s].gNames = this->gNames;
		result[s].logPM = axioms::logPaascheMat(axioms::dotProds(result[s].P, result[s].X));
	}

	return result;

}

std::vector<dt::tradestat> dt::strata::getClasses(const dt::splitData& splits) const {
	std::vector<std::string> cNames(splits[0].inters.size());
	return this->getClasses(splits, cNames);
}