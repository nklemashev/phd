#include "dataload.hpp"
#include "axioms.hpp"
#include <fstream>
#include <vector>
#include <string>
#include "strfuns.hpp"
#include <cstdlib>
#include <algorithm>

dt::matrix dataload::loadDataTable(const std::string& fName)
{
	std::ifstream fin(fName);
	dt::matrix result;
	std::vector<double> row;
	std::string str;
	std::vector<std::string> strs;
	unsigned int i;

	if (!fin)
	{
		throw "Unable to open file " + fName;
	}

	std::getline(fin, str);
	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, ',');
		row.clear();
		for (i = 1; i < strs.size(); i++) {
			row.push_back(std::stod(strs[i]));
		}
		result.push_back(row);
	}
	fin.close();
	return result;
}

std::vector<std::string> dataload::loadNames(const std::string& fName) {
	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> result;
	
	if (!fin)
	{
		throw "Unable to open file " + fName;
	}

	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		result.push_back(str);
	}
	fin.close();
	return result;
}

std::vector<double> dataload::loadTimes(const std::string& fName) {
	std::ifstream fin(fName);
	std::string str;
	std::vector<double> result;
	
	if (!fin)
	{
		throw "Unable to open file " + fName;
	}

	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		result.push_back(std::stod(str));
	}
	fin.close();
	return result;
}

std::vector<axioms::group> dataload::extractGroups(const dt::matrix& P, 
		const dt::matrix& X, const std::vector<std::string>& cNames,
		const std::vector<std::string>& tNames, const std::string& fName) {
	// Extracts groups from TS given by P and X.
	// Lists of commodities in groups with groups' names are in file
	// with full address given in fName.
	// The file contains lines in the form of intervals|group_name.
	// Intervals are in form num[-num[,num[-num[...]]]].
	// num is the number of good (starting from 1, not 0)
	// Examples of intervals:
	//  1-4,5-9
	//  4,5-8,9,43

	std::vector<axioms::group> result;
	axioms::group gr;
	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> strs, strs2;
	int i, j, lb, rb, t;
	std::vector<int> inds;
	std::vector<double> rowP, rowX;

	if (!fin) {
		throw "Unable to open file " + fName;
	}
	for (t = 0; t < tNames.size(); t++) {
		gr.times.push_back(tNames[t]);
	}
	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, '|');
		if (strs.size() > 2) {
			throw "extractGroups::Unable to parse group line " + str + "; too many |'s";
		}
		gr.gName = strs[1];
		gr.P.clear();
		gr.X.clear();
		gr.cNames.clear();
		inds.clear();
		strs = strfuns::split(strs[0], ',');
		for (i = 0; i < strs.size(); i++) {
			strs2 = strfuns::split(strs[i], '-');
			if (strs2.size() > 2) {
				throw "extractGroups::Unable to parse interval " + strs[i] + "; too many -'s";
			}
			else if (strs2.size() == 2) {
				lb = std::stoi(strs2[0]);
				rb = std::stoi(strs2[1]);
				if (rb < lb) {
					throw "extractGroups::Left boundary is greater than the right one in interval " + strs[i];
				}
				for (j = lb; j <= rb; j++) {
					if (std::find(inds.begin(), inds.end(), j - 1) != inds.end()) {
						throw "extractGroups::Commodity " + std::to_string(j) + " is specified twice in " + str;
					}
					else {
						inds.push_back(j - 1);
					}
				}
			}
			else if (strs2.size() == 1) {
				lb = std::stoi(strs2[0]);
				if (std::find(inds.begin(), inds.end(), lb - 1) != inds.end()) {
					throw "extractGroups::Commodity " + strs2[0] + " is specified twice in " + str;
				}
				inds.push_back(lb - 1);
			}
			else {
				throw "extractGroups::Unable to parse empty interval in " + str;
			}
		}

		std::sort(inds.begin(), inds.end());
		rowP.clear();
		rowX.clear();
		for (i = 0; i < inds.size(); i++) {
			gr.cNames.push_back(cNames[inds[i]]);
			rowP.push_back(0.0);
			rowX.push_back(0.0);
		}
		for (t = 0; t < P.size(); t++) {
			for (i = 0; i < inds.size(); i++) {
				rowP[i] = P[t][inds[i]];
				rowX[i] = X[t][inds[i]];
			}
			gr.P.push_back(rowP);
			gr.X.push_back(rowX);
		}
		result.push_back(gr);
	}
	return result;
}

std::vector<dt::matrix> dataload::splitStrataIntoClasses(const std::vector<dt::matrix>& Xs, 
	const std::vector<std::vector<int>>& splits) {

	// split[t] -- split data for year t
	// Xs[m] consumption matrix for statum m.
	//
	// Example with three classes:
	// splits[t] = (i, j) such that
	// 0..(i - 1) -- poor class
	// i..(j - 1) -- middle class
	// j..(Xs.size() - 1) -- rich class
	// i \in [1, Xs.size() - 2]
	// j \in [i + 1, Xs.size() - 1]

	std::vector<dt::matrix> result;
	int i, j, k, t;
	std::vector<int> genSplit;
	std::vector<double> row;
	dt::matrix X;
	int T, N, S, C;
	double sum;

	S = Xs.size(); // Num of strata
	T = Xs[0].size();
	N = Xs[0][0].size();
	C = splits[0].size() + 1; // Num of classes

	// Start with consistency checks
	if (splits.size() != T) {
		throw std::string("Incorrect split data: number of splits differs from number of periods.");
	}
	for (t = 0; t < T; t++) {
		// Is there any split data?
		if (splits[t].size() == 0) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + " : empty split.";
		}
		// Is the number of classes the same for all periods?
		if (splits[t].size() != C - 1) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + " : number of splits for this period " +
				"differs from that for the first period.";
		}
		// What about strict monotonicity?
		for (i = 1; i < splits[t].size(); i++) {
			if (splits[t][i] < splits[t][i - 1]) {
				throw "Incorrect split data for period " + std::to_string(t + 1) + " : split indices are not strictly monotone.";
			}
		}
		// What about bounding restrictions?
		if (splits[t][0] <= 0) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + 
				" : the first index should be greater than zero.";
		}
		if (splits[t][splits[t].size() - 1] >= S) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + 
				" : the last index should be less than " + std::to_string(S) + '.';
		}
	}

	// Initialize row, genSplit, X and result
	for (i = 0; i < N; i++) {
		row.push_back(0.0);
	}
	for (t = 0; t < C; t++) {
		genSplit.push_back(0);
	}
	genSplit.push_back(S);
	for (t = 0; t < T; t++) {
		X.push_back(row);
	}
	for (j = 0; j < C; j++) {
		result.push_back(X);
	}
	row.clear();
	X.clear();

	// Splitting
	for (t = 0; t < T; t++) {
		// Generalized split data
		for (j = 0; j < C - 1; j++) {
			genSplit[j + 1] = splits[t][j];
		}		

		// Splitting itself
		for (j = 1; j < genSplit.size(); j++) {
			for (i = 0; i < N; i++) {
				sum = 0.0;
				for (k = genSplit[j - 1]; k < genSplit[j]; k++) {
					sum += Xs[k][t][i];
				}
				result[j - 1][t][i] = sum;
			}
		}
	}
	return result;
}

std::vector<axioms::group> dataload::splitStrataIntoClasses(const std::vector<dt::matrix>& Xs, 
		const dt::matrix& P, const std::vector<std::string> tNames, const std::vector<std::string> cNames,
		const std::vector<std::vector<int>>& splits)
{
	// split[t] -- split data for year t
	// Xs[m] consumption matrix for statum m.
	//
	// Example with three classes:
	// splits[t] = (i, j) such that
	// 0..(i - 1) -- poor class
	// i..(j - 1) -- middle class
	// j..(Xs.size() - 1) -- rich class
	// i \in [1, Xs.size() - 2]
	// j \in [i + 1, Xs.size() - 1]
	//
	// Ver 2 (why do we need price matrix and to return vector of groups)
	// Some classes may dissappear/appear with time.
	// Now functions accepts split data with first index being zero
	// (meaning no first class), last index being equal to the Xs.size()
	// (meaning no last class), and with several indices having same value
	// (meaning no intermediate classes. That is why we need prices (there are
	// less price data for classes which dissappear in some periods) and
	// to return vector of groups. Notice that the number of split numbers
	// still should be the same for each split line.

	std::vector<axioms::group> result;
	int i, j, k, t;
	std::vector<int> genSplit;
	std::vector<double> row;
	axioms::group G;
	int T, N, S, C;
	double sum;

	S = Xs.size(); // Num of strata
	T = Xs[0].size();
	N = Xs[0][0].size();
	C = splits[0].size() + 1; // Num of classes

	// Start with consistency checks
	if (splits.size() != T) {
		throw std::string("Incorrect split data: number of splits differs from number of periods.");
	}
	for (t = 0; t < T; t++) {
		// Is there any split data?
		if (splits[t].size() == 0) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + " : empty split.";
		}
		// Is the number of classes the same for all periods?
		if (splits[t].size() != C - 1) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + " : number of splits for this period " +
				"differs from that for the first period.";
		}
		// What about monotonicity?
		for (i = 1; i < splits[t].size(); i++) {
			if (splits[t][i] < splits[t][i - 1]) {
				throw "Incorrect split data for period " + std::to_string(t + 1) + " : split indices are not monotone.";
			}
		}
		// What about bounding restrictions?
		if (splits[t][0] < 0) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + 
				" : the first index should be nonnegative.";
		}
		if (splits[t][splits[t].size() - 1] >= S) {
			throw "Incorrect split data for period " + std::to_string(t + 1) + 
				" : the last index should be less than or equal to " + std::to_string(S) + '.';
		}
	}

	// Initialize row, genSplit, G, and result
	for (i = 0; i < N; i++) {
		row.push_back(0.0);
	}
	G.cNames = cNames;
	for (t = 0; t < C; t++) {
		genSplit.push_back(0);
		G.gName = "Class " + strfuns::int2str(t + 1);
		result.push_back(G);
	}
	genSplit.push_back(S);

	// Splitting
	for (t = 0; t < T; t++) {
		// Generalized split data
		for (j = 0; j < C - 1; j++) {
			genSplit[j + 1] = splits[t][j];
		}		

		// Splitting itself
		for (j = 1; j < genSplit.size(); j++) {
			if (genSplit[j] > genSplit[j - 1]) {
				for (i = 0; i < N; i++) {
					sum = 0.0;
					for (k = genSplit[j - 1]; k < genSplit[j]; k++) {
						sum += Xs[k][t][i];
					}
					row[i] = sum;
				}
			}
			result[j - 1].X.push_back(row);
			result[j - 1].P.push_back(P[t]);
			result[j - 1].times.push_back(tNames[t]);
		}
	}
	return result;
}

std::vector<std::vector<int>> dataload::loadClassSplitData(const std::string& fName) {
	// Split indices are given by string of indices separated by comma.
	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> strs;
	std::vector<int> splitData;
	std::vector<std::vector<int>> result;
	int i;

	if (!fin) {
		throw "Unable to open file " + fName;
	}

	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, ',');
		splitData.clear();
		for (i = 0; i < strs.size(); i++) {
			splitData.push_back(std::stoi(strs[i]));
		}
		result.push_back(splitData);
	}
	return result;
}

std::vector<dt::intervals> dataload::loadClassSplitDataNew(const std::string& fName) {
	// File with splits data should contain lines in form <interval>{,<interval>{,<interval>...}}
	// <interval> should be either in the form of lb-rb, where lb and rb defines 
	// percentages of population bounding the specific class, or empty string
	// which means empty interval.
	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> strs, strs2;
	dt::intervals splitData;
	std::vector<dt::intervals> result;
	int i;

	if (!fin) {
		throw "Unable to open file " + fName;
	}

	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, ',');
		splitData.clear();
		for (i = 0; i < strs.size(); i++) {
			if (strs[i].size() == 0) {
				splitData.push_back(dt::empty_interval);
				continue;
			}
			strs2 = strfuns::split(strs[i], '-');
			if (strs2.size() != 2) {
				throw dt::error("loadClassSplitData: Incorrect interval format in " + strs[i]);
			}
			splitData.push_back(dt::interval(std::stoi(strs2[0]), std::stoi(strs2[1])));
		}
		result.push_back(splitData);
	}
	return result;
}

dt::splitData dataload::loadClassSplitDataNew2(const std::string& fName) {
	// File with splits data should contain lines in form <period>:<interval>{,<interval>{,<interval>...}}
	// <interval> should be either in the form of lb-rb, where lb and rb defines 
	// percentages of population bounding the specific class, or empty string
	// which means empty interval.
	// <period> specifies period for which the split data comes after double dot.
	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> strs, strs2;
	dt::splitLine sLine;
	dt::splitData result;
	int i;

	if (!fin) {
		throw "Unable to open file " + fName;
	}

	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, ':');
		if (strs.size() != 2) {
			throw dt::error("loadClassSplitData: Incorrect split line format in " + str);
		}
		sLine.period = strs[0];
		sLine.inters.clear();
		strs = strfuns::split(strs[1], ',');
		for (i = 0; i < strs.size(); i++) {
			if (strs[i].size() == 0) {
				sLine.inters.push_back(dt::empty_interval);
				continue;
			}
			strs2 = strfuns::split(strs[i], '-');
			if (strs2.size() != 2) {
				throw dt::error("loadClassSplitData: Incorrect interval format in " + strs[i]);
			}
			sLine.inters.push_back(dt::interval(std::stoi(strs2[0]), std::stoi(strs2[1])));
		}
		result.push_back(sLine);
	}
	return result;
}

dt::lorenzData dataload::loadLorenzCurveData(const std::string& fName) {

	std::ifstream fin(fName);
	std::string str;
	std::vector<std::string> strs;
	dt::lorenzData result;

	if (!fin) {
		throw "loadLorenzCurveData: Unable to open file " + fName;
	}
	std::getline(fin, str);
	// Add zeros
	result.x.push_back(0.0);
	result.y.push_back(0.0);
	while (true) {
		std::getline(fin, str);
		if (fin.eof()) {
			break;
		}
		strs = strfuns::split(str, ',');
		result.x.push_back(std::stod(strs[0]));
		result.y.push_back(std::stod(strs[1]));
	}
	return result;
}
