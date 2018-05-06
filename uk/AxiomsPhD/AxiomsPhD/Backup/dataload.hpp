#ifndef DATALOAD_HPP
#define DATALOAD_HPP

#include "axioms.hpp"
#include <string>
#include <vector>

namespace dataload
{
	dt::matrix loadDataTable(const std::string& fName);
	std::vector<double> loadTimes(const std::string& fName);
	std::vector<std::string> loadNames(const std::string& fName);
	std::vector<axioms::group> extractGroups(const dt::matrix& P, 
		const dt::matrix& X, const std::vector<std::string>& cNames, 
		const std::vector<std::string>& tNames, const std::string& fName);
	std::vector<dt::matrix> splitStrataIntoClasses(const std::vector<dt::matrix>& Xs, 
		const std::vector<std::vector<int>>& splits);
	std::vector<axioms::group> splitStrataIntoClasses(const std::vector<dt::matrix>& Xs, 
		const dt::matrix& P, const std::vector<std::string> tNames, const std::vector<std::string> cNames,
		const std::vector<std::vector<int>>& splits);
	std::vector<std::vector<int>> loadClassSplitData(const std::string& fName); // Loads data for splits variable in splitStrataIntoClasses function.
	axioms::lorenzData loadLorenzCurveData(const std::string& fName);
}

#endif