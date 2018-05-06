#ifndef AXIOMS_HPP
#define AXIOMS_HPP

#include "dt.hpp"
#include <vector>
#include <algorithm>

namespace axioms 
{
	struct phkdIndices
	{
		std::vector<double> Q;
		std::vector<double> F;
		double omega;
	};

	struct lorenzData {
		std::vector<double> x;
		std::vector<double> y;
	};

	struct group {
		dt::matrix P; // Prices
		dt::matrix X; // Consumption
		std::string gName; // Group name
		std::vector<std::string> cNames; // Commodities names
		std::vector<std::string> times; // Periods names
	};

	const double zero = 0.000001;

	phkdIndices HARP(const dt::matrix& P, const dt::matrix& X, double logOmega);
	bool HARPTest(const dt::matrix& P, const dt::matrix& X, double logOmega);
	bool HARPTest(const dt::matrix& logPaascheMat, double logOmega);
	double HARPLogIrInd(const dt::matrix& P, const dt::matrix& X, 
		double logOlb, double logOub, double locAcc = 0.000001);
	dt::matrix dotProds(const dt::matrix& P, const dt::matrix& X);
	dt::matrix dotProds(const dt::goods& gds);
	dt::matrix logPaascheMat(const dt::matrix& dp, double logOmega = 0);
	dt::matrix transpose(const dt::matrix& mat);
	dt::matrix selectCols(const dt::matrix& mat, const std::string& cols); // Select rows from dt::matrix
	std::vector<int> intv2vec(const std::string& cols); // Convert string representation of commodity group to vector of goods indices
	double jiniIndex(const lorenzData& lcd);

	/*
	class dt::matrix
	{
		private std::vector<std::vector<double>*> *data;
		private int nr; //Number of rows
		private int nc; //Number of cols

		public dt::matrix(int r, int c)
		{
			if ((r <= 0) || (c <= 0))
			{
				throw "Bad size for matrix; File:" + __FILE__ + "; Line: " + __LINE__;
			}
			data = new std::vector<std::vector<double>*>(0);
			std::vector<double> vec;
			for (unsigned int i = 0; i < r, i++)
			{
				vec = new std::vector<double>(c);
				data->push_back(vec);
			}
			nr = r;
			nc = c;
		}

		public double& operator[](int r, int c)
	}
	*/
}

#endif //AXIOMS_HPP