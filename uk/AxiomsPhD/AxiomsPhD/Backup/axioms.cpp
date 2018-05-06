#include "axioms.hpp"
#include "strfuns.hpp"
#include <vector>
#include <string>
#include <cmath>

// Price and quantity matrices: row<->time, column<->commodity

axioms::phkdIndices axioms::HARP(const dt::matrix& P, 
								 const dt::matrix& X, double logOmega)
{
	dt::matrix dp = axioms::dotProds(P, X);
	dt::matrix logPM = axioms::logPaascheMat(dp, logOmega);
	bool integr = true;
	dt::matrix newLogPM;
	int T = P.size();
	int t, tau, i;
	std::vector<double> row(T);
	std::vector<double> Q(T);
	std::vector<double> F(T);
	double sum1, sum2;
	axioms::phkdIndices result;

	// Fill newLogPM with zeros
	for (t = 0; t < T; t++)
	{
		newLogPM.push_back(row);
	}

	// Idempotent degree T
	for (i = 1; i < T; i++)
	{
		for (t = 0; t < T; t++)
		{
			for (tau = 0; tau < T; tau++)
			{
				if (i % 2)
				{
					newLogPM[t][tau] = std::min(logPM[t][tau],
						logPM[t][i] + logPM[i][tau]);
					if (newLogPM[t][t] < -axioms::zero)
					{
						integr = false;
					}
				}
				else
				{
					logPM[t][tau] = std::min(newLogPM[t][tau],
						newLogPM[t][i] + newLogPM[i][tau]);
					if (logPM[t][t] < -axioms::zero)
					{
						integr = false;
					}
				}
			}
		}
	}

	// PHKD indices
	if (integr)
	{
		for (t = 0; t < T; t++)
		{
			if (T % 2)
			{
				sum1 = logPM[t][0];
			}
			else
			{
				sum1 = newLogPM[t][0];
			}
			for (tau = 0; tau < T; tau++)
			{
				if (T % 2)
				{
					sum2 = logPM[t][tau];
				}
				else 
				{
					sum2 = newLogPM[t][tau];
				}
				sum1 = std::min(sum1, sum2);
			}
			sum2 = std::exp(sum1);
			Q[t] = 1 / sum2;
			F[t] = sum2 * dp[t][t];
		}
		result.F = F;
		result.Q = Q;
		result.omega = std::exp(logOmega);
	}
	else
	{
		result.omega = -1;
	}

	return result;
}

bool axioms::HARPTest(const dt::matrix& logPaascheMat,
					  double logOmega) {
	int i, t, tau;
	dt::matrix logPM; // log Paasche matrix plus log(omega)
	dt::matrix newLogPM; // for idempotent degree
	int T = logPaascheMat.size();
	std::vector<double> row(T);

	// Add log(omega) to logPaasheMat
	for (t = 0; t < T; t++) {
		for (tau = 0; tau < T; tau++) {
			row[tau] = logPaascheMat[t][tau] + logOmega;
		}
		logPM.push_back(row);
	}
	// For allowing omega < 1
	for (t = 0; t < T; t++) {
		logPM[t][t] = 0;
	}

	// Fill newLogPM with zeros
	row.clear();
	for (tau = 0; tau < T; tau++) {
		row.push_back(0);
	}
	for (t = 0; t < T; t++)
	{
		newLogPM.push_back(row);
	}

	for (i = 1; i < T; i++)
	{
		for (t = 0; t < T; t++)
		{
			for (tau = 0; tau < T; tau++)
			{
				if (i % 2)
				{
					newLogPM[t][tau] = std::min(logPM[t][tau],
						logPM[t][i] + logPM[i][tau]);
					if (newLogPM[t][t] < -axioms::zero)
					{
						return false;
					}
				}
				else
				{
					logPM[t][tau] = std::min(newLogPM[t][tau],
						newLogPM[t][i] + newLogPM[i][tau]);
					if (logPM[t][t] < -axioms::zero)
					{
						return false;
					}
				}
			}
		}
	}

	return true;
}

bool axioms::HARPTest(const dt::matrix& P,
					  const dt::matrix& X, double logOmega)
{
	dt::matrix dp = axioms::dotProds(P, X);
	dt::matrix logPM = axioms::logPaascheMat(dp);
	return axioms::HARPTest(logPM, logOmega);
}

double axioms::HARPLogIrInd(const dt::matrix& P, const dt::matrix& X,
						 double logOlb, double logOub, 
						 double logAcc) {
	
	// result < logOlb -> HARP satisfied with log(omega) = logOlb
	// result > logOub -> HARP is not satisfied with log(omega) = logOub
	// result \in [logOlb, logOub] -> irrationality index found

	dt::matrix dp = axioms::dotProds(P, X);
	dt::matrix logPM = axioms::logPaascheMat(dp);
	double lb = logOlb;
	double ub = logOub;
	double mb;
	bool lbHARP, ubHARP, mbHARP;
	
	lbHARP = axioms::HARPTest(logPM, lb);
	ubHARP = axioms::HARPTest(logPM, ub);
	if (lbHARP && ubHARP) {
		return logOlb - 1;
	}
	else if (!(lbHARP || ubHARP)) {
		return logOub + 1;
	}

	while (std::abs(ub - lb) > logAcc) {
		mb = (ub + lb) / 2;
		mbHARP = axioms::HARPTest(logPM, mb);
		if (mbHARP) {
			ub = mb;
		}
		else {
			lb = mb;
		}
	}
	return ub;
}

dt::matrix axioms::dotProds(const dt::matrix& P, 
								const dt::matrix& X)
{
	// result[t,tau] = <P^t, X^{\tau}>
	int t, tau, i;
	int T = P.size();
	int N = P[0].size();
	dt::matrix result;
	std::vector<double> row(T);
	double sum;
	for (t = 0; t < T; t++)
	{
		for (tau = 0; tau < T; tau++)
		{
			sum = 0.0;
			for (i = 0; i < N; i++)
			{
				sum += P[t][i] * X[tau][i];
			}
			row[tau] = sum;
		}
		result.push_back(row);
	}
	return result;
}

dt::matrix axioms::dotProds(const dt::goods& gds) {
	// result[t,tau] = <P^t, X^{\tau}>
	int t, tau, i;
	int T = gds[0].p.size();
	int N = gds.size();
	dt::matrix result;
	std::vector<double> row(T);
	double sum;
	for (t = 0; t < T; t++)
	{
		for (tau = 0; tau < T; tau++)
		{
			sum = 0.0;
			for (i = 0; i < N; i++)
			{
				sum += gds[i].p[t] * gds[i].x[tau];
			}
			row[tau] = sum;
		}
		result.push_back(row);
	}
	return result;
}

dt::matrix axioms::logPaascheMat(
	const dt::matrix& dp, double logOmega)
{
	// result[t,tau] = \log(\frac{<P^{tau},X^{tau}>}{<P^t, X^{tau}>});
	int T = dp.size();
	dt::matrix result;
	std::vector<double> row(T);
	int t, tau;

	for (t = 0; t < T; t++)
	{
		for (tau = 0; tau < T; tau++)
		{
			if (dp[tau][tau] <= axioms::zero)
			{
				row[tau] = -10 + logOmega;
			}
			else
			{
				row[tau] = std::log(dp[tau][t] / dp[t][t]) + logOmega;
			}
		}
		result.push_back(row);
	}

	for (t = 0; t < T; t++)
	{
		result[t][t] = 0;
	}

	return result;
}

dt::matrix axioms::transpose(const dt::matrix& mat) {
	// result = mat^{T}
	int m = mat.size();
	int n = mat[0].size();
	int i, j;
	dt::matrix result;
	std::vector<double> row(m);

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			row[i] = mat[i][j];
		}
		result.push_back(row);
	}

	return result;
}

dt::matrix axioms::selectCols(const dt::matrix& mat, const std::string& cols) {
	// cols format: interval_1, interval_2, ..., interval_N
	// interval_i = lb-rb or lr
	// cols specify matrix indices (starting from 1) for matrix columns
	
	dt::matrix result;
	std::vector<std::string> strs1, strs2;
	std::string str;
	int i, j;
	int lb, rb;
	dt::matrix matT = axioms::transpose(mat);

	strs1 = strfuns::split(cols, ',');
	for (i = 0; i < strs1.size(); i++) {
		strs2 = strfuns::split(strs1[i], '-');
		if ((strs2.size() > 2) || (strs2.size() <= 0)) {
			throw "Incorrect columns interval " + strs1[i];
		}
		else if (strs2.size() == 2) {
			lb = std::stoi(strs2[0]);
			rb = std::stoi(strs2[1]);
			if (lb > rb) {
				throw "Left bound greater than the right one in interval " + strs1[i];
			}
			for (j = lb - 1; j < rb; j++) {
				result.push_back(matT[j]);
			}
		}
		else {
			lb = std::stoi(strs2[0]);
			result.push_back(matT[lb - 1]);
		}
	}

	return axioms::transpose(result);
}

std::vector<int> axioms::intv2vec(const std::string& cols) {
	std::vector<int> result;
	std::vector<std::string> strs1, strs2;
	std::string str;
	int i, j;
	int lb, rb;

	strs1 = strfuns::split(cols, ',');
	for (i = 0; i < strs1.size(); i++) {
		strs2 = strfuns::split(strs1[i], '-');
		if ((strs2.size() > 2) || (strs2.size() <= 0)) {
			throw "Incorrect columns interval " + strs1[i];
		}
		else if (strs2.size() == 2) {
			lb = std::stoi(strs2[0]);
			rb = std::stoi(strs2[1]);
			if (lb > rb) {
				throw "Left bound greater than the right one in interval " + strs1[i];
			}
			for (j = lb - 1; j < rb; j++) {
				result.push_back(j);
			}
		}
		else {
			lb = std::stoi(strs2[0]);
			result.push_back(lb - 1);
		}
	}

	return result;
}

double axioms::jiniIndex(const axioms::lorenzData& lcd) {
	double sum;
	int i;

	sum = 0.0;
	for (i = 2; i < lcd.x.size(); i++) {
		sum += (lcd.x[i] - lcd.x[i - 1]) * lcd.y[i];
	}
	return (0.5 - sum) * 2;
}