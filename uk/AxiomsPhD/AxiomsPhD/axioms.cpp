#include "axioms.hpp"
#include "strfuns.hpp"
#include <vector>
#include <string>
#include <cmath>

// Price and quantity matrices: row<->time, column<->commodity

dt::indices axioms::HARP(const dt::matrix& P, const dt::matrix& X, double logOmega)
{
	bool integr = true;
	dt::matrix dp = axioms::dotProds(P, X);
	dt::matrix logPMO = axioms::logPaascheMat(dp, logOmega);
	dt::matrix newLogPMO;
	int T = logPMO.size();
	int t, tau, i;
	std::vector<double> row(T);
	std::vector<double> Q(T);
	std::vector<double> F(T);
	double sum1, sum2;
	dt::indices result;

	// Fill newLogPMO with zeros
	for (t = 0; t < T; t++)
	{
		newLogPMO.push_back(row);
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
					newLogPMO[t][tau] = std::min(logPMO[t][tau],
						logPMO[t][i] + logPMO[i][tau]);
					if (newLogPMO[t][t] < -axioms::zero)
					{
						integr = false;
					}
				}
				else
				{
					logPMO[t][tau] = std::min(newLogPMO[t][tau],
						newLogPMO[t][i] + newLogPMO[i][tau]);
					if (logPMO[t][t] < -axioms::zero)
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
				sum1 = logPMO[t][0];
			}
			else
			{
				sum1 = newLogPMO[t][0];
			}
			for (tau = 0; tau < T; tau++)
			{
				if (T % 2)
				{
					sum2 = logPMO[t][tau];
				}
				else 
				{
					sum2 = newLogPMO[t][tau];
				}
				sum1 = std::min(sum1, sum2);
			}
			sum2 = std::exp(sum1);
			Q[t] = 1 / sum2;
			F[t] = sum2 * dp[t][t];
		}
		result.x = F;
		result.p = Q;
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

bool axioms::HARPTest(const dt::matrix& P, const dt::matrix& X, double logOmega)
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
/*
dt::matrix axioms::dotProds(const dt::goods& gds) {
	// result[t,tau] = <P^t, X^{\tau}>
	int t, tau, i;
	int T = gds[0]->p.size();
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
				sum += gds[i]->p[t] * gds[i]->x[tau];
			}
			row[tau] = sum;
		}
		result.push_back(row);
	}
	return result;
}
*/
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

double axioms::giniIndex(const dt::lorenzData& lcd) {
	double sum;
	int i;

	sum = 0.0;
	for (i = 2; i < lcd.x.size(); i++) {
		sum += (lcd.x[i] - lcd.x[i - 1]) * lcd.y[i];
	}
	return (0.5 - sum) * 2;
}

bool axioms::hav_sys_solvability(const dt::matrix& A, double r)
{
	// Checks if the system a_{ij} + x_j \leqslant r + x_i has positive solution in xs.

	int n = A.size();
	dt::matrix A_with_r, new_A_with_r;
	std::vector<double> row(n);
	int i, j, k;

	// Initialize A_with_r
	for (i = 0; i < n; i++)
	{
		row.clear();
		for (j = 0; j < n; j++)
		{
			row.push_back(A[i][j] - r);
		}
		row[i] = 0.0;
		A_with_r.push_back(row);
	}

	// Initialize new_A_with_r
	row.clear();
	for (i = 0; i < n; i++)
	{
		row.push_back(0.0);
	}
	for (i = 0; i < n; i++)
	{
		new_A_with_r.push_back(row);
	}

	// The test itself
	for (k = 1; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (k % 2)
				{
					new_A_with_r[i][j] = std::max(A_with_r[i][j], 
						A_with_r[i][k] + A_with_r[k][j]);
					if (new_A_with_r[i][i] > axioms::zero)
					{
						return false;
					}
				}
				else
				{
					A_with_r[i][j] = std::max(new_A_with_r[i][j], 
						new_A_with_r[i][k] + new_A_with_r[k][j]);
					if (A_with_r[i][i] > axioms::zero)
					{
						return false;
					}
				}
			}
		}
	}
	return true;
}

std::vector<double> axioms::hav_sys_solve(const dt::matrix& A, double r)
{
	// This attempts to solve the system a_{ij} + x_j \leqslant r + x_i in terms of xs.

	int n = A.size();
	dt::matrix A_with_r, new_A_with_r;
	std::vector<double> row(n);
	int i, j, k;
	bool solvable = true;
	std::vector<double> x(n);
	double max1;

	x.clear();

	// Initialize A_with_r
	for (i = 0; i < n; i++)
	{
		row.clear();
		for (j = 0; j < n; j++)
		{
			row.push_back(A[i][j] - r);
		}
		row[i] = 0.0;
		A_with_r.push_back(row);
	}

	// Initialize new_A_with_r
	row.clear();
	for (i = 0; i < n; i++)
	{
		row.push_back(0.0);
	}
	for (i = 0; i < n; i++)
	{
		new_A_with_r.push_back(row);
	}

	// Find A^{*n}
	for (k = 1; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (k % 2)
				{
					new_A_with_r[i][j] = std::max(A_with_r[i][j], 
						A_with_r[i][k] + A_with_r[k][j]);
					if (new_A_with_r[i][i] > axioms::zero)
					{
						// The system has no solution :(
						return x;
					}
				}
				else
				{
					A_with_r[i][j] = std::max(new_A_with_r[i][j], 
						new_A_with_r[i][k] + new_A_with_r[k][j]);
					if (A_with_r[i][i] > axioms::zero)
					{
						// The system has no solution :(
						return x;
					}
				}
			}
		}
	}
	
	// Find a solution.
	for (i = 0; i < n; i++)
	{
		if (n % 2)
		{
			max1 = A_with_r[i][0];
		}
		else
		{
			max1 = new_A_with_r[i][0];
		}
		for (j = 1; j < n; j++)
		{
			if (n % 2)
			{
				max1 = std::max(max1, A_with_r[i][j]);
			}
			else
			{
				max1 = std::max(max1, new_A_with_r[i][j]);
			}
		}
		x.push_back(max1);
	}
	return x;
}

double axioms::hav_sys_min_r(const dt::matrix& A, double lb, double rb, double eps)
{
	bool hav_left = axioms::hav_sys_solvability(A, lb);
	bool hav_right = axioms::hav_sys_solvability(A, rb);
	double mid_point = (rb + lb) / 2;
	bool hav_mid = axioms::hav_sys_solvability(A, mid_point);
	double diff;

	if (hav_left)
	{
		return lb;
	}
	else if (!hav_right)
	{
		return rb;
	}
	else
	{
		while(true)
		{
			if (hav_mid)
			{
				rb = mid_point;
			}
			else
			{
				lb = mid_point;
			}
			if (rb - lb < eps)
			{
				return rb;
			}
			mid_point = (lb + rb) * 0.5;
			hav_mid = axioms::hav_sys_solvability(A, mid_point);
		}
	}
}