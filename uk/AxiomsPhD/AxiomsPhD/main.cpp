#include "dataload.hpp"
#include "axioms.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include "strfuns.hpp"
#include <cmath>
#include <algorithm>

int tenPartsByIncome() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::vector<std::string> names;
	std::vector<double> periods;
	dt::matrix P;
	std::vector<dt::matrix> X;
	int i;
	double omega = 0;
	double omegaMax = 0;
	dt::matrix Q, F, QT, FT;
	dt::indices inds;
	std::ofstream report;

	try {
		report.open(path + "\\Ten parts by income\\report.txt");
		names = dataload::loadNames(path + "comNames.csv");
		periods = dataload::loadTimes(path + "perNums.csv");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		std::cout << "P: " << P.size() << "x" << P[0].size() << std::endl;
		for (i = 0; i <= 10; i++) {
			X.push_back(dataload::loadDataTable(path + "Ten parts by income\\Part" + 
				strfuns::int2str(i) + ".csv"));
			std::cout << "X[" << i << "]: " << X[i].size() << "x" << X[i][0].size() << std::endl;
			omega = std::exp(axioms::HARPLogIrInd(
				P, X[i], std::log(0.5), std::log(2)));
			report << "Group " << i << ": ";
			if (omega > 2) {
				std::cout << "Omega > 2" << std::endl;
				report << "Omega > 2" << std::endl;
			}
			/*else if (omega < 1) {
				std::cout << "Omega < 1" << std::endl;
				report << "Omega < 1" << std::endl;
			}*/
			else {
				std::cout << "Omega = " << omega << std::endl;
				report << "Omega = " << omega << std::endl;
				omegaMax = std::max(omegaMax, omega);
			}
		}
		omegaMax += axioms::zero;
		report << "Selected omega for further analysis: " << omegaMax << std::endl;
		for (i = 1; i <= 10; i++) {
			inds = axioms::HARP(P, X[i], std::log(omegaMax));
			if (inds.omega < 0) {
				throw "Unable to integrate group " + strfuns::int2str(i);
			}
			else {
				Q.push_back(inds.p);
				F.push_back(inds.x);
			}
		}

		FT = axioms::transpose(F);
		std::cout << "FT: " << FT.size() << "x" << FT[0].size() << std::endl;
		QT = axioms::transpose(Q);
		std::cout << "QT: " << QT.size() << "x" << QT[0].size() << std::endl;
		omega = std::exp(axioms::HARPLogIrInd(QT, FT, std::log(1), std::log(2)));
		std::cout << "Omega: " << omega << std::endl;
		report << "Omega for TS created from PHKD indices: " << omega << std::endl;
	}
	catch (std::string err) {
		std::cout << err << std::endl;
	}

	report.close();

	return 0;
}

int severalGroups() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::vector<std::string> namesRus;
	dt::matrix P, PS;
	dt::matrix X, XS;
	int i, j;
	double omega = 0;
	std::ofstream report;
	std::ifstream fin;
	std::vector<std::string> inds, gNames;
	std::string str;
	std::vector<std::string> strs;
	std::vector<int> indsVec;

	try {
		report.open(path + "Some groups\\report.txt");
		if (!report) {
			throw "Unable to open " + path + "Some groups\\report.txt";
		}
		fin.open(path + "Some groups\\Groups.txt");
		if (!fin) {
			throw "Unable to open " + path + "Some groups\\Groups.txt";
		}

		// Load group names and intervals
		while(true) {
			std::getline(fin, str);
			if (fin.eof()) {
				break;
			}
			strs = strfuns::split(str, '|');
			inds.push_back(strs[0]);
			gNames.push_back(strs[1]);
		}

		// Load data
		X = dataload::loadDataTable(path + "\\Ten parts by income\\Part0.csv");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		namesRus = dataload::loadNames(path + "comNamesRus.csv");

		// Compute omega for groups
		for (i = 0; i < inds.size(); i++) {
			PS = axioms::selectCols(P, inds[i]);
			XS = axioms::selectCols(X, inds[i]);
			omega = std::exp(axioms::HARPLogIrInd(PS, XS,
				std::log(1), std::log(2)));
			std::cout << "Group " << i + 1 << "omega : " 
				<< omega << std::endl;
			// Report
			report << gNames[i] << std::endl;
			indsVec = axioms::intv2vec(inds[i]);
			report << namesRus[indsVec[0]];
			for (j = 1; j < indsVec.size(); j++) {
				report << "; " << namesRus[indsVec[j]];
			}
			report << std::endl;
			report << "Omega: " << omega << std::endl << std::endl;
		}
	}
	catch (std::string err) {
		std::cout << "ERROR!" << err << std::endl;
	}

	report.close();
	return 0;
}

int dataForIndex() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::vector<std::string> names;
	std::vector<double> periods;
	dt::matrix P;
	dt::matrix X;
	std::ofstream foutTor(path + "BCC2008 for Index\\BCC2008.tor");
	std::ofstream foutSet(path + "BCC2008 for Index\\BCC2008.set");
	std::ofstream foutLng(path + "BCC2008 for Index\\BCC2008.lng");
	std::ofstream foutS(path + "BCC2008 for Index\\BCC2008.s");
	int T, N, t, i;

	try {
		names = dataload::loadNames(path + "comNamesRus.csv");
		periods = dataload::loadTimes(path + "perNums.csv");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		std::cout << "P: " << P.size() << "x" << P[0].size() << std::endl;
		X = dataload::loadDataTable(path + "Ten parts by income\\Part0.csv");
		std::cout << "X: " << X.size() << "x" << X[0].size() << std::endl;
		T = X.size();
		N = X[0].size();

		// tor
		foutTor << "Time" << std::endl;
		foutTor << T << std::endl;
		foutTor << "Prod" << std::endl;
		foutTor << N << std::endl;
		for (t = 0; t < T; t++) {
			for (i = 0; i < N; i++) {
				foutTor << P[t][i] << std::endl;
				foutTor << X[t][i] << std::endl;
			}
		}
		foutTor.close();

		// lng
		foutLng << "1" << std::endl;
		foutLng << "Все товары и услуги" << std::endl;
		foutLng << "1" << std::endl;
		foutLng << N << std::endl;
		foutLng << "$" << std::endl;
		foutLng.close();

		// set
		foutSet << T << std::endl;
		for (t = 0; t < T; t++) {
			foutSet << periods[t] << std::endl;
		}
		foutSet.close();

		// s
		for (i = 0; i < N; i++) {
			foutS << i + 1 << std::endl;
			foutS << names[i] << std::endl;
		}
		foutS << '$' << std::endl;
		foutS.close();
	}
	catch (std::string err) {
		std::cout << err << std::endl;
	}

	return 0;
}

int dutch() {
	std::vector<dt::matrix> P, X;
	int i;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\Dutch for Stata\\";
	double omega;

	for (i = 0; i <= 5; i++) {
		P.push_back(dataload::loadDataTable(path + "P" + strfuns::int2str(i) + ".csv"));
		std::cout << "P[" << i << "]: " << P[i].size() << "x" << P[i][0].size() << std::endl;
		X.push_back(dataload::loadDataTable(path + "X" + strfuns::int2str(i) + ".csv"));
		std::cout << "X[" << i << "]: " << X[i].size() << "x" << X[i][0].size() << std::endl;
		omega = std::exp(axioms::HARPLogIrInd(P[i], X[i], std::log(0.5), std::log(2), 0.0000001));
		std::cout << "Omega: " << omega << std::endl;
	}

	return 0;
}

int compareIndices1() {
	// Compare PHKD indices calculated in three ways:
	// 1) For whole TS
	// 2) For whole TS with aggregating from data on social strata
	// 3) For whole TS with aggregating from data on subgroups

	dt::indices w1, w2, w3, ind;
	std::ofstream report, reportOnIndices2, reportOnIndices3;
	dt::matrix P, X, PS, XS;
	dt::matrix Q, F;
	std::vector<std::string> namesRus;
	double omega;
	std::ifstream fin;
	std::vector<std::string> inds, gNames, strs;
	std::string str;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	int i, j, k;
	double omegaMax;
	std::vector<int> indsVec;
	std::ofstream fout;

	try {
		namesRus = dataload::loadNames(path + "comNamesRus.csv");
		report.open(path + "Compare indices 1\\report.txt");
		reportOnIndices2.open(path + "Compare indices 1\\reportOnIndices2.csv");
		reportOnIndices3.open(path + "Compare indices 1\\reportOnIndices3.csv");

		// Way 1
		report << "Way 1" << std::endl;
		X = dataload::loadDataTable(path + "Ten parts by income\\Part0.csv");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		omega = std::exp(axioms::HARPLogIrInd(P, X, 0, std::log(2)));
		if (omega < 1) {
			omega = 1;
		}
		else if (omega > 2) {
			throw "Unable to find omega for the whole group";
		}
		std::cout << "Way 1" << std::endl << "Omega = " << omega << std::endl;
		report << "Omega = " << omega << std::endl;
		w1 = axioms::HARP(P, X, std::log(omega));
		if (w1.omega < 0) {
			throw "Unable to integrate whole group with omega = " + std::to_string(omega);
		}
		report << std::endl;

		// Way 2
		report << "Way 2" << std::endl;
		std::cout << "Way 2" << std::endl;
		omegaMax = 0;
		for (i = 1; i <= 10; i++) {
			X = dataload::loadDataTable(path + "Ten parts by income\\Part" + strfuns::int2str(i) + ".csv");
			omega = std::exp(axioms::HARPLogIrInd(P, X, 0, std::log(2)));
			if (omega < 1) {
				omega = 1;
			}
			else if (omega > 2) {
				throw "Unable to find omega for the group " + std::to_string(i);
			}
			std::cout << "Group " << i << " omega = " << omega << std::endl;
			report << "Group " << i << " omega = " << omega << std::endl;
			omegaMax = std::max(omega, omegaMax);
		}
		std::cout << "Max omega = " << omegaMax << std::endl;
		report << "Max omega = " << omegaMax << std::endl;
		for (i = 1; i <= 10; i++) {
			X = dataload::loadDataTable(path + "Ten parts by income\\Part" + strfuns::int2str(i) + ".csv");
			ind = axioms::HARP(P, X, std::log(omegaMax));
			if (ind.omega < 0) {
				throw "Unable to integrate group " + std::to_string(i) + " with omega = " + std::to_string(omegaMax);
			}
			reportOnIndices2 << ind.p[0];
			for (k = 1; k < ind.p.size(); k++) {
				reportOnIndices2 << ',' << ind.p[k];
			}
			reportOnIndices2 << ",1" ;
			for (k = 1; k < ind.x.size(); k++) {
				reportOnIndices2 << ',' << ind.x[k] / ind.x[0];
			}
			reportOnIndices2 << std::endl;
			Q.push_back(ind.p);
			F.push_back(ind.x);
		}
		Q = axioms::transpose(Q);
		F = axioms::transpose(F);
		std::cout << "Q: " << Q.size() << 'x' << Q[0].size() << std::endl;
		std::cout << "F: " << F.size() << 'x' << F[0].size() << std::endl;
		omega = std::exp(axioms::HARPLogIrInd(Q, F, 0, std::log(2)));
		if (omega < 1) {
			omega = 1;
		}
		else if (omega > 2) {
			throw "Unable to find omega for the group created from trade statistics of stata";
		}
		report << "Omega for trade statistics of PHKD indices = " << omega << std::endl;
		w2 = axioms::HARP(Q, F, std::log(omega));
		if (w2.omega < 0) {
			throw "Unable to integrate whole group with omega = " + std::to_string(omega);
		}
		report << std::endl;

		// Way 3
		report << "Way 3" << std::endl;
		std::cout << "Way 3" << std::endl;
		omegaMax = 0;
		F.clear();
		Q.clear();
		fin.open(path + "Compare indices 1\\ComGroups.txt");
		if (!fin) {
			throw "Unable to open file " + path + "Compare indices 1\\ComGroups.txt";
		}
		while (true) {
			std::getline(fin, str);
			if (fin.eof()) {
				break;
			}
			strs = strfuns::split(str, '|');
			inds.push_back(strs[0]);
			gNames.push_back(strs[1]);
		}

		X = dataload::loadDataTable(path + "\\Ten parts by income\\Part0.csv");
		P = dataload::loadDataTable(path + "pallForStata.csv");

		// Compute omega for groups
		for (i = 0; i < inds.size(); i++) {
			PS = axioms::selectCols(P, inds[i]);
			XS = axioms::selectCols(X, inds[i]);
			omega = std::exp(axioms::HARPLogIrInd(PS, XS,
				std::log(1), std::log(2)));
			if (omega < 1) {
				omega = 1;
			}
			else if (omega > 2) {
				throw "Unable to find omega for the group " + std::to_string(i);
			}
			omegaMax = std::max(omega, omegaMax);
			std::cout << "Group " << i + 1 << "omega : " 
				<< omega << std::endl;
			// Report
			report << gNames[i] << std::endl;
			indsVec = axioms::intv2vec(inds[i]);
			report << namesRus[indsVec[0]];
			for (j = 1; j < indsVec.size(); j++) {
				report << "; " << namesRus[indsVec[j]];
			}
			report << std::endl;
			report << "Omega: " << omega << std::endl << std::endl;
		}

		for (i = 0; i < inds.size(); i++) {
			PS = axioms::selectCols(P, inds[i]);
			XS = axioms::selectCols(X, inds[i]);
			ind = axioms::HARP(PS, XS, std::log(omegaMax));
			if (ind.omega < 0) {
				throw "Unable to integrate group " + std::to_string(i + 1) + " with omega = " + std::to_string(omegaMax);
			}
			reportOnIndices3 << ind.p[0];
			for (k = 1; k < ind.p.size(); k++) {
				reportOnIndices3 << ',' << ind.p[k];
			}
			reportOnIndices3 << ",1";
			for (k = 1; k < ind.x.size(); k++) {
				reportOnIndices3 << ',' << ind.x[k] / ind.x[0];
			}
			reportOnIndices3 << std::endl;
			F.push_back(ind.x);
			Q.push_back(ind.p);
		}
		// Add chemist's goods
		XS = axioms::transpose(axioms::selectCols(X, "55"));
		PS = axioms::transpose(axioms::selectCols(P, "55"));
		reportOnIndices3 << '1';
		for (k = 1; k < PS[0].size(); k++) {
			reportOnIndices3 << ',' << PS[0][k] / PS[0][0];
		}
		reportOnIndices3 << ",1";
		for (k = 1; k < XS[0].size(); k++) {
			reportOnIndices3 << ',' << XS[0][k] / XS[0][0];
		}
		reportOnIndices3 << std::endl;
		F.push_back(XS[0]);
		Q.push_back(PS[0]);
		
		F = axioms::transpose(F);
		Q = axioms::transpose(Q);
		std::cout << "Q: " << Q.size() << 'x' << Q[0].size() << std::endl;
		std::cout << "F: " << F.size() << 'x' << F[0].size() << std::endl;
		omega = std::exp(axioms::HARPLogIrInd(Q, F, 0, std::log(2)));
		if (omega < 1) {
			omega = 1;
		}
		else if (omega > 2) {
			throw "Unable to find omega for the group created from the trade statistics of subgroups";
		}
		report << "Omega for trade statistics of PHKD indices = " << omega << std::endl;
		w3 = axioms::HARP(Q, F, std::log(omega));
		if (w3.omega < 0) {
			throw "Unable to integrate the group created from the trade statistics of subgroups with omega = " + std::to_string(omega);
		}

		// Output indices
		fout.open(path + "Compare indices 1\\indices.csv");
		if (!fout) {
			throw "Unable to open file " + path + "Compare indices 1\\indices.csv";
		}
		fout << "Q1,F1,Q2,F2,Q3,F3" << std::endl;
		for (i = 0; i < w1.x.size(); i++) {
			fout << w1.p[i] << ',' << w1.x[i] << ','
				<< w2.p[i] << ',' << w2.x[i] << ','
				<< w3.p[i] << ',' << w3.x[i] << ',' << std::endl;
		}
		fout.close();
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	report.close();
	return 0;

}

int point2_post() {
	std::ifstream report;
	std::ofstream strataOmegas;
	std::vector<std::string> strs, strsStrata, strsOmegas;
	std::string str;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	int i;

	try {
		report.open(path + "Point 2\\report.txt");
		if (!report) {
			throw "Unable to open file " + path + "Point 2\\report.txt";
		}
		strataOmegas.open(path + "Point 2\\strataOmegas.csv");
		if (!strataOmegas) {
			throw "Unable to open file " + path + "Point 2\\strataOmegas.csv";
		}
		strataOmegas << "Share 1,Share 2,Share 3,Omega 1,Omega 2,Omega 3, Omega All" << std::endl;
		while (true) {
			std::getline(report, str);
			if (report.eof()) {
				break;
			}
			strs = strfuns::split(str, '|');
			strsStrata = strfuns::split(strs[0], ',');
			strsOmegas = strfuns::split(strs[1], ',');
			strataOmegas << strsStrata[0];
			for (i = 1; i < strsStrata.size(); i++) {
				strataOmegas << ',' << strsStrata[i];
			}
			for (i = 0; i < strsOmegas.size(); i++) {
				strataOmegas << ',' << strsOmegas[i];
			}
			strataOmegas << ',' << strs[8];
			strataOmegas << std::endl;
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}
	return 0;
}

int point2() {
	int i, j, ii, tt, jj;
	dt::matrix P, Q, F;
	std::vector<dt::matrix> X3;
	std::vector<double> lOmega;
	std::vector<double> row, col;
	double lOmegaMax;
	dt::indices ind;
	std::ofstream report;
	std::vector<dt::matrix> Xs;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	double sum, llb = std::log(2);
	bool failed = false;

	try {
		report.open(path + "Point 2\\report.txt");
		if (!report) {
			throw "Unable to open file " + path + "Point 2\\report.txt";
		}

		P = dataload::loadDataTable(path + "pallForStata.csv");
		for (i = 1; i <= 10; i++) {
			Xs.push_back(dataload::loadDataTable(path + "Ten parts by income\\Part" + std::to_string(i) + ".csv"));
			std::cout << "X" << i << ": " << Xs[i - 1].size() << 'x' << Xs[i - 1][0].size() << std::endl;
		}
		for (i = 0; i < P[0].size(); i++) {
			row.push_back(0.0);
		}
		for (i = 0; i < 3; i++) {
			X3.push_back(Xs[0]);
			lOmega.push_back(0.0);
		}

		// 0..(i - 1) -- poor class
		// i..(j - 1) -- middle class
		// j..(Xs.size() - 1) -- rich class
		// i \in [1, Xs.size() - 2]
		// j \in [i + 1, Xs.size() - 1]
		for (i = 1; i < Xs.size() - 1; i++) {
			for (j = i + 1; j < Xs.size(); j++) {
				std::cout << "(i, j) = (" << i << ", " << j << ')' << std::endl;
				for (ii = 0; ii < 3; ii++) {
					X3[ii].clear();
				}
				// poor class
				for (tt = 0; tt < P.size(); tt++) {
					for (ii = 0; ii < P[0].size(); ii++) {
						sum = 0.0;
						for (jj = 0; jj < i; jj++) {
							sum += Xs[jj][tt][ii];
						}
						row[ii] = sum;
					}
					X3[0].push_back(row);
				}
				// middle class
				for (tt = 0; tt < P.size(); tt++) {
					for (ii = 0; ii < P[0].size(); ii++) {
						sum = 0.0;
						for (jj = i; jj < j; jj++) {
							sum += Xs[jj][tt][ii];
						}
						row[ii] = sum;
					}
					X3[1].push_back(row);
				}
				// rich class
				for (tt = 0; tt < P.size(); tt++) {
					for (ii = 0; ii < P[0].size(); ii++) {
						sum = 0.0;
						for (jj = j; jj < Xs.size(); jj++) {
							sum += Xs[jj][tt][ii];
						}
						row[ii] = sum;
					}
					X3[2].push_back(row);
				}

				// Computations
				for (ii = 0; ii < 3; ii++) {
					lOmega[ii] = axioms::HARPLogIrInd(P, X3[ii], 0, llb);
					if (lOmega[ii] < 0) {
						lOmega[ii] = 0;
					}
				}
				lOmegaMax = lOmega[0];
				for (ii = 1; ii < 3; ii++) {
					lOmegaMax = std::max(lOmegaMax, lOmega[ii]);
				}
				report << std::to_string(i * 10) << ',' << std::to_string((j - i) * 10) << ',' << std::to_string((10 - j) * 10) << '|';
				report << std::exp(lOmega[0]) << ',' << std::exp(lOmega[1]) << ',' << std::exp(lOmega[2]) << '|';
				Q.clear();
				F.clear();
				for (ii = 0; ii < 3; ii++) {
					if (lOmega[ii] > llb) {
						failed = true;
						report << '0';
						for (tt = 1; tt < P.size(); tt++) {
							report << ",0";
						}
						report << "|0";
						for (tt = 1; tt < P.size(); tt++) {
							report << ",0";
						}
						report << '|';
					}
					else {
						ind = axioms::HARP(P, X3[ii], lOmegaMax);
						Q.push_back(ind.p);
						F.push_back(ind.x);
						report << ind.p[0];
						for (tt = 1; tt < P.size(); tt++) {
							report << ',' << ind.p[tt];
						}
						report << '|' << ind.x[0];
						for (tt = 1; tt < P.size(); tt++) {
							report << ',' << ind.x[tt];
						}
						report << '|';
					}
				}
				if (!failed) {
					F = axioms::transpose(F);
					Q = axioms::transpose(Q);
					lOmegaMax = axioms::HARPLogIrInd(Q, F, 0, llb);
					if (lOmegaMax < 0) {
						lOmegaMax = 0;
					}
					report << std::exp(lOmegaMax) << '|';
					if (lOmegaMax > llb) {
						report << '0';
						for (tt = 1; tt < P.size(); tt++) {
							report << ",0";
						}
						report << "|0";
						for (tt = 1; tt < P.size(); tt++) {
							report << ",0";
						}
						report << "|";
					}
					else {
						ind = axioms::HARP(Q, F, lOmegaMax);
						report << ind.p[0];
						for (tt = 1; tt < P.size(); tt++) {
							report << ',' << ind.p[tt];
						}
						report << '|' << ind.x[0];
						for (tt = 1; tt < P.size(); tt++) {
							report << ',' << ind.x[tt];
						}
						report << '|';
					}
				}
				report << std::endl;
			}
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	return 0;
}

int point3() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	std::ofstream report;
	double logOmega, logOmegaMin;
	dt::matrix P, X, Pc, Xc;
	std::vector<std::string> times;
	int i, t, tt, tMin;

	try {
		report.open(path + "Point 3\\report.txt");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		times = dataload::loadNames(path + "perNums.csv");
		for (i = 1; i <= 10; i++) {
			std::cout << "Group " << i << ": omega = ";
			report << "Group " << i << ": omega = ";
			X = dataload::loadDataTable(path + "Ten parts by income\\Part" + std::to_string(i) + ".csv");
			logOmega = axioms::HARPLogIrInd(P, X, llb, lrb);
			std::cout << std::exp(logOmega) << std::endl;
			report << std::exp(logOmega) << std::endl;
			logOmegaMin = 100500;
			for (t = 0; t < times.size(); t++) {
				std::cout << "Exclude year " << times[t] << " -> omega = ";
				report << "Exclude year " << times[t] << " -> omega = ";
				Pc.clear();
				Xc.clear();
				for (tt = 0; tt < times.size(); tt++) {
					if (tt == t) {
						continue;
					}
					Pc.push_back(P[tt]);
					Xc.push_back(X[tt]);
				}
				logOmega = axioms::HARPLogIrInd(Pc, Xc, llb, lrb);
				if (logOmega < llb) {
					logOmega = llb;
				}
				else if (logOmega > lrb) {
					logOmega = std::log(100500);
				}
				std::cout << std::exp(logOmega) << std::endl;
				report << std::exp(logOmega) << std::endl;
				if (logOmega < logOmegaMin) {
					tMin = t;						
					logOmegaMin = logOmega;
				}
			}
			report << "Maximum decrease in omega is attained from the exclusion of year " << times[tMin] << std::endl << std::endl;
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	return 0;
}

int proportions() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	std::ofstream reportOnIndicesQ, reportOnIndicesF, report;
	std::vector<axioms::group> groups;
	dt::matrix P, X;
	int i, g, t;
	std::vector<std::string> namesRus, times;
	dt::indices indAll, ind;
	double lomegaMax, lomega;
	dt::matrix tmpMat;

	try {
		report.open(path + "Proportions for 100\\report.txt");
		P = dataload::loadDataTable(path + "pallForStata.csv");
		namesRus = dataload::loadNames(path + "comNamesRus.csv");
		times = dataload::loadNames(path + "perNums.csv");
		for (i = 1; i <= 100; i++) {
			std::cout << "Stratum " << i << std::endl;
			report << "Stratum " << i << std::endl;
			X = dataload::loadDataTable(path + "Hundred parts by income\\Part" + std::to_string(i) + ".csv");
			lomegaMax = axioms::HARPLogIrInd(P, X, llb, lrb);
			if (lomegaMax < llb) {
				lomegaMax = llb;
			}
			std::cout << "Omega for all goods = " << std::exp(lomegaMax) << std::endl;
			report << "Omega for all goods = " << std::exp(lomegaMax) << std::endl;
			groups = dataload::extractGroups(P, X, namesRus, times, path + "Proportions for 100\\ComGroups.txt");
			for (g = 0; g < groups.size(); g++) {
				if (groups[g].P[0].size() == 1) {
					continue;
				}
				lomega = axioms::HARPLogIrInd(groups[g].P, groups[g].X, llb, lrb);
				if (lomega < llb) {
					lomega = llb;
				}
				std::cout << "Omega for group " + std::to_string(g + 1) << " = " << std::exp(lomega) << std::endl;
				report << "Omega for group " + std::to_string(g + 1) << " = " << std::exp(lomega) << std::endl;
				lomegaMax = std::max(lomega, lomegaMax);
			}
			std::cout << std::endl;
			report << std::endl;

			indAll = axioms::HARP(P, X, lomegaMax);
			reportOnIndicesQ.open(path + "Proportions for 100\\QS" + std::to_string(i) + ".txt");
			reportOnIndicesF.open(path + "Proportions for 100\\FS" + std::to_string(i) + ".txt");
			reportOnIndicesQ << indAll.p[0];
			reportOnIndicesF << indAll.x[0];
			for (t = 1; t < P.size(); t++) {
				reportOnIndicesQ << ' ' << indAll.p[t];
				reportOnIndicesF << ' ' << indAll.x[t];
			}
			reportOnIndicesQ << std::endl;
			reportOnIndicesF << std::endl;
			for (g = 0; g < groups.size(); g++) {
				if (groups[g].P[0].size() == 1) {
					tmpMat = axioms::transpose(groups[g].P);
					ind.p = tmpMat[0];
					tmpMat = axioms::transpose(groups[g].X);
					ind.x = tmpMat[0];
					ind.omega = 0;
				}
				else {
					ind = axioms::HARP(groups[g].P, groups[g].X, lomegaMax);
				}
				reportOnIndicesQ << ind.p[0];
				reportOnIndicesF << ind.x[0];
				for (t = 1; t < P.size(); t++) {
					reportOnIndicesQ << ' ' << ind.p[t];
					reportOnIndicesF << ' ' << ind.x[t];
				}
				reportOnIndicesQ << std::endl;
				reportOnIndicesF << std::endl;
			}
			reportOnIndicesQ.close();
			reportOnIndicesF.close();
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}
	return 0;
}

int proportions3() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::string pathToStrataConsumption = path + "Hundred parts by income\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	double lOmegaMax, lOmega;
	std::ofstream reportOnIndicesQ, reportOnIndicesF, report;
	std::vector<axioms::group> groups;
	dt::matrix P, X, tmpMat;
	std::vector<dt::matrix> Xs, Xc;
	std::vector<std::vector<int>> splits;
	int i, g, t;
	int S = 100;
	std::vector<std::string> namesRus, times;
	dt::indices indAll, ind;

	try {
		report.open(path + "Proportions 3\\report.txt");
		P = dataload::loadDataTable(path + "pallForStata.csv");		
		namesRus = dataload::loadNames(path + "comNamesRus.csv");
		times = dataload::loadNames(path + "perNums.csv");
		// Load strata consumption data
		std::cout << "Loading consumption data for " << S << " strata..." << std::endl;
		for (i = 1; i <= S; i++) {
			Xs.push_back(dataload::loadDataTable(pathToStrataConsumption + "Part" + std::to_string(i) + ".csv"));
		}
		std::cout << "Loading consumption data for " << S << " strata: Complete" << std::endl;
		// Join strata into three classes
		splits = dataload::loadClassSplitData(path + "Proportions 3\\Splits.txt");
		Xc = dataload::splitStrataIntoClasses(Xs, splits);

		// Get proportions
		for (i = 0; i < Xc.size(); i++) {
			std::cout << "Class " << i + 1 << std::endl;
			report << "Class " << i + 1 << std::endl;
			groups = dataload::extractGroups(P, Xc[i], namesRus, times, path + "Proportions 3\\ComGroups.txt");
			
			// Omega for all goods
			lOmegaMax = axioms::HARPLogIrInd(P, Xc[i], llb, lrb);
			if (lOmegaMax < llb) {
				lOmegaMax = llb;
			}
			std::cout << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;
			report << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;

			// Omegas for groups
			for (g = 0; g < groups.size(); g++) {
				if (groups[g].cNames.size() > 1) {
					lOmega = axioms::HARPLogIrInd(groups[g].P, groups[g].X, llb, lrb);
					if (lOmega < llb) {
						lOmega = llb;
					}
					std::cout << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
					report << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
					lOmegaMax = std::max(lOmegaMax, lOmega);
				}
			}
			std::cout << std::endl;
			report << std::endl;

			// PHKD indices
			indAll = axioms::HARP(P, Xc[i], lOmegaMax);
			reportOnIndicesQ.open(path + "Proportions 3\\QS" + std::to_string(i) + ".txt");
			reportOnIndicesF.open(path + "Proportions 3\\FS" + std::to_string(i) + ".txt");
			reportOnIndicesQ << indAll.p[0];
			reportOnIndicesF << indAll.x[0];
			for (t = 1; t < P.size(); t++) {
				reportOnIndicesQ << ' ' << indAll.p[t];
				reportOnIndicesF << ' ' << indAll.x[t];
			}
			reportOnIndicesQ << std::endl;
			reportOnIndicesF << std::endl;
			for (g = 0; g < groups.size(); g++) {
				if (groups[g].P[0].size() == 1) {
					tmpMat = axioms::transpose(groups[g].P);
					ind.p = tmpMat[0];
					tmpMat = axioms::transpose(groups[g].X);
					ind.x = tmpMat[0];
					ind.omega = 0;
				}
				else {
					ind = axioms::HARP(groups[g].P, groups[g].X, lOmegaMax);
				}
				reportOnIndicesQ << ind.p[0];
				reportOnIndicesF << ind.x[0];
				for (t = 1; t < P.size(); t++) {
					reportOnIndicesQ << ' ' << ind.p[t];
					reportOnIndicesF << ' ' << ind.x[t];
				}
				reportOnIndicesQ << std::endl;
				reportOnIndicesF << std::endl;
			}
			reportOnIndicesQ.close();
			reportOnIndicesF.close();
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	return 0;
}

int proportions3Inf() {
	char c;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::string pathToStrataConsumption = path + "Hundred parts by income\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	double lOmegaMax, lOmega;
	std::ofstream reportOnIndicesQ, reportOnIndicesF, reportOnJini, report, toMatlab;
	std::vector<axioms::group> groups;
	dt::matrix P, X, tmpMat;
	std::vector<dt::matrix> Xs, Xc;
	std::vector<dt::intervals> splits;
	std::vector<std::vector<dt::simpleDInterval>> splitsJini;
	std::vector<dt::simpleDInterval> vecInter;
	dt::simpleDInterval inter(0.0, 0.0);
	int i, g, t;
	int S = 100;
	std::vector<std::string> namesRus, times;
	dt::indices indAll, ind;
	std::vector<dt::lorenzData> lcds, lcdsClasses;
	dt::matrix* jiniIndices;
	dt::strata strats, stratsClasses;

	try {
		std::cout << "Python: This is proportions3Inf" << std::endl;
		strats.P = dataload::loadDataTable(path + "pallForStata.csv");		
		strats.gNames = dataload::loadNames(path + "comNamesRus.csv");
		strats.times = dataload::loadNames(path + "perNums.csv");
		// Load strata consumption data
		std::cout << "Loading consumption data for " << S << " strata..." << std::endl;
		for (i = 1; i <= S; i++) {
			strats.Xs.push_back(dataload::loadDataTable(pathToStrataConsumption + "Part" + std::to_string(i) + ".csv"));
			strats.sNames.push_back("");
		}
		std::cout << "Loading consumption data for " << S << " strata: Complete" << std::endl;
		std::cout << "Loading data for Lorenz curve..." << std::endl;
		for (i = 75; i <= 99; i++) {
			lcds.push_back(dataload::loadLorenzCurveData(path + "Jini\\Year" + strfuns::int2str(i) + ".csv"));
		}
		std::cout << "Loading data for Lorenz curve: Complete" << std::endl;
		std::cout << "Python:OK" << std::endl;
		while (true) {
			std::cout << "proportions3Inf?" << std::endl;
			std::cin >> c;
			if (c == 'y') {
				report.open(path + "Proportions 3\\report.txt");
				// Join strata into three classes
				splits = dataload::loadClassSplitDataNew(path + "Proportions 3\\SplitsNew.txt");
				stratsClasses = strats.combineStrata(splits);
				// Save number of classes for Matlab
				toMatlab.open(path + "Proportions 3\\toMatlab.txt");
				toMatlab << splits[0].size() << std::endl;
				toMatlab.close();
				
				// Get Jini indices
				/*
				std::cout << "Computing Jini indices..." << std::endl;
				splitsJini.clear();
				vecInter.clear();
				inter.lb = 0.0;
				inter.rb = 0.0;
				vecInter.push_back(inter);
				for (t = 0; t < splits.size(); t++) {
					vecInter[0].rb = (double)splits[t][0].rb / 100;
					splitsJini.push_back(vecInter);
				}
				for (t = 0; t < splits.size(); t++) {
					for (i = 1; i < splits[0].size(); i++) {
						inter.lb = (double)(splits[t][i].lb - 1) / 100;
						inter.rb = (double)splits[t][i].rb / 100;
						splitsJini[t].push_back(inter);
					}
				}
				jiniIndices = new dt::matrix(strats.times.size(), splits[0].size());
				for (t = 0; t < splits.size(); t++) {
					std::cout << "Period " << t + 1 << " out of " << splits.size() << std::endl;
					lcdsClasses = lcds[t].split(splitsJini[t]);
					for (i = 0; i < splits[0].size(); i++) {
						std::cout << "\tClass " << i + 1 << " out of " << splits[0].size() << std::endl;
						(*jiniIndices)(t, i) = lcdsClasses[i].jiniIndex();
					}
				}
				reportOnJini.open(path + "Proportions 3\\JiniInds.txt");
				for (i = 0; i < splits[0].size(); i++) {
					reportOnJini << (*jiniIndices)(0, i);
					for (t = 1; t < splits.size(); t++) {
						reportOnJini << ',' << (*jiniIndices)(t, i);
					}
					reportOnJini << std::endl;
				}
				delete jiniIndices;
				reportOnJini.close();
				std::cout << "Computing Jini indices: Complete" << std::endl;
				*/
				// Get proportions
				for (i = 0; i < stratsClasses.Xs.size(); i++) {
					std::cout << "Class " << i + 1 << std::endl;
					report << "Class " << i + 1 << std::endl;
					groups = dataload::extractGroups(stratsClasses.P, stratsClasses.Xs[i], 
						stratsClasses.gNames, stratsClasses.times, path + "Proportions 3\\ComGroups.txt");
					
					// Omega for all goods
					lOmegaMax = axioms::HARPLogIrInd(stratsClasses.P, stratsClasses.Xs[i], llb, lrb);
					if (lOmegaMax < llb) {
						lOmegaMax = llb;
					}
					std::cout << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;
					report << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;
		
					// Omegas for groups
					for (g = 0; g < groups.size(); g++) {
						if (groups[g].cNames.size() > 1) {
							lOmega = axioms::HARPLogIrInd(groups[g].P, groups[g].X, llb, lrb);
							if (lOmega < llb) {
								lOmega = llb;
							}
							std::cout << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
							report << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
							lOmegaMax = std::max(lOmegaMax, lOmega);
						}
					}
					std::cout << std::endl;
					report << std::endl;
		
					// PHKD indices
					indAll = axioms::HARP(stratsClasses.P, stratsClasses.Xs[i], lOmegaMax);
					reportOnIndicesQ.open(path + "Proportions 3\\QS" + std::to_string(i) + ".txt");
					reportOnIndicesF.open(path + "Proportions 3\\FS" + std::to_string(i) + ".txt");
					reportOnIndicesQ << indAll.p[0];
					reportOnIndicesF << indAll.x[0];
					for (t = 1; t < indAll.p.size(); t++) {
						reportOnIndicesQ << ' ' << indAll.p[t];
						reportOnIndicesF << ' ' << indAll.x[t];
					}
					reportOnIndicesQ << std::endl;
					reportOnIndicesF << std::endl;
					for (g = 0; g < groups.size(); g++) {
						if (groups[g].P[0].size() == 1) {
							tmpMat = axioms::transpose(groups[g].P);
							ind.p = tmpMat[0];
							tmpMat = axioms::transpose(groups[g].X);
							ind.x = tmpMat[0];
							ind.omega = 0;
						}
						else {
							ind = axioms::HARP(groups[g].P, groups[g].X, lOmegaMax);
						}
						reportOnIndicesQ << ind.p[0];
						reportOnIndicesF << ind.x[0];
						for (t = 1; t < ind.p.size(); t++) {
							reportOnIndicesQ << ' ' << ind.p[t];
							reportOnIndicesF << ' ' << ind.x[t];
						}
						reportOnIndicesQ << std::endl;
						reportOnIndicesF << std::endl;
					}
					reportOnIndicesQ.close();
					reportOnIndicesF.close();
				}
				report.close();
				std::cout << "OK" << std::endl;
			}
			else if (c == 'n') {
				break;
			}
			else {
				continue;
			}
		}
	}
	catch (std::string err) {
		if (jiniIndices != nullptr) {
			delete jiniIndices;
		}
		std::cout << "ERROR: " << err << std::endl;
	}
	catch (dt::error err) {
		if (jiniIndices != nullptr) {
			delete jiniIndices;
		}
		std::cout << "ERROR: " << err.text << std::endl;
	}

	return 0;
}

int proportions2SequenceInfWithGroups() {
	char c;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::string pathToStrataConsumption = path + "Hundred parts by income\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	double lOmegaMax, lOmega;
	std::ofstream report, reportF, reportQ;
	std::vector<axioms::group> groups;
	std::vector<dt::tradestat> tradestats;
	dt::splitData splits; // Data on how to split strata into two classes
	int i, g, t;
	int S = 100; // Number of stata
	std::vector<std::string> namesRus, times;
	dt::indices indAll, ind;
	std::vector<double> expAll, exp;
	dt::strata strats;

	double sum;

	try {
		std::cout << "This is proportions2SequenceInf" << std::endl;
		//Load data to strata class.
		strats.P = dataload::loadDataTable(path + "pallForStata.csv"); // Price data
		strats.gNames = dataload::loadNames(path + "comNamesRus.csv"); // Names of commodities
		strats.times = dataload::loadNames(path + "perNums.csv"); // Names of periods
		// Load strata consumption data
		std::cout << "Loading consumption data for " << S << " strata..." << std::endl;
		for (i = 1; i <= S; i++) {
			strats.Xs.push_back(dataload::loadDataTable(pathToStrataConsumption + "Part" + std::to_string(i) + ".csv"));
			strats.sNames.push_back("");
		}
		std::cout << "Loading consumption data for " << S << " strata: Complete" << std::endl;		
		std::cout << "Python:OK" << std::endl;
		while (true) {
			std::cout << "proportions2SequenceInf?" << std::endl;
			std::cin >> c;
			if (c == 'y') {
				report.open(path + "Proportions 2 sequence\\report.txt");
				// Join strata into three classes
				splits = dataload::loadClassSplitDataNew2(path + "Proportions 2 sequence\\Splits.txt");
				tradestats = strats.getClasses(splits);
				
				// Get omegas
				for (i = 0; i < tradestats.size(); i++) {
					std::cout << "Class " << i + 1 << std::endl;
					report << "Class " << i + 1 << std::endl;
					groups = dataload::extractGroups(tradestats[i].P, tradestats[i].X, 
						tradestats[i].gNames, tradestats[i].times, path + "Proportions 2 sequence\\ComGroups.txt");
					
					// Omega for all goods
					lOmegaMax = axioms::HARPLogIrInd(tradestats[i].P, tradestats[i].X, llb, lrb);
					if (lOmegaMax < llb) {
						lOmegaMax = llb;
					}
					std::cout << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;
					report << "Omega for all goods = " << std::exp(lOmegaMax) << std::endl;
		
					// Omegas for groups
					for (g = 0; g < groups.size(); g++) {
						if (groups[g].cNames.size() > 0) {
							lOmega = axioms::HARPLogIrInd(groups[g].P, groups[g].X, llb, lrb);
							if (lOmega < llb) {
								lOmega = llb;
							}
							std::cout << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
							report << "Omega for group " << g + 1 << " = " << std::exp(lOmega) << std::endl;
							lOmegaMax = std::max(lOmegaMax, lOmega);
						}
					}
					std::cout << std::endl;
					report << std::endl;

					// Open files for saving demand indices, print times
					reportF.open(path + "Proportions 2 sequence\\F" + std::to_string(i + 1) + ".txt");
					reportQ.open(path + "Proportions 2 sequence\\Q" + std::to_string(i + 1) + ".txt");
					reportF << tradestats[i].times[0];
					reportQ << tradestats[i].times[0];
					for (t = 1; t < tradestats[i].times.size(); t++) {
						reportF << ' ' << tradestats[i].times[t];
						reportQ << ' ' << tradestats[i].times[t];
					}
					reportF << std::endl;
					reportQ << std::endl;

					// Print demand indices for all goods
					indAll = axioms::HARP(tradestats[i].P, tradestats[i].X, lOmegaMax);
					reportF << indAll.x[0];
					reportQ << indAll.p[0];
					for (t = 1; t < tradestats[i].times.size(); t++) {
						reportF << ' ' << indAll.x[t];
						reportQ << ' ' << indAll.p[t];
					}
					reportF << std::endl;
					reportQ << std::endl;

					// Print expenditures for trade groups
					for (g = 0; g < groups.size(); g++) {
						ind = axioms::HARP(groups[g].P, groups[g].X, lOmegaMax);
						reportF << ind.x[0];
						reportQ << ind.p[0];
						for (t = 1; t < tradestats[i].times.size(); t++) {
							reportF << ' ' << ind.x[t];
							reportQ << ' ' << ind.p[t];
						}
						reportF << std::endl;
						reportQ << std::endl;
					}
					reportF.close();
					reportQ.close();
				}
				report.close();
				std::cout << "OK" << std::endl;
			}
			else if (c == 'n') {
				break;
			}
			else {
				continue;
			}
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}
	catch (dt::error err) {
		std::cout << "ERROR: " << err.text << std::endl;
	}

	return 0;
}

int proportions2SequenceInf() {
	char c;
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	std::string pathToStrataConsumption = path + "Hundred parts by income\\";
	double lrb = std::log(2);
	double llb = std::log(0.5);
	double lOmega;
	std::ofstream report, reportF, reportQ;
	std::vector<dt::tradestat> tradestats;
	dt::splitData splits; // Data on how to split strata into two classes
	int i, g, t;
	int S = 100; // Number of stata
	dt::indices indAll;
	dt::strata strats;

	double sum;

	try {
		//Load data to strata class.
		strats.P = dataload::loadDataTable(path + "pallForStata.csv"); // Price data
		strats.gNames = dataload::loadNames(path + "comNamesRus.csv"); // Names of commodities
		strats.times = dataload::loadNames(path + "perNums.csv"); // Names of periods
		// Load strata consumption data
		std::cout << "Loading consumption data for " << S << " strata..." << std::endl;
		for (i = 1; i <= S; i++) {
			strats.Xs.push_back(dataload::loadDataTable(pathToStrataConsumption + "Part" + std::to_string(i) + ".csv"));
			strats.sNames.push_back("");
		}
		std::cout << "Loading consumption data for " << S << " strata: Complete" << std::endl;

		while (true) {
			std::cout << "proportions2SequenceInf?" << std::endl;
			std::cin >> c;
			if (c == 'y') {
				report.open(path + "Proportions 2 sequence PhD\\report.txt");
				// Split strata into two classes
				splits = dataload::loadClassSplitDataNew2(path + "Proportions 2 sequence PhD\\Splits.txt"); // Load split scheme
				tradestats = strats.getClasses(splits); // Split strata into two classes according to the split scheme
		
				// Get omegas
				for (i = 0; i < tradestats.size(); i++) {
					std::cout << "Class " << i + 1 << std::endl;
					report << "Class " << i + 1 << std::endl;
			
					// Compute irrationality indices
					lOmega = axioms::HARPLogIrInd(tradestats[i].P, tradestats[i].X, llb, lrb);
					if (lOmega < llb) {
						lOmega = llb;
					}
					std::cout << "Omega for all goods = " << std::exp(lOmega) << std::endl;
					report << "Omega for all goods = " << std::exp(lOmega) << std::endl;

					std::cout << std::endl;
					report << std::endl;

					// Open files for saving indices, print times
					reportF.open(path + "Proportions 2 sequence PhD\\F" + std::to_string(i + 1) + ".txt");
					reportQ.open(path + "Proportions 2 sequence PhD\\Q" + std::to_string(i + 1) + ".txt");
					reportF << tradestats[i].times[0];
					reportQ << tradestats[i].times[0];
					for (t = 1; t < tradestats[i].times.size(); t++) {
						reportF << ' ' << tradestats[i].times[t];
						reportQ << ' ' << tradestats[i].times[t];
					}
					reportF << std::endl;
					reportQ << std::endl;

					// Print indices for all goods
					indAll = axioms::HARP(tradestats[i].P, tradestats[i].X, lOmega);
					reportF << indAll.x[0];
					reportQ << indAll.p[0];
					for (t = 1; t < tradestats[i].times.size(); t++) {
						reportF << ' ' << indAll.x[t];
						reportQ << ' ' << indAll.p[t];
					}
					reportF << std::endl;
					reportQ << std::endl;

					reportF.close();
					reportQ.close();
				}
				report.close();
				std::cout << "OK" << std::endl;
			}
			else if (c == 'n') {
				break;
			}
			else {
				continue;
			}
		}
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	return 0;
}

int jini() {
	std::string path = "D:\\PhD\\Work with BCC2008 data\\";
	dt::lorenzData lcd;
	int i;
	std::ofstream report;
	double jiniInd;

	try {
		report.open(path + "Jini\\report.txt");
		if (!report) {
			throw "Unable to open file " + path + "Jini\\report.txt";
		}
		for (i = 75; i <= 99; i++) {
			std::cout << "Year 19" << i << std::endl;
			lcd = dataload::loadLorenzCurveData(path + "Jini\\Year" + strfuns::int2str(i) + ".csv");
			jiniInd = axioms::giniIndex(lcd);
			report << jiniInd << std::endl;
		}
		report.close();
	}
	catch (std::string err) {
		std::cout << "ERROR: " << err << std::endl;
	}

	return 0;
}

bool is_in(int elem, const std::vector<int>& lst)
{
	for (int i = 0; i < lst.size(); i++)
	{
		if (lst[i] == elem)
		{
			return true;
		}
	}
	return false;
}

void exchange_lambdas_one_line(const std::string& line, std::ofstream& fout)
{
	std::vector<std::string> ss = strfuns::split(line, ',');
	std::string date = ss[0];
	std::string cur;
	int n_cur = 16;
	int n_cur_left;
	double max_prod, min_prod, next_prod;
	int i_max, i_min, j_max, j_min;
	dt::matrix A(n_cur, n_cur);
	dt::matrix A_cut;
	std::vector<double> row;
	int ind;
	int i, j;
	std::vector<int> bad_cur_inds;
	double r;
	std::vector<double> lambdas;

	ind = 1;

	// Get the matrix of cross-rates 
	for (i = 0; i < n_cur; i++)
	{
		for (j = 0; j < n_cur; j++)
		{
			if (i == j)
			{
				A[i][i] = 1.0;
			}
			else
			{
				cur = ss[ind++];
				if (cur == "")
				{
					if (!is_in(j, bad_cur_inds))
					{
						if (j > i)
						{
							bad_cur_inds.push_back(j);
						}
					}
					A[i][j] = -1.0;
				}
				else
				{
					A[i][j] = std::atof(cur.c_str());
				}
			}
		}
	}

	std::cout << date << std::endl;
	// Test symmetry of error
	std::cout << "Test symmetry of error" << std::endl;
	for (i = 0; i < n_cur; i++)
	{
		for (j = i + 1; j < n_cur; j++)
		{
			if (A[i][j] * A[j][i] < 0)
			{
				std::cout << "A[" << i << "][" << j << "] = " << A[i][j] 
				<< " A[" << j << "][" << i << "] = " << A[j][i] << std::endl;
				i = n_cur;
				break;
			}
		}
	}
	
	// Test reciprocal property
	std::cout << "Test reciprocal property" << std::endl;
	max_prod = A[0][1] * A[1][0];
	min_prod = max_prod;
	j_min = 1;
	j_max = 1;
	i_min = 0;
	i_max = 0;
	for (i = 0; i < n_cur; i++)
	{
		for (j = i + 1; j < n_cur; j++)
		{
			if ((A[i][j] == -1.0) || (A[j][i] == -1.0))
			{
				continue;
			}
			else
			{
				next_prod = A[i][j] * A[j][i];
				if (next_prod < min_prod)
				{
					min_prod = next_prod;
					i_min = i;
					j_min = j;
				}
				if (next_prod > max_prod)
				{
					max_prod = next_prod;
					i_max = i;
					j_max = j;
				}
			}
		}
	}
	std::cout << "max_prod = " << max_prod << " attained on i = " << i_max << " and j = " << j_max << std::endl;
	std::cout << "min_prod = " << min_prod << " attained on i = " << i_min << " and j = " << j_min << std::endl;

	// Calculate lambdas
	// Prepare the matrix
	n_cur_left = n_cur - bad_cur_inds.size();
	for (i = 0; i < n_cur_left; i++)
	{
		row.push_back(0.0);
	}
	for (i = 0; i < n_cur; i++)
	{
		if (is_in(i, bad_cur_inds))
		{
			continue;
		}
		else
		{
			row.clear();
			for (j = 0; j < n_cur; j++)
			{
				if (!is_in(j, bad_cur_inds))
				{
					row.push_back(std::log(A[i][j]));
				}
			}
			A_cut.push_back(row);
		}
	}

	// Irrationality index
	r = axioms::hav_sys_min_r(A_cut, -0.5, 15, 0.00000001);

	// Lambdas
	lambdas = axioms::hav_sys_solve(A_cut, r);
	for (i = 1; i < n_cur_left; i++)
	{
		lambdas[i] = std::exp(lambdas[i] - lambdas[0]);
	}
	lambdas[0] = 1.0;

	// Output
	ind = 0;
	fout << date;
	for (i = 0; i < n_cur; i++)
	{
		fout << ',';
		if (!is_in(i, bad_cur_inds))
		{
			fout << lambdas[ind++];
		}
	}
	fout << std::endl;
}

void exchange_lambdas()
{
	std::ifstream fin;
	std::ofstream fout;
	std::string wd_path = "D:\\PhD\\Shanghai 2015\\New data from Misha\\Nice new data\\";
	int i;
	std::string s;
	std::vector<std::string> ss;

	fin.open(wd_path + "cur_data_matched_merged_interp_trans.csv");
	fout.open(wd_path + "exchange_lambdas_usd.csv");
	std::getline(fin, s);
	ss = strfuns::split(s, ',');
	fout << "date,USD";
	for (i = 1; i < 16; i++)
	{
		fout << ',' << ss[i][3] << ss[i][4] << ss[i][5];
	}
	fout << std::endl;
	while(true)
	{
		std::getline(fin, s);
		if (fin.eof())
		{
			break;
		}
		exchange_lambdas_one_line(s, fout);
	}
	fout.close();
	fin.close();
}

int main(int argc, char** argv)
{
	char c;

	proportions2SequenceInf();

	std::cout << std::endl << "Work complete. Enter any key with letter." << std::endl;
	std::cin >> c;
	
	return 0;
}