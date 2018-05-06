#include "dt.hpp"
#include "axioms.hpp"

dt::tradestat::tradestat(const matrix& inP, const matrix& inX, const std::vector<std::string>& inTimes,
			const std::vector<std::string>& inNames, const std::string& inName) {
	P = inP;
	X = inX;
	times = inTimes;
	name = inName;
	logPM = axioms::logPaascheMat(axioms::dotProds(P, X));
}
/*
dt::tradestat dt::tradestat::select(const dt::intervals& gInts, const dt::intervals& tInts, const std::string& inName) {
	// gInts -- intervals for goods
	// tInts -- intervals for times
	// name -- name of the resulting trade statistics' group

	dt::tradestat result;
	int i, iInt, t, tInt;
	int N = gds.size();
	int T = gds[0]->p.size();
	dt::good* g;

	try {
		result.name = inName;
		for (tInt = 0; tInt < tInts.size(); tInt++) {
			for (t = tInts[tInt].lb - 1; t < tInts[tInt].rb; t++) {
				result.times.push_back(times[t]);
			}
		}
		for (iInt = 0; iInt < gInts.size(); iInt++) {
			for (i = gInts[iInt].lb - 1; i < gInts[iInt].rb; i++) {
				g = new dt::good();
				g->name = gds[i]->name;
				for (tInt = 0; tInt < tInts.size(); tInt++) {
					for (t = tInts[tInt].lb - 1; t < tInts[tInt].rb; t++) {
						g->p.push_back(gds[i]->p[t]);
						g->x.push_back(gds[i]->x[t]);
					}
				}
				result.gds.push_back(g);
			}
		}
		result.logPM = axioms::logPaascheMat(axioms::dotProds(result.gds));
	}
	catch(...) {
		if ((g != nullptr) && (g != result.gds[result.gds.size() - 1])) {
			delete g;
		}
		throw;
	}

	return result;
}

dt::tradestat dt::tradestat::selectGoods(const intervals& gInts, const std::string& inName) {
	// gInts -- intervals for goods
	// name -- name of the resulting trade statistics' group

	dt::intervals tInts;
	tInts.push_back(dt::interval(1, times.size()));
	return select(gInts, tInts, inName);
}

dt::tradestat dt::tradestat::selectTimes(const intervals& tInts, const std::string& inName) {
	// tInts -- intervals for times
	// name -- name of the resulting trade statistics' group

	dt::intervals gInts;
	gInts.push_back(dt::interval(1, gds.size()));
	return select(gInts, tInts, inName);
}
*/