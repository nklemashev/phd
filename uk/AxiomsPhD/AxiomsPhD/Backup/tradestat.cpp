#include "dt.hpp"
#include "axioms.hpp"

dt::tradestat::tradestat(const std::vector<good>& inGoods, const std::vector<std::string>& inTimes,
			const std::string& inName) {
	goods = inGoods;
	times = inTimes;
	name = inName;
	logPM = axioms::logPaascheMat(axioms::dotProds(inGoods));
}

dt::tradestat dt::tradestat::select(const dt::intervals& gInts, const dt::intervals& tInts, const std::string& inName) {
	// gInts -- intervals for goods
	// tInts -- intervals for times
	// name -- name of the resulting trade statistics' group

	dt::tradestat result;
	int i, iInt, t, tInt;
	int N = goods.size();
	int T = goods[0].p.size();
	dt::good g;

	result.name = inName;
	for (tInt = 0; tInt < tInts.size(); tInt++) {
		for (t = tInts[tInt].lb - 1; t < tInts[tInt].rb; t++) {
			result.times.push_back(times[t]);
		}
	}
	for (iInt = 0; iInt < gInts.size(); iInt++) {
		for (i = gInts[iInt].lb - 1; i < gInts[iInt].rb; i++) {
			g.p.clear();
			g.x.clear();
			g.name = goods[i].name;
			for (tInt = 0; tInt < tInts.size(); tInt++) {
				for (t = tInts[tInt].lb - 1; t < tInts[tInt].rb; t++) {
					g.p.push_back(goods[i].p[t]);
					g.x.push_back(goods[i].p[t]);
				}
			}
			result.goods.push_back(g);
		}
	}
	result.logPM = axioms::logPaascheMat(axioms::dotProds(result.goods);

	return result;
}

dt::tradestat dt::tradestat::selectGoods(const intervals& gInts, const std::string& inName) {
	// gInts -- intervals for goods
	// name -- name of the resulting trade statistics' group

	dt::tradestat result;
	int i, iInt
	int N = goods.size();

	result.name = inName;
	result.times = times;
	for (iInt = 0; iInt < gInts.size(); iInt++) {
		for (i = gInts[iInt].lb - 1; i < gInts[iInt].rb; i++) {
			result.goods.push_back(goods[i]);
		}
	}
	result.logPM = axioms::logPaascheMat(axioms::dotProds(result.goods);

	return result;
}

dt::tradestat dt::tradestat::selectTimes(const intervals& tInts, const std::string& inName) {
	// tInts -- intervals for times
	// name -- name of the resulting trade statistics' group

	dt::intervals gInts;
	gInts.push_back(dt::interval(1, goods.size()));
	return select(gInts, tInts, inName);
}