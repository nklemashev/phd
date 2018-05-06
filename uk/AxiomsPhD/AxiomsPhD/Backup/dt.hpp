#ifndef DT_HPP
#define DT_HPP

#include <vector>
#include <string>
#include "strfuns.hpp"

namespace dt {
	struct error {
		std::string text;

		error(const std::string& inText) {
			text = inText;
		}
	};
	struct interval {
		int lb, rb;

		interval(int inLb, int inRb) {
			lb = inLb;
			rb = inRb;
		}
		std::string toString() {
			return "[" + strfuns::int2str(lb) + ", " + strfuns::int2str(rb) + "]";
		}
		bool operator==(interval right) {
			return ((this->lb == right.lb) && (this->rb == right.rb));
		}
		bool operator> (interval right);
		bool operator>= (interval right);
		bool operator< (interval right);
		bool operator<= (interval right);
	};
	struct intervals : std::vector<interval> {
		void consistency(); // Converts two intervals [a, b], [b + 1, c] into one [a, c] if possible.
		std::string toString();
	};
	struct good {
		std::vector<double> p, x;
		std::string name;

		good(const std::vector<double>& inP, const std::vector<double>& inX, const std::string inName) {
			p = inP;
			x = inX;
			name = inName;
		}
	};
	typedef std::vector<good> goods;
	struct tradestat {
		std::vector<good> goods;
		std::vector<std::string> times;
		std::string name;
		std::vector<std::vector<double>> logPM;

		tradestat(const std::vector<good>& inGoods, const std::vector<std::string>& inTimes,
			const std::string& inName = "");
		tradestat select(const intervals& gInts, const intervals& tInts, const std::string& inName = ""); // Select goods (gInts) and periods (tInds)
		tradestat selectGoods(const intervals& gInts, const std::string& inName = "");
		tradestat selectTimes(const intervals& tInts, const std::string& inName = "");
	};
	typedef std::vector<std::vector<double>> matrix; // vector of rows
	struct indices : public good
	{
		double omega;
	};
	struct lorenzData {
		std::vector<double> x;
		std::vector<double> y;
	};
}

#endif