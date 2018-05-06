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
	struct matrix : public std::vector<std::vector<double>> {
		matrix(int nr = 0, int nc = 0);
		double operator()(int row, int col) const;
		double& operator()(int row, int col);
	};
	struct interval {
		int lb, rb;

		interval(int inLb, int inRb);
		std::string toString() const{
			return "[" + strfuns::int2str(lb) + ", " + strfuns::int2str(rb) + "]";
		}
		bool operator==(interval right) const{
			return ((this->lb == right.lb) && (this->rb == right.rb));
		}
		bool operator!=(interval right) const {
			return !(*this==right);
		}
		bool operator> (interval right) const;
		bool operator>= (interval right) const;
		bool operator< (interval right) const;
		bool operator<= (interval right) const;
	};
	const interval empty_interval(0, 0);
	struct interval_error : public error {
		int code;
		// 1 -- unable to compare two intervals.
		// 2 -- left boundary is greater than the right one.
		interval_error(const std::string& inText, int inCode = 0) : error(inText), code(inCode) {}
	};
	struct intervals : public std::vector<interval> {
		void compact(); // Converts two intervals [a, b], [b + 1, c] into one [a, c] if possible.
		std::string toString() const;
		bool checkCosistency() const;
	};
	struct splitLine {
		std::string period;
		intervals inters;

		splitLine() : period(""), inters() {}
	};
	typedef std::vector<splitLine> splitData;
	struct good {
		std::vector<double> p, x;
		std::string name;

		good() : p(0), x(0), name("") {
		}
		good(const std::vector<double>& inP, const std::vector<double>& inX, const std::string inName) {
			p = inP;
			x = inX;
			name = inName;
		}
	};
	class goods {
	private:
		std::vector<good*> data;
	public:
		goods(int n = 0) : data(n) {
		}
		void push_back(good* inGood) {
			data.push_back(inGood);
		}
		good* operator[](int ind) {
			return data[ind];
		}
		const good* operator[](int ind) const {
			return data[ind];
		}
		int size() const {
			return data.size();
		}
		/*void erase(int ind);*/
		~goods() {
			for (int i = 0; i < data.size(); i++) {
				delete(data[i]);
			}
		}
	};
	struct tradestat {
		matrix P, X;
		std::vector<std::string> times;
		std::string name;
		std::vector<std::string> gNames;
		matrix logPM;

		tradestat() : P(0), X(0), gNames(0), times(0), name(""), logPM(0) {
		}
		tradestat(const std::vector<goods*>& gds, std::vector<std::string> inTimes, std::string inName = ""); 
		tradestat(const matrix& inP, const matrix& inX, const std::vector<std::string>& inTimes,
			const std::vector<std::string>& inNames, const std::string& inName = "");
		tradestat select(const intervals& gInts, const intervals& tInts, const std::string& inName = ""); // Select goods (gInts) and periods (tInds)
		tradestat selectGoods(const intervals& gInts, const std::string& inName = "");
		tradestat selectTimes(const intervals& tInts, const std::string& inName = "");
	};
	struct strata {
		matrix P;
		std::vector<matrix> Xs;
		std::vector<std::string> times;
		std::vector<std::string> sNames;
		std::vector<std::string> gNames;
		std::string name;

		strata() : P(0), Xs(0), times(0), sNames(0), gNames(0), name("") {};
		strata(const matrix& inP, const std::vector<matrix>& inXs, 
			const std::vector<std::string>& inTimes, const std::vector<std::string>& inSNames,
			const std::vector<std::string>& inGNames, const std::string& inName = "");
		strata combineStrata(std::vector<intervals> splits, const std::vector<std::string>& sNames, const std::string& name = "") const;
		strata combineStrata(std::vector<intervals> splits, const std::string& name = "") const;
		std::vector<tradestat> getClasses(const splitData& splits,
			const std::vector<std::string>& cNames) const;
		std::vector<tradestat> getClasses(const dt::splitData& splits) const;
	};
	struct indices : public good
	{
		double omega;
	};
	struct simpleDInterval {
		double lb, rb;
		simpleDInterval(double inLb, double inRb) : lb(inLb), rb(inRb) {}
		std::string toString() const {
			return "[" + std::to_string(lb) + ", " + std::to_string(rb) + "]";
		}
	};
	struct lorenzData {
		std::vector<double> x;
		std::vector<double> y;

		lorenzData() : x(0), y(0) {}
		std::vector<lorenzData> split(const intervals& splits) const;
		std::vector<lorenzData> split(const std::vector<simpleDInterval>& splits) const;
		double giniIndex() const;
	};
}

#endif