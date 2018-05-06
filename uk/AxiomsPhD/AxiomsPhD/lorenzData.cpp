#include "dt.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include "axioms.hpp"

std::vector<dt::lorenzData> dt::lorenzData::split(const dt::intervals& splits) const {
	// splits[i] = [a, b]
	// a, b -- percentages of income for group i.
	int i, j, lInd, rInd, k;
	std::vector<dt::lorenzData> result;
	dt::lorenzData lcd;
	double xDiff, yDiff;
	rInd = 0;

	if(!splits.checkCosistency()) {
		throw dt::error("lorenzData::Incorrect splits data: the intervals " + splits.toString() + " are not sorted");
	}
	// Additional check -- intervals boundaries should be in [0, 1]
	for (i = 0; i < splits.size(); i++) {
		if ((splits[i].lb < 0) || (splits[i].rb > 1)) {
			throw dt::error("lorenzData::Incorrect splits data: the interval " + splits[i].toString() + " is not the subset of [0, 1]");
		}
	}
	for (i = 0; i < splits.size(); i++) {
		lcd.x.clear();
		lcd.y.clear();
		for (j = rInd; j < this->x.size(); j++) {
			if (this->x[j] >= splits[i].lb) {
				lInd = j;
				break;
			}
		}
		for (j = lInd; j < this->x.size(); j++) {
			if (this->x[j] > splits[i].rb) {
				rInd = j;
				break;
			}
		}
		xDiff = this->x[rInd] - this->x[lInd];
		if (std::abs(xDiff) < axioms::zero) {
			throw dt::error("lorenzData::Incorrect splits data: interval " + splits[i].toString() + " is too narrow.");
		}
		yDiff = this->y[rInd] - this->y[lInd];
		for (j = lInd; j <= rInd; j++) {
			lcd.x.push_back((this->x[j] - this->x[lInd]) / xDiff);
			lcd.y.push_back((this->y[j] - this->y[lInd]) / yDiff);
		}
		result.push_back(lcd);
	}
	return result;
}

std::vector<dt::lorenzData> dt::lorenzData::split(const std::vector<simpleDInterval>& splits) const {
	// splits[i] = [a, b]
	// a, b -- percentages of income for group i.
	int i, j, lInd, rInd, k;
	std::vector<dt::lorenzData> result;
	dt::lorenzData lcd;
	double xDiff, yDiff;
	rInd = 0;

	/*
	if(!splits.checkCosistency()) {
		throw dt::error("lorenzData::Incorrect splits data: the intervals " + splits.toString() + " are not sorted");
	}
	*/
	// Additional check -- intervals boundaries should be in [0, 1]
	for (i = 0; i < splits.size(); i++) {
		if ((splits[i].lb < 0) || (splits[i].rb > 1)) {
			throw dt::error("lorenzData::Incorrect splits data: the interval " + splits[i].toString() + " is not the subset of [0, 1]");
		}
	}
	for (i = 0; i < splits.size(); i++) {
		lcd.x.clear();
		lcd.y.clear();
		for (j = rInd; j < this->x.size(); j++) {
			if (this->x[j] >= splits[i].lb) {
				lInd = j;
				break;
			}
		}
		for (j = lInd; j < this->x.size(); j++) {
			if (this->x[j] > splits[i].rb) {
				rInd = j;
				break;
			}
		}
		if (j == this->x.size()) {
			rInd = this->x.size() - 1;
		}
		xDiff = this->x[rInd] - this->x[lInd];
		if (std::abs(xDiff) < axioms::zero) {
			throw dt::error("lorenzData::Incorrect splits data: interval " + splits[i].toString() + " is too narrow.");
		}
		yDiff = this->y[rInd] - this->y[lInd];
		for (j = lInd; j <= rInd; j++) {
			lcd.x.push_back((this->x[j] - this->x[lInd]) / xDiff);
			lcd.y.push_back((this->y[j] - this->y[lInd]) / yDiff);
		}
		result.push_back(lcd);
	}
	return result;
}

double dt::lorenzData::giniIndex() const {
	double sum;
	int i;

	sum = 0.0;
	for (i = 2; i < x.size(); i++) {
		sum += (x[i] - x[i - 1]) * y[i];
	}
	return (0.5 - sum) * 2;
}