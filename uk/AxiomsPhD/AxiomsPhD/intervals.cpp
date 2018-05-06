#include "dt.hpp"
#include <algorithm>
#include "strfuns.hpp"

void dt::intervals::compact() {
	dt::intervals::iterator i;
	
	std::sort(this->begin(), this->end());
	for (i = this->begin(); i != this->end() - 1; i++) {
		if ((i->rb + 1) == (*(i + 1)).lb) {
			i->rb = (*(i + 1)).rb;
			this->erase(i + 1);
		}
	}
}

std::string dt::intervals::toString() const {
	std::string result;
	dt::intervals::const_iterator i;

	if (this->size() == 0) {
		return "";
	}
	result = this->begin()->toString();
	for (i = this->begin() + 1; i != this->end(); i++) {
		result += ", " + i->toString();
	}
	return result;
}

bool dt::intervals::checkCosistency() const {
	// Checks if intervals are comparable and are sorted
	int i;

	try {
		for (i = 1; i < this->size(); i++) {
			if (this->operator[](i) <= this->operator[](i - 1)) {
				return false;
			}
		}
	}
	catch (dt::interval_error err) {
		if (err.code == 1) {
			return false;
		}
		else {
			throw dt::error(err.text);
		}
	}
	catch (...) {
		throw;
	}
}