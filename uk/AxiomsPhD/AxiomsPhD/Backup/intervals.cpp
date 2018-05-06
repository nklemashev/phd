#include "dt.hpp"
#include <algorithm>
#include "strfuns.hpp"

void dt::intervals::consistency() {
	dt::intervals::iterator i;
	
	std::sort(this->begin(), this->end());
	for (i = this->begin(); i != this->end() - 1; i++) {
		if ((i->rb + 1) == (*(i + 1)).lb) {
			i->rb = (*(i + 1)).rb;
			this->erase(i + 1);
		}
	}
}

std::string dt::intervals::toString() {
	std::string result;
	dt::intervals::iterator i;

	if (this->size() == 0) {
		return "";
	}
	result = this->begin()->toString();
	for (i = this->begin() + 1; i != this->end(); i++) {
		result += ", " + i->toString();
	}
	return result;
}