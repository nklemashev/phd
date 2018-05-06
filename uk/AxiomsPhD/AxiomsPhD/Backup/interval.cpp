#include "dt.hpp"
#include <algorithm>

bool dt::interval::operator> (dt::interval right) {
	if (this->rb > right.lb) {
		return true;
	}
	else if (right.rb >= this->lb) {
		return false;
	}
	else {
		throw dt::error("interval::operator>: Unable to compare " + this->toString() + " with " + right.toString());
	}
}

bool dt::interval::operator>= (dt::interval right) {
	return ((*this == right) || (*this > right));
}

bool dt::interval::operator< (dt::interval right) {
	return !(*this >= right);
}

bool dt::interval::operator<= (dt::interval right) {
	return !(*this > right);
}