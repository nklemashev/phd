#include "dt.hpp"
#include <algorithm>

dt::interval::interval(int inLb, int inRb) : lb(inLb), rb(inRb) {
	if (inLb > inRb) {
		throw interval_error("Interval: left boubary is greater than the right one", 2);
	}
}

bool dt::interval::operator> (dt::interval right) const{
	if (this->rb > right.lb) {
		return true;
	}
	else if (right.rb >= this->lb) {
		return false;
	}
	else {
		throw dt::interval_error("interval::operator>: Unable to compare " + this->toString() + " with " + right.toString(), 1);
	}
}

bool dt::interval::operator>= (dt::interval right) const{
	return ((*this == right) || (*this > right));
}

bool dt::interval::operator< (dt::interval right) const{
	return !(*this >= right);
}

bool dt::interval::operator<= (dt::interval right) const{
	return !(*this > right);
}