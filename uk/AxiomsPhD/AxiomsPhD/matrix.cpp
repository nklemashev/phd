#include "dt.hpp"

dt::matrix::matrix(int nr, int nc) {
	if ((nr < 0) || (nc < 0)) {
		throw dt::error("Matrix size data contain negative value");
	}
	std::vector<double> row;
	int i;

	for (i = 0; i < nc; i++) {
		row.push_back(0.0);
	}
	for (i = 0; i < nr; i++) {
		this->push_back(row);
	}
}

double dt::matrix::operator()(int row, int col) const {
	return this->operator[](row).operator[](col);
}

double& dt::matrix::operator()(int row, int col){
	return this->operator[](row).operator[](col);
}