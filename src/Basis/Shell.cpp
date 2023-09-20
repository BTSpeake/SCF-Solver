#include "Shell.h"

#include <cmath>
#include <numbers>

Shell::Shell(int i, int l, double px, double py, double pz) 
	: _atom(i), _l(l) {
	_position = Vector3(px, py, pz);
}

Shell::~Shell() {}

void Shell::addExponent(const double exp, const double cc) {
	_exponents.push_back(exp);
	_coefficients.push_back(cc);
}

int Shell::getAngularMomentum() const {
	return _l;
}

const std::vector<double>& Shell::getExponents() const {
	return _exponents;
}

const double& Shell::getExponent(int i) const {
	return _exponents[i];
}

int Shell::nExponents() const {
	return _exponents.size();
}

const std::vector<double>& Shell::getCoefficients() const {
	return _coefficients;
}

const double& Shell::getCoefficient(int i) const {
	return _coefficients[i];
}

int Shell::nCoefficients() const {
	return _coefficients.size();
}

int Shell::getSphericalStart() const {
	return _si;
}

int Shell::nSpherical() const {
	return (2 * _l) + 1;
}

int Shell::getCartesianStart() const {
	return _ci;
}

int Shell::nCartesian() const {
	return (_l + 1) * (_l + 2) / 2;
}

const Vector3& Shell::getPosition() const {
	return _position;
}

int Shell::getlmn(int k, int i) const {
	return _lmn[(k * 3) + i];
}
const double& Shell::getNormal(int i) const {
	return _normal[i];
}

void Shell::calculateProperties(int si, int ci) {
	_si = si; 
	_ci = ci; 
	calculateGaussianNormal();
	calculateCanonicallmn();
}

void Shell::calculateGaussianNormal() {
	double pi = std::numbers::pi;
	double c = std::pow((2.0 / pi), 0.75) * std::pow(2.0, _l);
	for (int i = 0; i < nExponents(); i++) {
		_coefficients[i] *= c;
		_coefficients[i] *= std::pow(_exponents[i], ((0.5 * _l) + 0.75));
	}

	double sn = std::pow(pi, 1.5) / std::pow(2.0, _l);
	double dn = 0.0; 
	// Parallize???
	for (int i = 0; i < nExponents(); i++) {
		double a = sn * _coefficients[i] * _coefficients[i];
		double b = std::pow((_exponents[i], _exponents[i]), (_l + 1.5));
		dn += (a / b);
	}
	dn = std::sqrt(dn);
	for (double i : _coefficients) {
		i /= dn;
	}
	for (int i = 0; i <= _l; i++) {
		for (int j = 0; j <= i; j++) {
			int ln = dfac(2 * (_l - i) - i);
			ln *= dfac(2 * (i - j) - 1);
			ln *= dfac(2 * j - 1);
			double norm = 1.0 / std::sqrt(ln);
			_normal.push_back(norm);
		}
	}
}

void Shell::calculateCanonicallmn() {
	_lmn.reserve((((_l + 1) * (_l + 2) / 2) * 3));
	for (int i = 0; i <= _l; i++) {
		for (int j = 0; j <= i; j++) {
			_lmn.push_back(_l - i);
			_lmn.push_back(i - j);
			_lmn.push_back(j);
		}
	}
}

//inline??
int Shell::dfac(int n) {
	int q = 1;
	while (n > 0) {
		q *= n; 
		n -= 2;
	}
	return q;
}