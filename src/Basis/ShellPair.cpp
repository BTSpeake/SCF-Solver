#include "ShellPair.h"

#include <cmath>
#include <numbers>
#include <iostream>

ShellPair::ShellPair(Shell* shli, Shell* shlj, bool spherical)
	: _shli(shli), _shlj(shlj)
{

	int la = shli->getAngularMomentum();
	l = la + shlj->getAngularMomentum();
	np = (l + 1) * (l + 2) * (l + 3) / 6;
	np -= ((la) * (la + 1) * (la + 2) / 6);

	if (shli->nExponents() * shlj->nExponents() == 1) {
		k = 1;
		calculatePrimitive();
	}
	else {
		calculateContracted();
	}

	std::cout << sa[0] << ", " << sb[0] << std::endl;
}

ShellPair::~ShellPair() {
	delete[] r; 
	delete[] u; 
	delete[] s;

	delete[] sa;
	delete[] sb;
}

void ShellPair::calculatePrimitive() {

	r = new std::complex<double>[3];
	u = new std::complex<double>[1];
	s = new std::complex<double>[1];
	
	k1[0] = _shlj->_kpw[0] - _shli->_kpw[0];
	k1[1] = _shlj->_kpw[1] - _shli->_kpw[1];
	k1[2] = _shlj->_kpw[2] - _shli->_kpw[2];

	std::complex<double> k1s = (k1[0] * k1[0]) + (k1[1] * k1[1]) + (k1[2] * k1[2]);
	if (k1s != 0.0) {
		occ = false;
	}

	double rab[3];
	rab[0] = _shlj->getPosition()[0] - _shli->getPosition()[0];
	rab[1] = _shlj->getPosition()[1] - _shli->getPosition()[1];
	rab[2] = _shlj->getPosition()[2] - _shli->getPosition()[2];
	double rsq = (rab[0] * rab[0]) + (rab[1] * rab[1]) + (rab[2] * rab[2]);
	if (rsq != 0.0) {
		occ = false;
	}

	sa = new std::complex<double>[1];
	sa[0] = 2.0 * _shli->getExponent(0);
	sb = new std::complex<double>[1];
	sb[0] = 2.0 * _shlj->getExponent(0);
	s[0] = 1.0 / (sa[0] + sb[0]);

	r[0] = ((sa[0] * _shli->getPosition()[0]) + (sb[0] * _shlj->getPosition()[0])) * s[0];
	r[1] = ((sa[0] * _shli->getPosition()[1]) + (sb[0] * _shlj->getPosition()[1])) * s[0];
	r[2] = ((sa[0] * _shli->getPosition()[2]) + (sb[0] * _shlj->getPosition()[2])) * s[0];
	std::complex<double> k1p{ (k1[0] * r[0]) + (k1[1] * r[1]) + (k1[2] * r[2]) };

	std::complex<double> j1(0.0, 1.0);
	double pi = std::numbers::pi;
	u[0] = std::pow(2.0 * pi * s[0], 1.5);
	u[0] *= std::exp(-0.5 * sa[0] * sb[0] * s[0] * rsq);
	u[0] *= std::exp(-0.5 * s[0] * k1s - j1 * k1p);
	u[0] *= _shli->getCoefficient(0) * _shlj->getCoefficient(0);

	r[0] -= j1 * s[0] * k1[0];
	r[1] -= j1 * s[0] * k1[1];
	r[2] -= j1 * s[0] * k1[2];

}

void ShellPair::calculateContracted() {

}