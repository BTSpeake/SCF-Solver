#include "ShellPair.h"

#include <cmath>
#include <numbers>

#include <iostream>

ShellPair::ShellPair(Shell* shli, Shell* shlj, bool spherical)
	: _shli(shli), _shlj(shlj)
{

	int la = shli->getAngularMomentum();
	_l = la + shlj->getAngularMomentum();
	np = (_l + 1) * (_l + 2) * (_l + 3) / 6;
	np -= ((la) * (la + 1) * (la + 2) / 6);

	_k = shli->nExponents() * shlj->nExponents();

	if (_k == 1) {
		calculatePrimitive();
	}
	else {
		calculateContracted();
	}

	std::cout << "New Shell Pair" << std::endl;
	std::cout << _u[0] << std::endl;

}

ShellPair::~ShellPair() {
	delete[] _r; 
	delete[] _u; 
	delete[] _s;

	delete[] _sa;
	delete[] _sb;
}


const int ShellPair::getli() const {
	return _shli->getAngularMomentum();
}
const int ShellPair::getlj() const {
	return _shlj->getAngularMomentum();
}
int ShellPair::getk() const {
	return _k;
}
std::complex<double> ShellPair::getrp(int i) const {
	return _r[i];
}
double ShellPair::getri(int i) const {
	return _shli->getPosition()[i];
}
double ShellPair::getrj(int i) const {
	return _shlj->getPosition()[i];
}
std::complex<double> ShellPair::getu(int i) const {
	return _u[i];
}
std::complex<double> ShellPair::gets(int i) const {
	return _s[i];
}
int ShellPair::nCartesiani() const {
	return _shli->nCartesian();
}
int ShellPair::nCartesianj() const {
	return _shlj->nCartesian();
}
int ShellPair::getlmni(int k, int i) const {
	return _shli->getlmn(k, i);
}
int ShellPair::getlmnj(int k, int i) const {
	return _shlj->getlmn(k, i);
}
const double& ShellPair::getNormali(int i) const {
	return _shli->getNormal(i);
}
const double& ShellPair::getNomralj(int i) const {
	return _shlj->getNormal(i);
}
int ShellPair::getCartesianStarti() const {
	return _shli->getCartesianStart();
}
int ShellPair::getCartesianStartj() const {
	return _shlj->getCartesianStart();
}


void ShellPair::calculatePrimitive() {

	_r = new std::complex<double>[3];
	_u = new std::complex<double>[1];
	_s = new std::complex<double>[1];
	
	_k1[0] = _shlj->_kpw[0] - _shli->_kpw[0];
	_k1[1] = _shlj->_kpw[1] - _shli->_kpw[1];
	_k1[2] = _shlj->_kpw[2] - _shli->_kpw[2];

	std::complex<double> k1s = (_k1[0] * _k1[0]) + (_k1[1] * _k1[1]) + (_k1[2] * _k1[2]);
	if (k1s != 0.0) {
		_occ = false;
	}

	double rab[3];
	rab[0] = _shlj->getPosition()[0] - _shli->getPosition()[0];
	rab[1] = _shlj->getPosition()[1] - _shli->getPosition()[1];
	rab[2] = _shlj->getPosition()[2] - _shli->getPosition()[2];
	double rsq = (rab[0] * rab[0]) + (rab[1] * rab[1]) + (rab[2] * rab[2]);
	if (rsq != 0.0) {
		_occ = false;
	}

	_sa = new std::complex<double>[1];
	_sa[0] = 2.0 * _shli->getExponent(0);
	_sb = new std::complex<double>[1];
	_sb[0] = 2.0 * _shlj->getExponent(0);
	_s[0] = 1.0 / (_sa[0] + _sb[0]);

	_r[0] = ((_sa[0] * _shli->getPosition()[0]) + (_sb[0] * _shlj->getPosition()[0])) * _s[0];
	_r[1] = ((_sa[0] * _shli->getPosition()[1]) + (_sb[0] * _shlj->getPosition()[1])) * _s[0];
	_r[2] = ((_sa[0] * _shli->getPosition()[2]) + (_sb[0] * _shlj->getPosition()[2])) * _s[0];
	std::complex<double> k1p{ (_k1[0] * _r[0]) + (_k1[1] * _r[1]) + (_k1[2] * _r[2]) };

	std::complex<double> j1(0.0, 1.0);
	double pi = std::numbers::pi;
	_u[0] = std::pow(2.0 * pi * _s[0], 1.5);
	_u[0] *= std::exp(-0.5 * _sa[0] * _sb[0] * _s[0] * rsq);
	_u[0] *= std::exp(-0.5 * _s[0] * k1s - j1 * k1p);
	_u[0] *= _shli->getCoefficient(0) * _shlj->getCoefficient(0);

	_r[0] -= j1 * _s[0] * _k1[0];
	_r[1] -= j1 * _s[0] * _k1[1];
	_r[2] -= j1 * _s[0] * _k1[2];

}

void ShellPair::calculateContracted() {

	_k1[0] = _shlj->_kpw[0] - _shli->_kpw[0];
	_k1[1] = _shlj->_kpw[1] - _shli->_kpw[1];
	_k1[2] = _shlj->_kpw[2] - _shli->_kpw[2];

	std::complex<double> k1s = (_k1[0] * _k1[0]) + (_k1[1] * _k1[1]) + (_k1[2] * _k1[2]);
	if (k1s != 0.0) {
		_occ = false;
	}

	double rab[3];
	rab[0] = _shlj->getPosition()[0] - _shli->getPosition()[0];
	rab[1] = _shlj->getPosition()[1] - _shli->getPosition()[1];
	rab[2] = _shlj->getPosition()[2] - _shli->getPosition()[2];
	double rsq = (rab[0] * rab[0]) + (rab[1] * rab[1]) + (rab[2] * rab[2]);
	if (rsq != 0.0) {
		_occ = false;
	}

	_r = new std::complex<double>[3 * _k];
	_u = new std::complex<double>[_k];
	_sa = new std::complex<double>[_k];
	_sb = new std::complex<double>[_k];
	_s = new std::complex<double>[_k];

	int ki = 0;
	std::complex<double> j1(0.0, 1.0);
	double pi = std::numbers::pi;

	for (int i = 0; i < _shli->nExponents(); i++) {
		if (ki == _k) {
			//std::cout << "ERROR:: Contracted Shell Pair numbering error k = " << _k << "  ki = " << ki << std::endl;
			break;
		}
		_sa[ki] = 2.0 * _shli->getExponent(i);
		for (int j = 0; j < _shlj->nExponents(); j++) {
			if (ki == _k) {
				//std::cout << "ERROR:: Contracted Shell Pair numbering error k = " << _k << "  ki = " << ki << std::endl;
				break;
			}
			_sb[ki] = 2.0 * _shli->getExponent(i);
			_s[ki] = 1.0 / (_sa[ki] + _sb[ki]);

			_r[(3 * ki) + 0] = ((_sa[ki] * _shli->getPosition()[0]) + (_sb[ki] * _shlj->getPosition()[0])) * _s[ki];
			_r[(3 * ki) + 1] = ((_sa[ki] * _shli->getPosition()[1]) + (_sb[ki] * _shlj->getPosition()[1])) * _s[ki];
			_r[(3 * ki) + 2] = ((_sa[ki] * _shli->getPosition()[2]) + (_sb[ki] * _shlj->getPosition()[2])) * _s[ki];
			std::complex<double> k1p{ (_k1[0] * _r[0]) + (_k1[1] * _r[1]) + (_k1[2] * _r[2]) };

			_u[0] = std::pow(2.0 * pi * _s[ki], 1.5);
			_u[0] *= std::exp(-0.5 * _sa[ki] * _sb[ki] * _s[ki] * rsq);
			_u[0] *= std::exp(-0.5 * _s[ki] * k1s - j1 * k1p);
			_u[0] *= _shli->getCoefficient(i) * _shlj->getCoefficient(j);

			_r[(3 * ki) + 0] -= j1 * _s[0] * _k1[0];
			_r[(3 * ki) + 1] -= j1 * _s[0] * _k1[1];
			_r[(3 * ki) + 2] -= j1 * _s[0] * _k1[2];

			ki++;
		}
	}

}

void ShellPair::calculateHrrMatrix() {

}