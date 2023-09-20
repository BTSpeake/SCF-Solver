#include "Basis.h"

#include <cmath>

Basis::Basis(int* funcMap, int lmax, int* atomMap, int nShells, double* pExp, double* cc, double* r, bool contracted, bool spherical)
	: _lmax(lmax), _contracted(contracted), _spherical(spherical)
{
	if (_contracted) {
		createContractedShells(funcMap, atomMap, nShells, pExp, cc, r);
	}
	else {
		createPrimitiveShells(funcMap, atomMap, nShells, pExp, cc, r);
	}

	if (_spherical) {
		setC2S();
	}
}

Basis::~Basis() {
	for (auto shl : _orbShells) {
		delete shl;
	}
	_orbShells.clear();
	if (_spherical) delete[] _c2sh;
}

int Basis::nShells() const {
	return _orbShells.size();
}

Shell* Basis::getShell(int i) const {
	return _orbShells[i];
}

int Basis::nSpherical() const {
	return _nSph;
}

int Basis::nCartesian() const {
	return _nCar;
}
int Basis::nFunctions() const {
	if (_spherical) { return _nSph; }
	else { return _nCar; }
}

bool Basis::isSpherical() const {
	return _spherical;
}

void Basis::createContractedShells(int* funcMap, int* atomMap, int& nShls, double* pExp, double* cc, double* r) {
	int nPrim = 0;
	for (int i = 0; i < nShls; i++) {
		nPrim += funcMap[i];
	}
	
	for (int l = 0; l <= _lmax; l++) {
		int iVal = 0;
		int lOffset = l * nPrim;
		for (int i = 0; i < nShls; i++) {
			int ri = (i * 3);
			Shell* shl = new Shell(atomMap[i], l, r[ri], r[ri + 1], r[ri + 2]);

			for (int j = 0; j < funcMap[i]; j++) {
				if (std::abs(cc[lOffset + iVal + j]) > 1e-10) {
					shl->addExponent(pExp[iVal + j], cc[lOffset + iVal + j]);
				}
			}
			if (shl->nExponents() > 0) {
				shl->calculateProperties(_nSph, _nCar);
				_orbShells.push_back(shl);
				_nSph += 2 * l + 1;
				_nCar += (l + 1) * (l + 2) / 2;
			}
			else {
				delete shl;
			}
			iVal += funcMap[i];
		}
	}
}

void Basis::createPrimitiveShells(int* funcMap, int* atomMap, int& nShls, double* pExp, double* cc, double* r) {
	int nPrim = 0;
	for (int i = 0; i < nShls; i++) {
		nPrim += funcMap[i];
	}

	for (int l = 0; l <= _lmax; l++) {
		int iVal = 0;
		int lOffset = l * nPrim;
		for (int i = 0; i < nShls; i++) {
			int ri = (i * 3);
			for (int j = 0; j < funcMap[i]; j++) {
				Shell* shl = new Shell(atomMap[i], l, r[ri], r[ri + 1], r[ri + 2]);
				if (std::abs(cc[lOffset + iVal + j]) > 1e-10) {
					shl->addExponent(pExp[iVal + j], 1.0);
				}
				if (shl->nExponents() > 0) {
					shl->calculateProperties(_nSph, _nCar);
					_orbShells.push_back(shl);
					_nSph += 2 * l + 1;
					_nCar += (l + 1) * (l + 2) / 2;
				}
				else {
					delete shl;
				}
			}
			iVal += funcMap[i];
		}
	}
}

void Basis::setC2S() {
	_lmax *= 2;
	int lmx = _lmax + 1;
	int ns = lmx + lmx + 1;
	int nc = (lmx + 1) * (lmx + 2) / 2;
	_c2shDim[1] = ns;
	_c2shDim[0] = nc;
	_c2sh = new double[lmx * nc * ns];
	for (int i = 0; i < (lmx * nc * ns); i++) {
		_c2sh[i] = 0.0;
	}
	double* fac = new double[ns];
	fac[0] = 1.0;
	for (int l = 1; l < ns; l++) {
		fac[l] = l * fac[l - 1];
	}
	for (int l = 0; l < lmx; l++) {
		for (int m = 0; m <= l; m++) {
			int lx = l - m;
			for (int n = 0; n <= m; n++) {
				int ly = m - n;
				int p = m * ((m + 1) / 2) + n;
				for (int q = -l; q <= l; q++) {
					_c2sh[(l * (nc * ns)) + (p * ns) + l + q] = C2Sco(lx, ly, n, q, fac);
				}
			}
		}
	}
	delete[] fac;
}

double Basis::C2Sco(int& lx, int& ly, int& lz, int& m, double* fac) {
	int l = lx + ly + lz;

	if (l == 0) {
		return 1.0;
	}
	else if (l == 1) {
		if (lx == 1 && m == -1) { return 1.0; }
		else if (ly == 1 && m == 0) { return 1.0; }
		else if (lz == 1 && m == 1) { return 1.0; }
		else { return 0.0; }
	}

	double sqrt2 = std::sqrt(2.0);
	int ma = std::abs(m);
	int imx = (l - ma) / 2;
	int j = (l - ma - lz) / 2;

	if ((2 * j) != (l - ma - lz)) { return 0.0; }

	double a;
	a = (fac[2 * lx] * fac[2 * ly] * fac[2 * lz] * fac[l]) / (fac[lx] * fac[ly] * fac[lz] * fac[2 * l]);
	a = std::sqrt(a);

	double b;
	b = fac[l - ma] / fac[l + ma];
	b = std::sqrt(b);
	b /= (std::pow(2.0, l) * fac[l]);

	double t{ 0.0 };

	if (lx == 2 && m == 2) {

	}

	for (int i = 0; i <= imx; i++) {
		if (j >= 0 && j <= i) {

			double c = (fac[l] / (fac[i] * fac[l - i]));
			c *= fac[2 * l - 2 * i];
			c *= (std::pow(-1.0, i) / fac[l - ma - 2 * i]);
			c *= (fac[i] / (fac[j] * fac[i - j]));

			for (int k = 0; k <= j; k++) {
				double e{ 0.0 };
				int lxm2k = lx - 2 * k;
				if (lxm2k >= 0 && lxm2k <= ma) {

					double d = (fac[j] / (fac[k] * fac[j - k]));
					d *= (fac[ma] / (fac[lx - 2 * k] * fac[ma + 2 * k - lx]));

					if (m == 0 && (lx % 2) == 0) {
						int p{ k - lx / 2 };
						e = std::pow(-1.0, p);
					}
					else if (m > 0 && (std::abs(ma - lx) % 2) == 0) {
						int p{ (2 * k + ma - lx) / 2 };
						e = sqrt2 * std::pow(-1.0, p);
					}
					else if (m < 0 && (std::abs(ma - lx) % 2) == 1) {
						int p{ (2 * k + ma - lx) / 2 };
						e = sqrt2 * std::pow(-1.0, p);
					}
					t += c * d * e;
				}
			}
		}
	}

	double co = a * b * t;
	return co;
}
