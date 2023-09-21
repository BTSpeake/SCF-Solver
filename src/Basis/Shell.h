#pragma once 

#include "../Utils/Vector3.h"
#include <vector>
#include <complex>

class Shell {
public:
	Shell(int i, int l, double px, double py, double pz);
	~Shell();

	//void setPosition(const double x, const double y, const double z);
	void addExponent(const double exp, const double cc);

	int getAngularMomentum() const;

	const std::vector<double>& getExponents() const;
	const double& getExponent(int i) const;
	int nExponents() const; 
	const std::vector<double>& getCoefficients() const;
	const double& getCoefficient(int i) const;
	int nCoefficients() const;

	int getSphericalStart() const;
	int nSpherical() const;
	int getCartesianStart() const;
	int nCartesian() const;

	const Vector3& getPosition() const;

	int getlmn(int k, int i) const;
	const double& getNormal(int i) const;

	void calculateProperties(int si, int ci);

	std::complex<double> _kpw[3]{ 0.0, 0.0, 0.0 };

private:
	void calculateGaussianNormal();
	void calculateCanonicallmn();
	int dfac(int n);


	int _l;
	int _atom;  
	int _si, _ci;
	std::vector<double> _exponents;
	std::vector<double> _coefficients;
	std::vector<double> _normal;
	std::vector<int> _lmn;
	Vector3 _position;

};
