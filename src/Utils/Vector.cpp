#include "Vector3.h"

Vector3::Vector3()
	: _x(0.0), _y(0.0), _z(0.0) {}

Vector3::Vector3(double x, double y, double z)
	: _x(x), _y(y), _z(z) {}

Vector3::~Vector3() {}

double Vector3::operator[](const int& i) const {
	switch (i) {
	case 0:
		return _x;
	case 1:
		return _y;
	case 2:
		return _z;
	}
	//raise an exception if not handled properly
}

Vector3 Vector3::operator-(const Vector3& rhs) const {
	double x = _x - rhs[0];
	double y = _y - rhs[1];
	double z = _z - rhs[2];
	return Vector3(x, y, z);
}

Vector3 Vector3::operator+(const Vector3& rhs) const {
	double x = _x + rhs[0];
	double y = _y + rhs[1];
	double z = _z + rhs[2];
	return Vector3(x, y, z);
}

Vector3& Vector3::operator=(const Vector3 rhs) {
	_x = rhs[0];
	_y = rhs[1];
	_z = rhs[2];
	return *this;
}
Vector3& Vector3::operator+=(const Vector3& rhs) {
	Vector3 tmp = (*this) + rhs;
	(*this) = tmp;
	return *this;
}
Vector3& Vector3::operator-=(const Vector3& rhs) {
	Vector3 tmp = (*this) - rhs;
	(*this) = tmp;
	return *this;
}

void Vector3::scalarMultiply(const double& s) {
	_x *= s;
	_y *= s;
	_z *= s;
}

void Vector3::scalarDivide(const double& s) {
	_x /= s;
	_y /= s;
	_z /= s;
}