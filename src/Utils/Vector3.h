#pragma once 

class Vector3 {
public:
	Vector3();
	Vector3(double x, double y, double z);
	~Vector3();
	double operator[](const int& i) const;
	Vector3 operator-(const Vector3& rhs) const;
	Vector3 operator+(const Vector3& rhs) const;
	Vector3& operator=(const Vector3 rhs);
	Vector3& operator+=(const Vector3& rhs);
	Vector3& operator-=(const Vector3& rhs);
	void scalarMultiply(const double& s);
	void scalarDivide(const double& s);

private:
	double _x, _y, _z;
};