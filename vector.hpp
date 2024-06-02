#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>


class Vector {

public:
    Vector() = default;

    Vector(const std::vector<double>& data);
    // OPERATORS
    Vector operator+(const Vector& rhs) const;
    Vector operator-(const Vector& rhs) const;
    Vector operator*(const double& factor) const;
    friend Vector operator*(const double& factor, const Vector& vec);

    double Dot(const Vector& rhs) const;

    // GETTERS 
    int GetDimensions() const;
    std::vector<double> GetData() const;

    double Magnitude() const;
    // NORMALIZATION METHODS
    void NormalizeInPlace();
    Vector GetNormalizedCopy() const;

    friend std::ostream& operator<<(std::ostream& output, const Vector& vector);

private:
    int nDimensions_ = 0;
    std::vector<double> data_;
};


#endif