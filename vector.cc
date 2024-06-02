#include "vector.hpp"


Vector::Vector(const std::vector<double>& data) {
    nDimensions_ = static_cast<int>(data.size());
    data_ = data;
}

/***************
 * OPERATORS
*/
Vector Vector::operator+(const Vector& rhs) const {
    if (nDimensions_ != rhs.nDimensions_) { throw std::invalid_argument("Different dimensions"); }
    Vector resultant(data_);
    for (int i = 0; i < nDimensions_; ++i) {
        resultant.data_[i] = data_[i] + rhs.data_[i];
    }
    return resultant;
}

Vector Vector::operator-(const Vector& rhs) const {
    if (nDimensions_ != rhs.nDimensions_) { throw std::invalid_argument("Different dimensions"); }
    Vector resultant(data_);
    for (int i = 0; i < nDimensions_; ++i) {
        resultant.data_[i] = data_[i] - rhs.data_[i];
    }
    return resultant;
}

Vector Vector::operator*(const double& factor) const {
    Vector result(data_);
    for (int i = 0; i < nDimensions_; ++i) {
        result.data_[i] *= factor;
    }
    return result;
}

Vector operator*(const double& factor, const Vector& vec) {
    Vector result(vec.data_);
    for (int i = 0; i < vec.nDimensions_; ++i) {
        result.data_[i] *= factor;
    }
    return result;
}

int Vector::GetDimensions() const {
    return nDimensions_;
}

double Vector::Dot(const Vector& rhs) const {
    if (nDimensions_ != rhs.nDimensions_) { throw std::invalid_argument("Different dimensions"); }
    double result = 0.0;
    for (int i = 0; i < nDimensions_; ++i) {
        result += (data_[i] * rhs.data_[i]);
    }
    return result;
}

double Vector::Magnitude() const {
    if (nDimensions_ == 0) { throw std::invalid_argument("Vector DNE"); }
    if (nDimensions_ == 1) { return data_[0]; }
    double result = 0.0;
    for (int i = 0; i < nDimensions_; ++i) {
        result += data_[i] * data_[i];
    }
    return sqrt(result);
}

/* *********************
NORMALIZATION 
*********************/
void Vector::NormalizeInPlace() {
    double magnitude = this->Magnitude();
    for (int i = 0; i < nDimensions_; ++i) {
        data_[i] /= magnitude;
    }
}

Vector Vector::GetNormalizedCopy() const {
    Vector normalized(data_);
    double norm = this->Magnitude();
    return normalized * (static_cast<double>(1) / norm);
}

std::vector<double> Vector::GetData() const {
    return data_;
}

std::ostream& operator<<(std::ostream& output, const Vector& vector) {
    for (int i = 0; i < vector.nDimensions_; ++i) {
        output << vector.data_[i] << " ";
    }
    output << std::endl;
    return output;
}

