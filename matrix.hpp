#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "vector.hpp"

constexpr double kTolerance = 1e-9;
class Matrix {
public:
  // Declare Constructors
  Matrix() = default;
  Matrix(int rows, int cols);
  Matrix(int rows, int cols, std::vector<double>& data);
  Matrix(const Matrix& matrix);

  // Destructor
  ~Matrix();

  // Copy assignment
  Matrix& operator=(const Matrix& rhs);


  // Overload comparison operator to check if 2 matrices are equivalent; method thatr tests if 2 matrices are equal given a specific tolerance
  friend bool operator==(const Matrix& lhs, const Matrix& rhs);
  bool Compare(const Matrix& matrix, double tolerance);

  // Matrix operations defined as overloaded operators
  friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
  friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
  friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
  friend Matrix operator*(const double& scalar, const Matrix& rhs);
  friend Matrix operator*(const Matrix& lhs, const double& scalar);
  friend Vector operator*(const Matrix& matrix, const Vector& vector);


  // Configuration
  void SetToIdentity();
  void SetToUpperTriangular();
  void ReducedEchelon();


  // Helper methods to get column/row vectors at a specified index (0 indexed)
  double* GetColumnVec(int col) const;
  double* GetRowVec(int row) const;
  // Overloaded stream insertion operator to display matrices 
  friend std::ostream& operator<<(std::ostream& output, const Matrix& matrix);

  // Getters and Setters
  void SetElement(double data, int row, int col);
  int GetCols() const;
  int GetRows() const;
  double** GetMatrix() const;

  bool IsSquare();
  bool CloseEnough(double num1, double num2) const;
  bool IsEchelon() const;

  // Elementary row ops.
  void RowSwap(int i, int j);
  void RowReplacement(int i, int j, double factor);
  void ScaleRow(int i, double factor);
  // RETURNS A MATRIX THAT REPRESENTS THE 2 MATRICES COMBINED
  Matrix Join(const Matrix& matrix2);
  
  // COMPUTE MATRIX INVERSE / TRANSPOSE
  Matrix Invert();
  Matrix Transpose() const;

  // COMPUTE DETERMINANT
  double Determinant() const;

  // SOLVE SYSTEM OF EQUATIONS
  Vector Solve(const Vector& vector);

private:
  double Dot(double* vec1, double* vec2, int size);
  void Deallocate();
  int nRows_ = 0;
  int nCols_ = 0;
  double** matrix_ptr_ = nullptr;
};

#endif
/***********************************************
 * TODO: FINISH FUNCTION TO SOLVE SYSTEM OF EQUATIONS  AND ROWECHELON FORM
*/