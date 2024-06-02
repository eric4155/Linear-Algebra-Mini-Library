#include "matrix.hpp"
#include "vector.hpp"
#include <iostream>
#include <cassert>
#include <vector>

void TestMatrix() {
  // Test 1: Create a 3x3 matrix and initialize it
  Matrix mat1(3, 3);
  mat1.SetElement(1.0, 0, 0);
  mat1.SetElement(2.0, 0, 1);
  mat1.SetElement(3.0, 0, 2);
  mat1.SetElement(4.0, 1, 0);
  mat1.SetElement(5.0, 1, 1);
  mat1.SetElement(6.0, 1, 2);
  mat1.SetElement(7.0, 2, 0);
  mat1.SetElement(8.0, 2, 1);
  mat1.SetElement(9.0, 2, 2);
  std::cout << "Test 1: Initial 3x3 matrix\n" << mat1;

  // Test 2: Copy constructor
  Matrix mat2(mat1);
  std::cout << "Test 2: Copy of the initial matrix\n" << mat2;

  // Test 3: Copy assignment
  Matrix mat3;
  mat3 = mat1;
  std::cout << "Test 3: Copy assignment of the initial matrix\n" << mat3;

  // Test 4: Row Swap
  mat1.RowSwap(0, 2);
  std::cout << "Test 4: Matrix after swapping row 0 and row 2\n" << mat1;

  // Test 5: Row Replacement
  mat1.RowReplacement(1, 2, 0.5);
  std::cout << "Test 5: Matrix after replacing row 2 with row 2 + 0.5*row 1\n" << mat1;

  // Test 6: Scale Row
  mat1.ScaleRow(1, 2.0);
  std::cout << "Test 6: Matrix after scaling row 1 by 2.0\n" << mat1;

  // Test 7: Set to Identity
  mat1.SetToIdentity();
  std::cout << "Test 7: Matrix set to identity\n" << mat1;

  // Test 8: Matrix addition
  Matrix mat4(3, 3);
  mat4.SetElement(1.0, 0, 0);
  mat4.SetElement(2.0, 0, 1);
  mat4.SetElement(3.0, 0, 2);
  mat4.SetElement(4.0, 1, 0);
  mat4.SetElement(5.0, 1, 1);
  mat4.SetElement(6.0, 1, 2);
  mat4.SetElement(7.0, 2, 0);
  mat4.SetElement(8.0, 2, 1);
  mat4.SetElement(9.0, 2, 2);
  Matrix mat5 = mat4 + mat4;
  std::cout << "Test 8: Matrix after addition\n" << mat5;

  // Test 9: Matrix subtraction
  Matrix mat6 = mat4 - mat4;
  std::cout << "Test 9: Matrix after subtraction\n" << mat6;

  // Test 10: Scalar multiplication
  Matrix mat7 = 2.0 * mat4;
  std::cout << "Test 10: Matrix after scalar multiplication\n" << mat7;

  // Test 11: Matrix multiplication
  Matrix mat8(3, 2);
  mat8.SetElement(1.0, 0, 0);
  mat8.SetElement(2.0, 0, 1);
  mat8.SetElement(3.0, 1, 0);
  mat8.SetElement(4.0, 1, 1);
  mat8.SetElement(5.0, 2, 0);
  mat8.SetElement(6.0, 2, 1);
  Matrix mat9 = mat4 * mat8;
  std::cout << "Test 11: Matrix after matrix multiplication\n" << mat9;

  // Test 12: Get Row Vector
  double* row_vec = mat4.GetRowVec(1);
  std::cout << "Test 12: Row vector of row 1\n";
  for (int i = 0; i < mat4.GetCols(); ++i) {
    std::cout << row_vec[i] << " ";
  }
  std::cout << "\n";
  delete[] row_vec;

  // Test 13: Get Column Vector
  double* col_vec = mat4.GetColumnVec(1);
  std::cout << "Test 13: Column vector of column 1\n";
  for (int i = 0; i < mat4.GetRows(); ++i) {
    std::cout << col_vec[i] << " ";
  }
  std::cout << "\n";
  delete[] col_vec;

  // Test 14: Matrix Join
  Matrix mat10(3, 1);
  mat10.SetElement(10.0, 0, 0);
  mat10.SetElement(20.0, 1, 0);
  mat10.SetElement(30.0, 2, 0);
  Matrix mat11 = mat4.Join(mat10);
  std::cout << "Test 14: Matrix after joining with another matrix\n" << mat11;

  // Test 15: Equality operator
  bool are_equal = (mat4 == mat4);
  std::cout << "Test 15: Check if a matrix is equal to itself: " << (are_equal ? "True" : "False") << "\n";
}

int main() {
  // TestMatrix();
  /*
  std::vector<double> data = {1,2,3,4};
  Matrix test(2,2,data);
  Matrix test_inverse = test.Invert();
  Matrix identity(2,2);
  identity.SetToIdentity();
  assert(test * test_inverse == identity);

  std::vector<double> data2 = {1,2,-1,-2,0,1,1,-1,0};
  Matrix test2(3,3,data2);
  Matrix test2_inverse = test2.Invert();
  Matrix identity2(3,3);
  identity2.SetToIdentity();
  assert(test2 * test2_inverse == identity2);

  std::vector<double> data3 = {1,1,0,1,0,3,8,1,1,2,0,0,2,2,3,1};
  Matrix test3(4,4, data3);
  Matrix test3_inverse = test3.Invert();
  Matrix identity3(4,4);
  identity3.SetToIdentity();
  assert(test3 * test3_inverse == identity3);
  std::cout << "Passed tests"; */
  
  /*std::vector<double> data1 = {
    2, -1, 0, 0, 0,
    -1, 2, -1, 0, 0,
    0, -1, 2, -1, 0,
    0, 0, -1, 2, -1,
    0, 0, 0, -1, 2
  };
  Matrix test1(5, 5, data1);
  Matrix test1_inverse = test1.Invert();
  Matrix identity1(5, 5);
  identity1.SetToIdentity();
  assert(test1 * test1_inverse == identity1);
  std::vector<double> data2 = {
    1, 0, 0, 0, 0, 0,
    0, 2, 0, 0, 0, 0,
    0, 0, 3, 0, 0, 0,
    0, 0, 0, 4, 0, 0,
    0, 0, 0, 0, 5, 0,
    0, 0, 0, 0, 0, 6
  };
  Matrix test2(6, 6, data2);
  Matrix test2_inverse = test2.Invert();
  Matrix identity2(6, 6);
  identity2.SetToIdentity();
  assert(test2 * test2_inverse == identity2);
  std::vector<double> data3 = {
    1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 1
  };
  Matrix test3(7, 7, data3);
  Matrix test3_inverse = test3.Invert();
  Matrix identity3(7, 7);
  identity3.SetToIdentity();
  assert(test3 * test3_inverse == identity3);

  std::vector<double> data4 = {
    2, 1, 0, 0, 0, 0, 0, 0,
    1, 2, 1, 0, 0, 0, 0, 0,
    0, 1, 2, 1, 0, 0, 0, 0,
    0, 0, 1, 2, 1, 0, 0, 0,
    0, 0, 0, 1, 2, 1, 0, 0,
    0, 0, 0, 0, 1, 2, 1, 0,
    0, 0, 0, 0, 0, 1, 2, 1,
    0, 0, 0, 0, 0, 0, 1, 2
  };
  Matrix test4(8, 8, data4);
  Matrix test4_inverse = test4.Invert();
  Matrix identity4(8, 8);
  identity4.SetToIdentity();
  assert(test4 * test4_inverse == identity4);

  std::vector<double> data5 = {
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };
  Matrix test5(10, 10, data5);
  Matrix test5_inverse = test5.Invert();
  Matrix identity5(10, 10);
  identity5.SetToIdentity();
  assert(test5 * test5_inverse == identity5); */
  /* std::vector<double> data = {1,2,3,4};
  Matrix test(2,2,data);
  std::cout << "Detemrinant 1: " << test.Determinant() << std::endl;

  std::vector<double> data1 = {1,2,-1,-2,0,1,1,-1,0};
  Matrix test2(3,3, data1);
  std::cout << "Determinant 2: " << test2.Determinant() << std::endl;*/
  // std::vector<double> data10 = {100,0,0,1,0,0,0,200,1};
  // Matrix test10(3,3, data10);
  // std::vector<double> data = {100,1,201};
  // Vector test1(data);
  // std::cout << test10.Solve(test1);
  //std::cout << test1.Magnitude() << std::endl;
  //std::vector<double> data1 = {1,1,1};
  //Vector test2(data1);
  //std::cout << test1.Dot(test2) << std::endl;
  std::vector<double> test100 = {123,90,100,40,93,21,90,23,101};
  Matrix test101(3,3, test100);
  Matrix inverse = test101.Invert();
  Matrix identity(3,3);
  identity.SetToIdentity();
  assert(test101 * inverse == identity);
  std::cout << "Correct" << std::endl;

  std::vector<double> data101 = {101, 200, 300, 431};
  Matrix test200(2,2,data101);
  std::cout << test200 * test200.Invert() << std::endl;

  return 0;
}
