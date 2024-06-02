#include "matrix.hpp"

// ************************* 
// HELPER FUNCTION TO DEALLOCATE MEMORY
void Matrix::Deallocate() {
  for (int row = 0; row < nRows_; ++row) {
    delete[] matrix_ptr_[row];
  }
  delete[] matrix_ptr_;
  matrix_ptr_ = nullptr;
}

/* ********************************************* 
CONSTRUCTORS
*************************************************/
Matrix::Matrix(int rows, int cols): nRows_(rows), nCols_(cols) {
  matrix_ptr_ = new double*[nRows_];
  for (int i = 0; i < nRows_; ++i) {
    auto* col = new double[nCols_];
    for (int j = 0; j < nCols_; ++j) {
      col[j] = 0;
    }
    matrix_ptr_[i] = col;
  }
}

// Constructor that inits a matrix w the specified data
Matrix::Matrix(int rows, int cols, std::vector<double>& data):
    nRows_(rows), nCols_(cols) {
  size_t data_idx = 0;
  matrix_ptr_ = new double*[nRows_];
  for (int row = 0; row < nRows_; ++row) {
    auto* new_col = new double[nCols_];
    for (int col = 0; col < nCols_; ++col) {
      new_col[col] = data[data_idx];
      data_idx++;
    }
    matrix_ptr_[row] = new_col;
  }
}

// Copy contructor implementing a deep copy policy on the matrix
Matrix::Matrix(const Matrix& matrix):
    nRows_(matrix.nRows_), nCols_(matrix.nCols_) {
  matrix_ptr_ = new double*[nRows_];
  for (int row = 0; row < nRows_; ++row) {
    auto* new_col = new double[nCols_];
    for (int col = 0; col < nCols_; ++col) {
      new_col[col] = matrix.matrix_ptr_[row][col];
    }
    matrix_ptr_[row] = new_col;
  }
}

/*************************
 * DESTRUCTOR
 * **************************
*/
Matrix::~Matrix() {
  nRows_ = 0;
  nCols_ = 0;
  Deallocate();
}


/************************************
 * OVERLOADED OPERATORS 
 * ***********************************
*/

// COMPARISON FUNCTION THAT TELLS WHETHER OR NOT 2 MATRICES ARE THE SAME GIVEN A TOLERANCE
bool Matrix::Compare(const Matrix& matrix, double tolerance) {
  if (nCols_ != matrix.nCols_ || nRows_ != matrix.nRows_) { return false; }
  for (int row = 0; row < nRows_; ++row) {
    for (int col = 0; col < nCols_; ++col) {
      if (fabs(matrix_ptr_[row][col] - matrix.matrix_ptr_[row][col]) > tolerance) { return false; }
    }
  }
  return true;
}

// Copy assignment operator
Matrix& Matrix::operator=(const Matrix& rhs) {
  if (this == &rhs) {
    return *this;
  }
  this->Deallocate();
  nRows_ = rhs.nRows_;
  nCols_ = rhs.nCols_;
  matrix_ptr_ = new double*[nRows_];
  for (int row = 0; row < nRows_; ++row) {
    matrix_ptr_[row] = new double[nCols_];
    for (int col = 0; col < nCols_; ++col) {
      matrix_ptr_[row][col] = rhs.matrix_ptr_[row][col];
    }
  }
  return *this;
}

// Addition operator
Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.nRows_ != rhs.nRows_ || lhs.nCols_ != rhs.nCols_) {
    throw std::invalid_argument("Matrices must have identical dimensions");
  }
  Matrix to_return(lhs.nRows_, lhs.nCols_);
  for (int row = 0; row < lhs.nRows_; ++row) {
    for (int col = 0; col < lhs.nCols_; ++col) {
      to_return.matrix_ptr_[row][col] = lhs.matrix_ptr_[row][col] + rhs.matrix_ptr_[row][col];
    }
  }
  return to_return;
}

// Subtracting operator
Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.nRows_ != rhs.nRows_ || lhs.nCols_ != rhs.nCols_) {
    throw std::invalid_argument("Matrices must have identical dimensions");
  }
  Matrix to_return(lhs.nRows_, lhs.nCols_);
  for (int row = 0; row < lhs.nRows_; ++row) {
    for (int col = 0; col < lhs.nCols_; ++col) {
      to_return.matrix_ptr_[row][col] = lhs.matrix_ptr_[row][col] - rhs.matrix_ptr_[row][col];
    }
  }
  return to_return;
}

// Equality operator
bool operator==(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.nRows_ != rhs.nRows_ || lhs.nCols_ != rhs.nCols_) {
    return false;
  }
  for (int row = 0; row < lhs.nRows_; ++row) {
    for (int col = 0; col < lhs.nCols_; ++col) {
      if (!(lhs.CloseEnough(lhs.matrix_ptr_[row][col],rhs.matrix_ptr_[row][col]))) {
        return false;
      }
    }
  }
  return true;
}

// Matrix scalar mult.
Matrix operator*(const Matrix& lhs, const double& scalar) {
  if (lhs.matrix_ptr_ == nullptr) {
    throw std::invalid_argument("Matrix does not exist");
  }
  Matrix to_return(lhs.nRows_, lhs.nCols_);
  for (int row = 0; row < lhs.nRows_; ++row) {
    for (int col = 0; col < lhs.nCols_; ++col) {
      to_return.matrix_ptr_[row][col] = (lhs.matrix_ptr_[row][col] * scalar);
    }
  }
  return to_return;
}
// SCALAR MATRIX MULT
Matrix operator*(const double& scalar, const Matrix& rhs) {
    if (rhs.matrix_ptr_ == nullptr) {
        throw std::invalid_argument("Matrix does not exist");
    }
    Matrix to_return(rhs.nRows_, rhs.nCols_);
    for (int row = 0; row < rhs.nRows_; ++row) {
        for (int col = 0; col < rhs.nCols_; ++col) {
            to_return.matrix_ptr_[row][col] = rhs.matrix_ptr_[row][col] * scalar;
        }
    }
    return to_return;
}

// Matrix matrix multiplication
Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
  // Check that dimensions are satisfied
  if (lhs.nCols_ != rhs.nRows_) {
    throw std::invalid_argument("LHS columns must match RHS rows");
  }
  Matrix new_matrix(lhs.nRows_, rhs.nCols_);
  // Loop through the new dimensions of the matrix
  for (int row = 0; row < lhs.nRows_; ++row) {
    for (int col = 0; col < rhs.nCols_; ++col) {
      double* col_vec = rhs.GetColumnVec(col);
      for (int lhs_row = 0; lhs_row < lhs.nRows_; ++lhs_row) {
        // Compute dot product of each corresponding row and column vector
        double* row_vec = lhs.GetRowVec(lhs_row);
        new_matrix.matrix_ptr_[lhs_row][col] = new_matrix.Dot(row_vec, col_vec, lhs.nCols_);
        delete[] row_vec;
      }
      delete[] col_vec;
    }
  }
  return new_matrix;
}

// Matrix Vector mult.
Vector operator*(const Matrix& matrix, const Vector& vector) {
  if (matrix.nCols_ != vector.GetDimensions()) { throw std::invalid_argument("Invalid dimensions"); }
  std::vector<double> result_vec;
  for (int row = 0; row < matrix.GetRows(); ++row) {
    double cumulative_sum = 0.0;
    for (int col = 0; col < matrix.GetCols(); ++col) {
      cumulative_sum += (vector.GetData()[col] * matrix.GetMatrix()[row][col]);
    }
    result_vec.push_back(cumulative_sum);
  }
  Vector vec(result_vec);
  return vec;
}

/* ********************************
 HELPER FUNCTIONS TO GET ROW/COLUMN VECTORS AND COMPUTE DOT PRODUCT
 ***********************************
*/
// Get Column vector from a matrix given the 0 indexed column index
double* Matrix::GetColumnVec(int col) const {
  auto* column_vec = new double[nRows_];
  for (int row = 0; row < nRows_; ++row) {
    column_vec[row] = matrix_ptr_[row][col];
  }
  return column_vec;
}

// Get row vector given the 0 indexed row index
double* Matrix::GetRowVec(int row) const {
  auto* row_vec = new double[nCols_];
  for (int col = 0; col < nCols_; ++col) {
    row_vec[col] = matrix_ptr_[row][col];
  }
  return row_vec;
}
// COMPUTE DOT PRODUCT
double Matrix::Dot(double* vec1, double* vec2, int size) {
  double total = 0;
  for (int i = 0; i < size; ++i) {
    total += (vec1[i] * vec2[i]);
  }
  return total;
}



/******************************************************
 * GETTERS AND SETTERS 
 * ****************************************************
*/
void Matrix::SetElement(double data, int row, int col) {
  if (row > nRows_ - 1 || col > nCols_ - 1) {
    throw std::invalid_argument("Invalid dimensions");
  }
  matrix_ptr_[row][col] = data;
}

int Matrix::GetCols() const {
  return nCols_;
}

int Matrix::GetRows() const {
  return nRows_;
}


double** Matrix::GetMatrix() const {
  return matrix_ptr_;
}

/* ***************************
STREAM INSERTION OPERATOR TO DISPLAY MATRIX
*******************************
*/
std::ostream& operator<<(std::ostream& output, const Matrix& matrix) {
    for (int row = 0; row < matrix.nRows_; ++row) {
        for (int col = 0; col < matrix.nCols_; ++col) {
            output << matrix.matrix_ptr_[row][col] << " ";
        }
        output << std::endl;
    }
    output << std::endl;
    return output;
}

/**********************
 * SET TO IDENTITY
 * ***************************
*/
void Matrix::SetToIdentity() {
  if (!IsSquare()) { throw std::runtime_error("Only square matrices can be brought to the identity matrix"); }
  for (int row = 0; row < nRows_; ++row) {
    for (int col = 0; col < nCols_; ++col) {
      if (row == col) {
        matrix_ptr_[row][col] = 1;
      } else {
        matrix_ptr_[row][col] = 0;
      }
    }
  }
}

/*******************
 * CLOSE ENOUGH FUNCTION
*/
bool Matrix::CloseEnough(double num1, double num2) const {
  num1 = static_cast<double>(num1);
  num2 = static_cast<double>(num2);
  return std::abs(num1 - num2) < kTolerance;
}

/*********************
 * IS SQUARE FUNCTION
 * *******************
*/
bool Matrix::IsSquare() {
  return nCols_ == nRows_;
}

/*******************
 * CHECKS IF MATRIX IS IN ECHELON FORM 
 * ********************
*/
bool Matrix::IsEchelon() const {
  double cumulative_sum = 0.0;
  for (int row = 1; row < nRows_; ++row) {
    for (int col = 0; col < row; ++col) {
      cumulative_sum += matrix_ptr_[row][col];
    }
  }
  return CloseEnough(cumulative_sum, 0.0);
}

/**********************************
 * ROW OPERATIONS
 * ****************************
*/
void Matrix::RowSwap(int i, int j) {
  if (i < 0 || i > nRows_ - 1 || j < 0 || j > nRows_ - 1) { throw std::invalid_argument("Invalid rows"); }
  for (int col = 0; col < nCols_; ++col) {
    std::swap(matrix_ptr_[i][col], matrix_ptr_[j][col]);
  }
}
// Row replacement with 0 indexed rows; also the first parameter is the row that is scaled and added to the second parameter row
void Matrix::RowReplacement(int i, int j, double factor) {
  if (i < 0 || i > nRows_ - 1 || j < 0 || j > nRows_ - 1) { throw std::invalid_argument("Invalid rows"); }
  for (int col = 0; col < nCols_; ++col) {
    matrix_ptr_[j][col] += factor * matrix_ptr_[i][col];
  }
} 

void Matrix::ScaleRow(int i, double factor) {
  if (i < 0 || i > nRows_ - 1) { throw std::invalid_argument("Invalid row"); }
  for (int col = 0; col < nCols_; ++col) {
    matrix_ptr_[i][col] *= factor;
  }
}
 
/*********************
 * JOIN
 * ********************
*/

Matrix Matrix::Join(const Matrix& matrix2) {
  if (nRows_ != matrix2.nRows_) { throw std::invalid_argument("Must have same # of rows"); }
  Matrix combined_mat(nRows_, nCols_ + matrix2.nCols_);
  for (int row = 0; row < nRows_; ++row) {
    for (int col = 0; col < combined_mat.nCols_; ++col){
      if (col < nCols_) { 
        combined_mat.matrix_ptr_[row][col] = matrix_ptr_[row][col]; 
      } else {
        combined_mat.matrix_ptr_[row][col] = matrix2.matrix_ptr_[row][col - nCols_];
      }
    }
  }
  return combined_mat;
}

/******************************
 * MATRIX INVERSE
 * ******************************
*/
Matrix Matrix::Invert() {
  if (nRows_ != nCols_) { throw std::invalid_argument("Cannot invert non-square matrix"); }
  // Create ID matrix to augment with; create augmented matrix
  Matrix identity(nRows_, nCols_); identity.SetToIdentity();
  Matrix augmented = this->Join(identity);
  // Begin iterating through each column (or row)
  for (int col = 0; col < nRows_; ++col) {
    // Check if diagonal element of column is zero
    if (std::abs(augmented.matrix_ptr_[col][col]) < kTolerance) {
      bool pivot_found = false;
      // If diag is zero, loop down the column to find a non zero pivot 
      for (int pivot = col + 1; pivot < nRows_; ++pivot) {
        if (std::abs(augmented.matrix_ptr_[pivot][col]) > kTolerance) {
          augmented.RowSwap(col, pivot);
          pivot_found = true;
          break;
        }
      }
      // If no non zero pivot is found in the column, the matrix is singular
      if (!pivot_found) { throw std::invalid_argument("Matrix is signular"); }
    }
    // Scale current row s.t diagonal element is 1
    double factor = static_cast<double>(1) / augmented.matrix_ptr_[col][col];
    augmented.ScaleRow(col, factor);
    // Eliminate all other rows, making sure to skip the pivot row
    for (int row = 0; row < nCols_; ++row) {
      if (row == col) { continue; }
      augmented.RowReplacement(col, row, -1 * augmented.matrix_ptr_[row][col]);
    }
  }
  std::vector<double> to_return;
  for (int row = 0; row < augmented.nRows_; ++ row) {
    for (int col = nCols_; col < nCols_ * 2; ++col) {
      to_return.push_back(augmented.matrix_ptr_[row][col]);
    }
  }
  identity.Deallocate();
  augmented.Deallocate();
  Matrix inverse(nRows_, nCols_, to_return);
  return inverse;
}

/*************************
 * DETERMINANT
 * ***********************  
*/
double Matrix::Determinant() const {
  if (nRows_ != nCols_) { throw std::invalid_argument("Cannot determine determinant of non square matrix"); }
  // Clone matrix to obtain a matrix we can change
  Matrix clone = *this;
  // Track changes to determinant from scaling a row or row swap
  int change_tracker = 0;
  for (int col = 0; col < nCols_; ++col) {
    // Get non zero entries in the diagonal; if 0 is currently in the diagonal, find another row to swap with
    if (std::abs(clone.matrix_ptr_[col][col]) < kTolerance) {
      bool pivot_found = false;
      for (int row = col + 1; row < nRows_; ++row) {
        if (std::abs(clone.matrix_ptr_[row][col]) > kTolerance) {
          clone.RowSwap(col, row);
          change_tracker++;
          pivot_found = true;
          break;
        }
      }
      // No row with non zero in the current diagonal position was found, thus matrix is singular and determinant is 0
      if (!pivot_found) { return 0; }
    }
    // Iterate through all rows and row replace in order to make every entry below current disagonal 0, effectivly turning the matrix into upper triangular form
    for (int row = col + 1; row < nRows_; ++row) {
      clone.RowReplacement(col, row, -1 * (clone.matrix_ptr_[row][col] / clone.matrix_ptr_[col][col]));
    }
  }
  // Iterate through diagonal to compute determinant
  double determinant = 1.0;
  for (int col = 0; col < nCols_; ++col) {
    determinant *= clone.matrix_ptr_[col][col];
  }
  // If odd number of row swaps were performed, multiply determinant by -1 
  if (change_tracker % 2 != 0) { determinant *= -1; }
  clone.Deallocate();
  return determinant;
}

Matrix Matrix::Transpose() const {
  // Init transpose with dimensions
  Matrix transpose(nCols_,nRows_);
  // Loop through the rows of the matrix, cols of transpose
  for (int row = 0; row < nRows_; ++row) {
    // Get the current row vector which wiull be the column vector of the transpose at the same idx
    double* row_vec = this->GetRowVec(row);
    // Loop down the current column of the transpose at insert the corresponding entries from the row vec
    for (int transpose_row = 0; transpose_row < nCols_; ++transpose_row) {
      transpose.matrix_ptr_[transpose_row][row] = row_vec[transpose_row];
    }
    delete[] row_vec;
  }
  return transpose;
}
/******************************
 * FINISH RREF AND SOLVE METHODS AT A LATER DATE
 * ***************************************
*/
void Matrix::ReducedEchelon() {
  for (int col = 0; col < nCols_ - 1; ++col) {
    if (std::abs(matrix_ptr_[col][col]) < kTolerance) {
      bool pivot_found = false;
      for (int pivot_idx = col + 1; pivot_idx < nRows_; ++pivot_idx) {
        if (std::abs(matrix_ptr_[pivot_idx][col]) > kTolerance) {
          pivot_found = true;
          RowSwap(col, pivot_idx);
          break;
        }
      }
      if (!pivot_found) { throw std::invalid_argument("Matrix has no solution"); }
    }
    ScaleRow(col, static_cast<double>(1) / matrix_ptr_[col][col]);
    for (int row = 0; row < nRows_; ++row) {
      if (col == row) { continue; }
      RowReplacement(col, row, -1 * matrix_ptr_[row][col]);
    }
  }
} 

Vector Matrix::Solve(const Vector& vector) {
  if (nCols_ != nRows_) { throw std::invalid_argument("Unique solution cannot be found"); }
  int vector_dimensions = vector.GetDimensions();

  std::vector<double> vector_data = vector.GetData();

  Matrix vector_as_matrix(vector_dimensions, 1, vector_data);

  Matrix augmented = this->Join(vector_as_matrix);

  augmented.ReducedEchelon();
  std::cout << augmented << std::endl;

  std::vector<double> solution;
  
  for (int row = 0; row < nRows_; ++row) {
    solution.push_back(augmented.matrix_ptr_[row][nCols_]);
  }
  Vector result(solution);
  vector_as_matrix.Deallocate(); augmented.Deallocate();
  return result;
}



