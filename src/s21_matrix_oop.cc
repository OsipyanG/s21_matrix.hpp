#include "s21_matrix_oop.h"

#include <cmath>
#include <limits>

//______________________________________________________________________________
// GETTERS AND SETTERS

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows_ == rows) return;

  if (rows < 1) {
    throw std::invalid_argument(
        "SettingRowsError: The number of rows cannot be less than 1");
  }
  S21Matrix tmp(rows, cols_);
  int rows_range;
  if (rows > rows_) {
    rows_range = rows_;
  } else {
    rows_range = rows;
  }
  for (int i = 0; i < rows_range; ++i) {
    for (int j = 0; j < cols_; ++j) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  RemoveMatrix();
  *this = tmp;
}

void S21Matrix::SetCols(int cols) {
  if (cols_ == cols) return;
  if (cols < 1) {
    throw std::invalid_argument(
        "SettingColsError: The number of cols cannot be less than 1");
  }
  S21Matrix tmp(rows_, cols);
  int cols_range;
  if (cols > cols_) {
    cols_range = cols_;
  } else {
    cols_range = cols;
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_range; ++j) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  RemoveMatrix();
  *this = tmp;
}
//______________________________________________________________________________

// CONSTRUCTORS
S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  if (rows_ < 1 || cols_ < 1) {
    throw std::out_of_range(
        "CreationError: The number of rows or cols cannot be less than 1");
  }
  AlocateMemory();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(nullptr) {
  AlocateMemory();
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

//______________________________________________________________________________

// Destructors
S21Matrix::~S21Matrix() { RemoveMatrix(); }
//______________________________________________________________________________

// Operations
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Different dimensions of the matrices");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}
//
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Different dimensions of the matrices.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) >=
          std::numeric_limits<double>::epsilon()) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "To multiply matrices, the number of columns of the first matrix "
        "must "
        "be equal to the number of rows of the second matrix");
  }
  S21Matrix tmp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      tmp.matrix_[i][j] = 0;
      for (int k = 0; k < cols_; k++) {
        tmp.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = tmp;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix tmp(cols_, rows_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      tmp.matrix_[j][i] = matrix_[i][j];
    }
  }
  return tmp;
}

double S21Matrix::Determinant() {
  double result = 1;
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");
  int size = rows_;

  S21Matrix temp_matrix = S21Matrix(*this);

  for (int i = 0; i < size; ++i) {
    int k = i;
    for (int j = k + 1; j < size; ++j)
      if (fabs(temp_matrix.matrix_[k][i]) < fabs(temp_matrix.matrix_[j][i]))
        k = j;

    std::swap(temp_matrix.matrix_[i], temp_matrix.matrix_[k]);

    if (i != k) result = -result;
    result *= temp_matrix.matrix_[i][i];

    for (int j = i + 1; j < size; ++j) {
      temp_matrix.matrix_[i][j] /= temp_matrix.matrix_[i][i];
    }

    for (int j = 0; j < size; ++j) {
      if (j != i && fabs(temp_matrix.matrix_[j][i]) >=
                        std::numeric_limits<double>::epsilon())
        for (k = i + 1; k < size; k++)
          temp_matrix.matrix_[j][k] -=
              temp_matrix.matrix_[j][i] * temp_matrix.matrix_[i][k];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");

  S21Matrix complementsMatrix = S21Matrix(rows_, cols_);

  if (rows_ == 1 && cols_ == 1) {
    complementsMatrix.matrix_[0][0] = 1;
    return complementsMatrix;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      complementsMatrix.matrix_[i][j] = CalcComplement(i, j);
    }
  }
  return complementsMatrix;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");
  double determinant = Determinant();
  if (fabs(determinant) < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("The determinant is zero!");
  S21Matrix adjacent_matrix = CalcComplements().Transpose();
  adjacent_matrix.MulNumber(1 / determinant);
  return adjacent_matrix;
}

//______________________________________________________________________________
//
//// Operators
S21Matrix& S21Matrix::operator=(const S21Matrix& b) {
  if (this == &b) return *this;
  RemoveMatrix();
  rows_ = b.rows_;
  cols_ = b.cols_;
  AlocateMemory();
  for (int i = 0; i < rows_; ++i) {
    std::copy(b.matrix_[i], b.matrix_[i] + b.cols_, matrix_[i]);
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& b) {
  SumMatrix(b);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& b) {
  SubMatrix(b);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& b) {
  MulMatrix(b);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double& b) {
  MulNumber(b);
  return *this;
}

S21Matrix operator+(const S21Matrix& m1, const S21Matrix& m2) {
  S21Matrix tmp(m1);
  tmp.SumMatrix(m2);
  return tmp;
}
S21Matrix operator-(const S21Matrix& m1, const S21Matrix& m2) {
  S21Matrix tmp(m1);
  tmp.SubMatrix(m2);
  return tmp;
}
S21Matrix operator*(const S21Matrix& m1, const S21Matrix& m2) {
  S21Matrix tmp(m1);
  tmp.MulMatrix(m2);
  return tmp;
}
S21Matrix operator*(const S21Matrix& m1, const double& number) {
  S21Matrix tmp(m1);
  tmp.MulNumber(number);
  return tmp;
}

S21Matrix operator*(const double& number, const S21Matrix& m1) {
  S21Matrix tmp(m1);
  tmp.MulNumber(number);
  return tmp;
}
double& S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("InvalidIndexError: Index is out of range");
  }
  return matrix_[row][col];
}
double S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("InvalidIndexError: Index is out of range");
  }
  return matrix_[row][col];
}

//______________________________________________________________________________
// PRIVATE FUNCTIONS

void S21Matrix::AlocateMemory() {
  matrix_ = new double*[rows_];

  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
  }
}

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::RemoveMatrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}
double S21Matrix::CalcComplement(int i, int j) {
  S21Matrix minorMatrix = S21Matrix(rows_ - 1, cols_ - 1);
  for (int k = 0; k < rows_; ++k) {
    for (int l = 0; l < cols_; ++l) {
      if (k != i && l != j) {
        int new_i = k > i ? k - 1 : k;
        int new_j = l > j ? l - 1 : l;
        minorMatrix.matrix_[new_i][new_j] = matrix_[k][l];
      }
    }
  }
  double result = minorMatrix.Determinant() * pow(-1, i + j);
  return result;
}
