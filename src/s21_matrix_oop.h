#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;  // Rows and columns
  double** matrix_;
  void AlocateMemory();
  void CopyMatrix(const S21Matrix& other);
  void RemoveMatrix();
  double CalcComplement(int i, int j);

 public:
  // Getters and setters

  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows);
  void SetCols(int cols);

  // Constructors
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;

  // Destructor
  ~S21Matrix();

  // Operators
  bool operator==(const S21Matrix& b) { return EqMatrix(b); }
  S21Matrix& operator=(const S21Matrix& b);
  S21Matrix& operator+=(const S21Matrix& b);
  S21Matrix& operator-=(const S21Matrix& b);
  S21Matrix& operator*=(const S21Matrix& b);
  S21Matrix& operator*=(const double& b);

  double& operator()(int row, int col);
  double operator()(int row, int col) const;
  // Matrix operations
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();
};

std::ostream& operator<<(std::ostream& os, const S21Matrix& d);

S21Matrix operator+(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator-(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator*(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator*(const S21Matrix& m1, const double& number);
S21Matrix operator*(const double& number, const S21Matrix& m1);
#endif  // S21_MATRIX_OOP_H
