#include "../s21_matrix_test.h"

TEST(functions_suite, inverse_function_test_1) {
  const int rowsCount = 3;
  const int columnsCount = 3;
  double arr1[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  double expectedValues[3][3] = {
      {0, -0.1, 0.2}, {-0.25, 0.35, 0.05}, {0.5, -0.2, -0.1}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(rowsCount, columnsCount);
  fill_matrix(expectedResult, expectedValues);

  S21Matrix result = A.InverseMatrix();

  EXPECT_TRUE(result == expectedResult);
}

TEST(functions_suite, inverse_function_test_2) {
  const int rowsCount = 1;
  const int columnsCount = 1;
  double arr1[1][1] = {{52}};

  double expectedValues[1][1] = {{1.0 / 52}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(rowsCount, columnsCount);
  fill_matrix(expectedResult, expectedValues);

  S21Matrix result = A.InverseMatrix();

  EXPECT_TRUE(result == expectedResult);
}

TEST(functions_suite, inverse_function_test_3) {
  const int rowsCount = 3;
  const int columnsCount = 4;
  double arr1[3][4] = {{1, 2, 3, 4}, {0, 4, 2, 5}, {5, 2, 1, 6}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  EXPECT_THROW(A.InverseMatrix(), std::logic_error);
}

TEST(functions_suite, inverse_function_test_4) {
  const int rowsCount = 2;
  const int columnsCount = 2;
  double arr1[2][2] = {{1, 1}, {1, 1}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  EXPECT_THROW(A.InverseMatrix(), std::runtime_error);
}
