
#include "../s21_matrix_test.h"

TEST(functions_suite, calc_complements_function_test_1) {
  const int rowsCount = 3;
  const int columnsCount = 3;
  double arr1[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  double expectedValues[3][3] = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(rowsCount, columnsCount);
  fill_matrix(expectedResult, expectedValues);

  S21Matrix result = A.CalcComplements();
  EXPECT_TRUE(expectedResult == result);
}

TEST(functions_suite, calc_complements_function_test_2) {
  const int rowsCount = 3;
  const int columnsCount = 4;
  double arr1[3][4] = {
      {10.4, 8.53, -9.321, 0}, {3.0, -10.6, 5.6, 0}, {55, -53.7, 123, 0}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  EXPECT_THROW(A.CalcComplements(), std::logic_error);
}

TEST(functions_suite, calc_complements_function_test_3) {
  const int rowsCount = 1;
  const int columnsCount = 1;
  double arr1[1][1] = {{3.235153}};

  double expectedValues[1][1] = {{1}};

  S21Matrix A = S21Matrix(rowsCount, columnsCount);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(rowsCount, columnsCount);
  fill_matrix(expectedResult, expectedValues);

  S21Matrix result = A.CalcComplements();
  EXPECT_TRUE(expectedResult == result);
}
