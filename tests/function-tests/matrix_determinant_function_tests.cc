
#include "../s21_matrix_test.h"
TEST(functions_suite, determinant_function_test_1) {
  double arr1[1][1] = {{10}};

  double expectedDeterminant = 10;

  S21Matrix A = S21Matrix(1, 1);
  fill_matrix(A, arr1);

  EXPECT_DOUBLE_EQ(A.Determinant(), expectedDeterminant);
}

TEST(functions_suite, determinant_function_test_2) {
  double arr1[2][2] = {{10, 7}, {4, -2}};

  double expectedDeterminant = -48;

  S21Matrix A = S21Matrix(2, 2);
  fill_matrix(A, arr1);

  EXPECT_DOUBLE_EQ(A.Determinant(), expectedDeterminant);
}

TEST(functions_suite, determinant_function_test_3) {
  double arr1[3][3] = {
      {10.4, 8.53, -9.321}, {3.0, -10.6, 5.6}, {55, -53.7, 123}};

  double expectedDeterminant = -14884.8919;

  S21Matrix A = S21Matrix(3, 3);
  fill_matrix(A, arr1);

  EXPECT_DOUBLE_EQ(A.Determinant(), expectedDeterminant);
}

TEST(functions_suite, determinant_function_test_4) {
  double arr1[3][4] = {
      {10.4, 8.53, -9.321, 0}, {3.0, -10.6, 5.6, 0}, {55, -53.7, 123, 0}};

  S21Matrix A = S21Matrix(3, 4);
  fill_matrix(A, arr1);

  EXPECT_THROW(A.Determinant(), std::logic_error);
}