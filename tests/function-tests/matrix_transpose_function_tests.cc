#include "../s21_matrix_test.h"

TEST(functions_suite, transpose_function_test_1) {
  double arr1[3][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};

  double expectedValues[4][3] = {{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}};

  S21Matrix A = S21Matrix(3, 4);
  fill_matrix(A, arr1);

  S21Matrix expected_result = S21Matrix(4, 3);
  fill_matrix(expected_result, expectedValues);

  S21Matrix result = A.Transpose();
  EXPECT_TRUE(expected_result == result);
}
