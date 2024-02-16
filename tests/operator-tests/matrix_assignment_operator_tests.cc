#include "../s21_matrix_test.h"

TEST(operators_suite, assignment_operator_test_1) {
  double arr1[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  double arr2[4][4] = {{1, 2, 4, 5},
                       {1, 5, 5, 6},
                       {4.1434, 631.31636, 135.315, 136.148689},
                       {1.5, 616.527, 16.13, 6426.64246}};

  S21Matrix A = S21Matrix(3, 3);
  fill_matrix(A, arr1);

  S21Matrix result = S21Matrix(4, 4);
  fill_matrix(result, arr2);

  S21Matrix expectedResult = S21Matrix(3, 3);
  fill_matrix(expectedResult, arr1);

  result = A;
  EXPECT_TRUE(expectedResult == result);
  result(0, 0) = 0;
  EXPECT_NE(result(0, 0), expectedResult(0, 0));
}

TEST(operators_suite, assignment_operator_test_2) {
  double arr1[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  S21Matrix A = S21Matrix(3, 3);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(3, 3);
  fill_matrix(expectedResult, arr1);

  S21Matrix result;
  result = A;
  EXPECT_TRUE(expectedResult == result);
}

TEST(operators_suite, assignment_operator_test_3) {
  double arr1[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};

  S21Matrix A = S21Matrix(3, 3);
  fill_matrix(A, arr1);

  S21Matrix expectedResult = S21Matrix(3, 3);
  fill_matrix(expectedResult, arr1);

  S21Matrix result = A;
  EXPECT_TRUE(expectedResult == result);
}
