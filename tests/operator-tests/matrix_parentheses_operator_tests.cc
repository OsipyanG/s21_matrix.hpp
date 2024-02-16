#include "../s21_matrix_test.h"

TEST(operators_suite, parentheses_operators_test_1) {
  S21Matrix matrix = S21Matrix();
  matrix(1, 1) = 10.123;
  EXPECT_DOUBLE_EQ(((const S21Matrix)matrix)(1, 1), 10.123);
}
TEST(operators_suite, parentheses_operator_test_2) {
  S21Matrix matrix = S21Matrix(3, 4);
  EXPECT_THROW(matrix(0, 4) = 10.132, std::out_of_range);
  const S21Matrix const_matrix = S21Matrix(4, 3);
  EXPECT_THROW(const_matrix(0, 4), std::out_of_range);
}

TEST(operators_suite, parentheses_operator_test_3) {
  S21Matrix matrix = S21Matrix(4, 3);
  EXPECT_THROW(matrix(4, 0) = 10.132, std::out_of_range);
  const S21Matrix const_matrix = S21Matrix(4, 3);
  EXPECT_THROW(const_matrix(4, 0), std::out_of_range);
}