#include "../s21_matrix_test.h"
TEST(constructors_suite, default_constructor_1) {
  S21Matrix matrix = S21Matrix();
  EXPECT_EQ(matrix.GetRows(), 3);
  EXPECT_EQ(matrix.GetCols(), 3);
}
TEST(constructors_suite, constructor_with_matrix_size_1) {
  S21Matrix matrix = S21Matrix(5, 6);
  EXPECT_EQ(matrix.GetRows(), 5);
  EXPECT_EQ(matrix.GetCols(), 6);
}
TEST(constructors_suite, constructor_with_matrix_size_2) {
  EXPECT_THROW(S21Matrix(-1, -3), std::out_of_range);
}
TEST(constructors_suite, copy_constructor_1) {
  double arr[4][3] = {
      {1.52, 2.135, 3.315},
      {4.316316136, 5.136136, 6.31531},
      {7.136, 8.1614, 9.64314},
      {3.341134, 10.136133, 136.36136},

  };
  int rows = 4;
  int cols = 3;
  S21Matrix matrix = S21Matrix(rows, cols);
  fill_matrix(matrix, arr);

  S21Matrix expected_matrix = S21Matrix(rows, cols);
  fill_matrix(expected_matrix, arr);

  S21Matrix result = S21Matrix(matrix);

  EXPECT_TRUE(result == expected_matrix);
}

TEST(constructors_suite, move_constructor_1) {
  double arr[4][3] = {
      {1.52, 2.135, 3.315},
      {4.316316136, 5.136136, 6.31531},
      {7.136, 8.1614, 9.64314},
      {3.341134, 10.136133, 136.36136},

  };
  int rows = 4;
  int cols = 3;
  S21Matrix matrix = S21Matrix(rows, cols);
  fill_matrix(matrix, arr);

  S21Matrix expected_matrix = S21Matrix(rows, cols);
  fill_matrix(expected_matrix, arr);

  S21Matrix result = S21Matrix(std::move(matrix));

  EXPECT_TRUE(result == expected_matrix);
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetCols(), 0);
}