#include "../s21_matrix_test.h"
TEST(accessors_suite, get_rows_test_1) {
  S21Matrix result = S21Matrix(4, 4);
  EXPECT_TRUE(result.GetRows() == 4);
}
TEST(accessors_suite, get_cols_test_1) {
  S21Matrix result = S21Matrix(4, 4);
  EXPECT_TRUE(result.GetCols() == 4);
}