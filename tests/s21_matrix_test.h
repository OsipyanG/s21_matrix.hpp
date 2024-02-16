#ifndef S21_MATRIX_TEST_H
#define S21_MATRIX_TEST_H
#include <gtest/gtest.h>

#include "../src/s21_matrix_oop.h"
#define fill_matrix(m, arr_values)          \
  for (int i = 0; i < m.GetRows(); i++) {   \
    for (int j = 0; j < m.GetCols(); j++) { \
      m(i, j) = arr_values[i][j];           \
    }                                       \
  }
#endif
