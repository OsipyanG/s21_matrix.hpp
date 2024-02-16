cmake_minimum_required(VERSION 3.10)

project(Matrix++ VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_STATIC_LIBRARY_PREFIX "")
set(CMAKE_BUILD_TYPE Release)

if(CMAKE_BUILD_TYPE EQUAL "Debug")
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

add_library(s21_matrix_oop STATIC src/s21_matrix_oop.cc)
target_compile_options(s21_matrix_oop PRIVATE -Wall -Werror -Wextra -Wpedantic)

find_package(GTest REQUIRED)
include(GoogleTest)
enable_testing()

add_executable(
        tests
        src/s21_matrix_oop.cc
        tests/other-tests/matrix_mutators_tests.cc
        tests/unit_tests.cc
        tests/other-tests/matrix_constructors_tests.cc
        tests/other-tests/matrix_accessor_tests.cc
        tests/operator-tests/matrix_assignment_operator_tests.cc
        tests/operator-tests/matrix_mult_matrix_assignment_operator_tests.cc
        tests/operator-tests/matrix_mult_matrix_operator_tests.cc
        tests/operator-tests/matrix_mult_number_assignment_operator_tests.cc
        tests/operator-tests/matrix_mult_number_operator_tests.cc
        tests/operator-tests/matrix_parentheses_operator_tests.cc
        tests/operator-tests/matrix_sub_assignment_operator_tests.cc
        tests/operator-tests/matrix_sub_operator_tests.cc
        tests/operator-tests/matrix_sum_assignment_operator_tests.cc
        tests/operator-tests/matrix_sum_operator_tests.cc
        tests/function-tests/matrix_calc_complements_function_tests.cc
        tests/function-tests/matrix_determinant_function_tests.cc
        tests/function-tests/matrix_eq_function_tests.cc
        tests/function-tests/matrix_inverse_function_tests.cc
        tests/function-tests/matrix_mult_matrix_function_tests.cc
        tests/function-tests/matrix_mult_number_function_tests.cc
        tests/function-tests/matrix_sub_function_tests.cc
        tests/function-tests/matrix_sum_function_tests.cc
        tests/function-tests/matrix_transpose_function_tests.cc
        tests/operator-tests/matrix_mult_number_operator_tests.cc)
add_executable(my_test src/s21_matrix_oop.cc)
target_compile_options(tests PRIVATE -Wall -Werror -Wextra -Wpedantic)

target_link_libraries(tests GTest::gtest_main s21_matrix_oop)

gtest_discover_tests(tests)
