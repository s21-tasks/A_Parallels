#include "matrix.h"
#include "blas.h"
#include <gtest/gtest.h>
#include "../utility/m_random.h"

using namespace s21;

TEST(Matrix, Constructor) {
  Matrix<int> m(2, 3);
  EXPECT_EQ(2, m.GetRows());
  EXPECT_EQ(3, m.GetCols());
}

TEST(Matrix, Geters) {
    Matrix<int> m(2, 3, [] { return Random::Int(0, 100); });
    EXPECT_EQ(m(1, 2), m[1][2]);
}

TEST(Matrix, mul1) {
    Matrix<int> a(10, 10, [] { return Random::Int(0, 100); });
    Matrix<int> b(10, 10, [] { return Random::Int(0, 100); });
    auto c1 = a * b;
    Matrix<int> c2(10, 10);
    Matrix<int>::MulParallel(a, b, c2, 4);
    EXPECT_EQ(c1, c2);
}

TEST(Matrix, mul_blas) {
    Matrix<float> a(10, 10, [] { return Random::Int(0, 100); });
    Matrix<float> b(10, 10, [] { return Random::Int(0, 100); });
    auto c1 = a * b;
    Matrix<float> c2(10, 10);
    BLAS<float>::Mul(a, b, c2);
    EXPECT_EQ(c1, c2);
}

