#include <gtest/gtest.h>

#include "../../sub/utility/m_random.h"

#ifdef BASIC_WINOGRAD
#include "../winograd.h"
#define PREFIX Basic
#define FUNCTIONAL_CLASS BASIC__functional
#define ERROR_CLASS BASIC__error
#elif defined PARALLEL_WINOGRAD
#include "../winograd_parallel.h"
#define PREFIX Parallel
#define FUNCTIONAL_CLASS PARALLEL__functional
#define ERROR_CLASS PARALLEL__error
#endif

using namespace s21;

template <class T, class W>
void Helper(unsigned int N, T random_min, T random_max, const W &w) {
  Matrix<T> A(N, N, [&] { return Random::Easy<T>::R(random_min, random_max); });
  Matrix<T> B(N, N, [&] { return Random::Easy<T>::R(random_min, random_max); });
  Matrix<T> C1(N, N);
  auto C2 = A * B;
  w.Mul(A, B, C1);
  Matrix<T>::fp_compare_precision = 1e-5;
  EXPECT_EQ(C1, C2);
}

template <class T>
void Functional(unsigned int N, T random_min = static_cast<T>(-10),
                T random_max = static_cast<T>(10)) {
  Helper<T>(N, random_min, random_max, PREFIX::Winograd<T>(N));
}

#define __ONE_SIZE_TEST(T, n) \
  TEST(FUNCTIONAL_CLASS, __##n##x##n) { Functional<T>(n); }

#define __ERROR_TEST(N, M, K, L)                             \
  TEST(ERROR_CLASS, __##N##x##M##_dot_##K##x##L) {           \
    Matrix<float> A(N, M);                                   \
    Matrix<float> B(K, L);                                   \
    Matrix<float> C(N, L);                                   \
    EXPECT_ANY_THROW(PREFIX::Winograd<float>::Mul(A, B, C)); \
  }

__ONE_SIZE_TEST(double, 2)
__ONE_SIZE_TEST(long double, 3)
__ONE_SIZE_TEST(int, 4)
__ONE_SIZE_TEST(short, 5)
__ONE_SIZE_TEST(unsigned short, 6)
__ONE_SIZE_TEST(long, 7)
__ONE_SIZE_TEST(long long, 8)
__ONE_SIZE_TEST(unsigned int, 9)
__ONE_SIZE_TEST(unsigned long, 10)
__ONE_SIZE_TEST(unsigned long long, 11)
__ONE_SIZE_TEST(char, 12)
__ONE_SIZE_TEST(unsigned char, 13)
__ONE_SIZE_TEST(float, 14)

TEST(FUNCTIONAL_CLASS, __from_15x15_to_64x64) {
  for (unsigned int k = 15; k <= 64; ++k) {
    Functional<double>(k, Random::Normal(0.0, 15.5),
                       Random::Normal(-0.5, 1579.9));
  }
}

TEST(FUNCTIONAL_CLASS, __from_65x65_to_128x128) {
  for (unsigned int k = 65; k <= 128; ++k) {
    Functional<int>(k);
  }
}

__ERROR_TEST(5, 5, 6, 6)
__ERROR_TEST(5, 6, 6, 5)
__ERROR_TEST(5, 6, 6, 7)
__ERROR_TEST(5, 6, 7, 8)
__ERROR_TEST(1, 2, 3, 4)
__ERROR_TEST(1, 1, 1, 1)
__ERROR_TEST(0, 0, 0, 0)
