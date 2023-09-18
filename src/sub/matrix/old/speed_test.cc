#include <vector>

#include "../utility/random.h"
#include "../utility/str.h"
#include "../utility/time.h"
#include "blas/blas_matrix.h"
#include "openCL/cl_matrix.h"

static void one_output(int count, const std::string &name, int64_t time) {
  std::cout << name << ": " << time
            << ", 1 iter: " << (double)time / (double)count << " ms\n";
}

static void output(const std::string &name, int count,
                   std::pair<int64_t, int64_t> res) {
  std::cout << "\n" << name << " count = " << count << ":\n";
  one_output(count, "blas", res.first);
  one_output(count, "ocl ", res.second);
}

template <class Type>
static void ConstructorVector(const std::string &name, int count, int cols,
                              int rows) {
  std::vector<Type> vec(cols * rows);
  for (auto &i : vec) {
    i = s21::Random::Normal<Type>(0, 1000);
  }
  output(s21::Str::Fill<' '>("Constructor vector", name, '(', cols, 'x', rows,
                             ')'),
         count,
         s21::Time::Compare([&] { BLAS::Matrix<Type>(cols, rows, vec); },
                            [&] { CL::Matrix<Type>(cols, rows, vec); }, count));
}

template <class Type>
static void Constructor(const std::string &name, int count, int cols,
                        int rows) {
  output(s21::Str::Fill<' '>("Constructor", name, '(', cols, 'x', rows, ')'),
         count,
         s21::Time::Compare([&] { BLAS::Matrix<Type>(cols, rows); },
                            [&] { CL::Matrix<Type>(cols, rows); }, count));
}

template <class Type>
static void MulTest(const std::string &name, int m, int n, int k, int count) {
  BLAS::Matrix<Type> BMres(m, n);
  std::vector<Type> v1(m * k);
  for (auto &i : v1) {
    i = s21::Random::Normal<Type>(0, 1000);
  }
  std::vector<Type> v2(k * n);
  for (auto &i : v2) {
    i = s21::Random::Normal<Type>(0, 1000);
  }
  CL::Matrix<Type> CMres(m, n);
  output(s21::Str::Fill<' '>("Mul", name, m, 'x', k, "mul", k, 'x', n), count,
         s21::Time::Compare(
             [&] {
               BLAS::Arithmetic::Mul(BLAS::Matrix<Type>(m, k, v1),
                                     BLAS::Matrix<Type>(k, n, v2), BMres);
             },
             [&] {
               CL::Arithmetic::Mul(CL::Matrix<Type>(m, k, v1),
                                   CL::Matrix<Type>(k, n, v2), CMres);
             }));
}

int main() {
  Constructor<float>("float", 1, 10000, 10000);
  Constructor<double>("double", 1, 10000, 10000);
  ConstructorVector<float>("float", 1, 1000, 10000);
  ConstructorVector<double>("double", 1, 1000, 10000);
  Constructor<float>("float", 100, 10000, 1000);
  Constructor<double>("double", 100, 10000, 1000);
  ConstructorVector<float>("float", 100, 1000, 1000);
  ConstructorVector<double>("double", 100, 1000, 1000);
  Constructor<float>("float", 10000, 100, 5);
  Constructor<double>("double", 10000, 50, 10);
  ConstructorVector<float>("float", 10000, 3, 40);
  ConstructorVector<double>("double", 10000, 4, 60);

  MulTest<float>("float", 5, 5, 5, 10000);
  MulTest<float>("float", 50, 50, 50, 2000);
  MulTest<float>("float", 50, 50, 50, 2000);
  MulTest<float>("float", 500, 500, 500, 300);
  MulTest<float>("float", 1000, 1000, 1000, 30);
  MulTest<float>("float", 1000, 10000, 100, 30);
  MulTest<float>("float", 10000, 10000, 100, 10);
  MulTest<float>("float", 1000, 1000, 1000, 10);
  MulTest<float>("float", 1, 10000, 10000, 10);
  MulTest<float>("float", 10000, 1, 10000, 10);
  MulTest<float>("float", 1, 1, 10000000, 10);
  MulTest<float>("float", 1, 10000000, 1, 10);

  return 0;
}