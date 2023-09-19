

#include "../matrix/blas/blas_matrix.h"
#include "../utility/m_random.h"
#include "../utility/m_time.h"
#include "m_thread.h"

int main() {
  int thread_count = 4;

  BLAS::Matrix<float> M1D(
      100, 500, [&] { return s21::Random::Normal<float>(0.0, 15.0); });
  BLAS::Matrix<float> M1C(M1D);
  BLAS::Matrix<float> M1(M1D);
  BLAS::Matrix<float> M2(1, 500,
                         [&] { return s21::Random::Normal<float>(0.0, 5.0); });

  ThreadManager TM(thread_count);

  auto cr = s21::Time::Compare<s21::Time::mcs>(
      10000,
      [&] {
        for (int k = 0; k < M1.GetRows(); ++k) {
          for (int g = 0; g < M1.GetCols(); ++g) {
            M1(k, g) -= 2 * M2(0, g);
          }
        }
      },
      [&] {
        ThreadManager::DispThreads(
            thread_count,
            ThreadManager::LoopThreads(thread_count, M1D.GetRows(), [&](int k) {
              for (int g = 0; g < M1D.GetCols(); ++g) {
                M1D(k, g) -= 2 * M2(0, g);
              }
            }));
      },
      [&] {
        TM.LoopExecute(M1C.GetRows(), [&](int k) {
          for (int g = 0; g < M1C.GetCols(); ++g) {
            M1C(k, g) -= 2 * M2(0, g);
          }
        });
      });

  for (auto i : cr) {
    std::cout << i << '\n';
  }

  std::cout << (M1 == M1D) << (M1 == M1C) << '\n';

  return 0;
}
