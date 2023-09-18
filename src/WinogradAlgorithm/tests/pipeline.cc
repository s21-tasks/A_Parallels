#include "tests.h"

using namespace s21;

template <class T, class Unit>
void Test(unsigned int lines, unsigned int matrix_size, unsigned int repeat) {
  static int counter = 0;
  std::vector<Matrix<T>> tasks(repeat, Matrix<T>(matrix_size));
  for (auto &t : tasks) t.Fill(RS<T>);

  std::vector<Matrix<T>> workers(lines, Matrix<T>(matrix_size));
  for (auto &w : workers) w.Fill(RS<T>);

  PipeLine::Winograd<T> pipe(workers, matrix_size);
  Basic::Winograd<T> basic(matrix_size);
  Parallel::Winograd<T> par(matrix_size);

  unsigned int k1 [[maybe_unused]] = 0, k2 [[maybe_unused]] = 0,
                  k3 [[maybe_unused]] = 0;
  auto result = Time::Compare<Unit>(
      repeat, [&] { pipe.Execute(tasks[k1++]); },
      [&] {
        for (auto &w : workers) {
          basic.Mul(tasks[k2], w, tasks[k2]);
        }
        ++k2;
      },
      [&] {
        for (auto &w : workers) {
          par.Mul(tasks[k3], w, tasks[k3]);
        }
        ++k3;
      });

  std::cout << "\033[1;32m";
  std::cout << "\nTest " << ++counter << " {Matrix size = " << matrix_size
            << ". Lines = " << lines << "}:\n";
  std::cout << "\033[1;37m";
  std::cout << '\t' << result[0] << ' ' << Time::Prefix<Unit>::val
            << "\t - Pipeline Winograd\n";
  std::cout << '\t' << result[1] << ' ' << Time::Prefix<Unit>::val
            << "\t - Basic Winograd\n";
  std::cout << '\t' << result[2] << ' ' << Time::Prefix<Unit>::val
            << "\t - Parallel\n";
}

int main() {
  // Test<float, Time::ns>(4, 8, 50000);
  // Test<float, Time::ns>(4, 25, 5000);
  // Test<float, Time::mcs>(4, 64, 5000);
  // Test<float, Time::mcs>(4, 78, 500);
  // Test<float, Time::mcs>(4, 128, 500);
  // Test<float, Time::mcs>(4, 133, 500);
  // Test<float, Time::mcs>(4, 211, 500);
  // Test<float, Time::mcs>(4, 255, 200);
  // Test<float, Time::mcs>(4, 256, 150);
  // Test<float, Time::ms>(4, 512, 50);
  // Test<float, Time::ms>(4, 1024, 3);

  std::vector<Matrix<float>> tasks(100, Matrix<float>(5));
  for (auto &t : tasks) t.Fill(RS<float>);

  std::vector<Matrix<float>> workers(4, Matrix<float>(5));
  for (auto &w : workers) w.Fill(RS<float>);

  PipeLine::Winograd<float> pipe(workers, 5);
  Basic::Winograd<float> basic(5);
  // Parallel::Winograd<float> par(5);

  for (auto &t : tasks) {
    pipe.Execute(t);
    // basic.Mul(t, workers[0], t);
  }

  // sleep
  std::this_thread::sleep_for(std::chrono::seconds(5));

  return 0;
}
