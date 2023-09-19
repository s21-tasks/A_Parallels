#include "../../sub/thread/m_thread.h"
#include "tests.h"

using namespace s21;

template <class T, class Unit>
void Test(unsigned int lines, unsigned int matrix_size, unsigned int repeat) {
  static int counter = 0;
  std::vector<Matrix<T>> tasks(repeat, Matrix<T>(matrix_size));
  for (auto &t : tasks) t.Fill(RS<T>);

  std::vector<Matrix<T>> workers(lines, Matrix<T>(matrix_size));
  for (auto &w : workers) w.Fill(RS<T>);

  std::cout << "\033[1;32m"
            << "\nTest " << ++counter << " {Matrix size = " << matrix_size
            << ". Lines = " << lines << "}:\n"
            << "\033[1;37m";

  auto time_point = std::chrono::high_resolution_clock::now();
  {
    PipeLine::Winograd<T> W(workers, matrix_size);
    for (auto &t : tasks) {
      W.Execute(t);
    }
  }
  std::cout << '\t' << Time::Duration<Unit>(time_point) / repeat << ' '
            << Time::Prefix<Unit>::val << "\t - Pipeline Winograd\n";

  time_point = std::chrono::high_resolution_clock::now();
  {
    Parallel::Winograd<T> W(matrix_size);
    for (auto &t : tasks) {
      for (auto &w : workers) {
        W.Execute(t, w, t);
      }
    }
  }
  std::cout << '\t' << Time::Duration<Unit>(time_point) / repeat << ' '
            << Time::Prefix<Unit>::val << "\t - Parallel Winograd\n";

  time_point = std::chrono::high_resolution_clock::now();
  {
    Basic::Winograd<T> W(matrix_size);
    for (auto &t : tasks) {
      for (auto &w : workers) {
        W.Execute(t, w, t);
      }
    }
  }
  std::cout << '\t' << Time::Duration<Unit>(time_point) / repeat << ' '
            << Time::Prefix<Unit>::val << "\t - Basic Winograd\n";

  time_point = std::chrono::high_resolution_clock::now();
  {
    ThreadManager TM(4);
    std::shared_ptr<Basic::Winograd<T>> W_ptr =
        std::make_shared<Basic::Winograd<T>>(matrix_size);

    TM.LoopExecute(tasks.size(), [W_ptr, &workers, &tasks](int k) {
      for (auto &w : workers) {
        W_ptr->Execute(tasks[k], w, tasks[k]);
      }
    });
  }
  std::cout << '\t' << Time::Duration<Unit>(time_point) / repeat << ' '
            << Time::Prefix<Unit>::val
            << "\t - Basic Winograd (ThreadManager)\n";
}

int main() {
  Test<float, Time::ns>(4, 4, 100000);
  Test<float, Time::ns>(4, 5, 100000);
  Test<float, Time::ns>(4, 6, 100000);
  Test<float, Time::ns>(4, 7, 100000);
  Test<float, Time::ns>(4, 8, 100000);
  // Test<float, Time::ns>(4, 9, 100000);
  // Test<float, Time::ns>(4, 10, 100000);
  // Test<float, Time::ns>(4, 12, 100000);
  // Test<float, Time::ns>(4, 14, 100000);
  Test<float, Time::ns>(4, 16, 100000);

  return 0;
}
