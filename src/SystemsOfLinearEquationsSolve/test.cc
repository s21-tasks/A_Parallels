#include "SLAEGauss.h"
#include "matrix.h"

#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>

#include "../sub/thread_pool/thread_pool.h"

void gauss_test() {
  s21::Matrix<double> m{{2,-1,0},
                        {-1,1,4},
                        {1,2,3}};
  std::vector<double> v{0,13,14};
  s21::SLAEGauss gauss(m,v,2);
  auto sol = gauss.UsualExecute();
  for (double & a : sol) {
    std::cout << a << " ";
  }
  std::cout << std::endl;

  auto sol2 = gauss.ParallelExecute();
  for (double & a : sol2) {
    std::cout << a << " ";
  }
  std::cout << std::endl;

}

void pool_test_func(const int number = 0, int sleep = 0) {
  std::cout << "func " << number << " start running at " << std::this_thread::get_id() << std::endl;
  std::this_thread::sleep_for(std::chrono::seconds(sleep));
  std::cout << "func " << number << " end running" << std::endl;
}


int main() {
//  s21::ThreadPool pool(4);
//  int num = 10;
//  for (int i = 0; i < num; ++i) {
//    int r = i;
//    pool.AddTask([i] { return pool_test_func(i, 0); });
//  }
//  pool.WaitForComplete();
//  std::cout << "all end" << std::endl;
  for (int i = 0; i < 50; ++i) {
    gauss_test();
    std::cout << i << std::endl;
  }
  return 0;
}