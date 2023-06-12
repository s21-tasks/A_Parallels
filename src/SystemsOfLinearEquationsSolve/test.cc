#include "SLAEGauss.h"
#include "matrix.h"

#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>

int main() {
  s21::Matrix<double> m{{1,2,3},
                        {3,5,7},
                        {1,3,4}};
  std::vector<double> v{3,0,1};
  s21::SLAEGauss gauss(m,v);
  auto sol = gauss.UsualExecute();
  for (auto a = sol.begin();a!=sol.end();++a) {
    std::cout << *a << " ";
  }
  std::cout << std::endl;
  return 0;
}