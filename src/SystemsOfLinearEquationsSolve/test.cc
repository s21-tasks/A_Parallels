#include "SLAEGauss.h"
#include "matrix.h"

#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>

template<typename T>
bool compare_results(s21::Matrix<T> m, std::vector<T> v, std::vector<T> result) {
  s21::SLAEGauss gauss(m,v);
  std::vector<T> sol = gauss.UsualExecute();
    if (sol != result)
      return false;
    for (int i = 1; i <= 10; ++i) {
      sol = gauss.ParallelExecute(i);
      if (sol != result)
        return false;
    }
  return true;
}

void gauss_test() {
  s21::Matrix<double> m{{2,-1,0},
                        {-1,1,4},
                        {1,2,3}};
  std::vector<double> v{0,13,14};
  std::vector<double> result{1,2,3};
}

void gauss_test_4x4() {
  s21::Matrix<double> m{{5,3,2,-8},
                        {1,1,1,1},
                        {3,5,1,4},
                        {4,2,3,1}};
  std::vector<double> v{1,0,0,3};

}


int main() {
  s21::Matrix<double> m{{2,-1,0},
                        {-1,1,4},
                        {1,2,3}};
  std::vector<double> v{0,13,14};
  std::vector<double> result{1,2,3};
  if (compare_results(m,v,result))
    std::cout  << "1st test ok"<< std::endl;

  s21::Matrix<double> m2{{5,3,2,-8},
                        {1,1,1,1},
                        {3,5,1,4},
                        {4,2,3,1}};
  std::vector<double> v2{1,0,0,3};
  std::vector<double> result2{7,-8,-5,6};
  if (compare_results(m2,v2,result2))
    std::cout  << "2nd test ok"<< std::endl;


  return 0;
}
