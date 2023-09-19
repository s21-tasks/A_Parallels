#pragma once

#include "../../sub/ConsoleInterface/ConsoleInterface.h"
#include "../../sub/utility/m_random.h"
#include "../../sub/utility/m_sstr.h"
#include "../../sub/utility/m_time.h"
#include "../SLAEGauss.h"

namespace s21 {

class GaussSpeedTest : public ConsoleInterface {
 public:
  GaussSpeedTest() {
    auto home_menu = AddMenu(
        {"loading the blas from console", "set random plas", "print blas",
         "solve blas using usual method", "solve  blas using parallel method",
         "compare speed from the methods"});

    home_menu->Connect(
        1, AddInput(new OneValueInput<int>(
               this, home_menu, [&](int &capacity) { RandomFill(capacity); },
               "capacity")));

    home_menu->Connect(2, home_menu, [&] { gauss.Print(); });

    home_menu->Connect(3, home_menu, [&] { UsualSolve(); });

    home_menu->Connect(
        4, AddInput(new OneValueInput<int>(
               this, home_menu,
               [&](int &threads_count) { ParallelSolve(threads_count); },
               "count of threads")));
    home_menu->Connect(5, home_menu, [&] { SpeedTest(); });
  }

  void SpeedTest() {
    for (int i = 1; i <= 10; ++i) {
      int64_t parallel_time =
          Time::Test([&] { return gauss.ParallelExecute(i); });
      std::cout << "Parralel time with " << i << " threads: " << parallel_time
                << std::endl;
    }
    int64_t usual_time = Time::Test([&] { return gauss.UsualExecute(); });
    std::cout << "Usual time: " << usual_time << std::endl;
  }

  void RandomFill(unsigned capacity) {
    Matrix<double> matr(capacity, capacity,
                        [&] { return Random::Normal<double>(-100, 100); });
    std::vector<double> consts(capacity);
    std::generate(consts.begin(), consts.end(),
                  [] { return Random::Uniform<double>(-100, 100); });
    gauss.set_equations(matr, consts);
  }

  void UsualSolve() {
    auto sol = gauss.UsualExecute();
    SStr::Print(sol);
  }

  void ParallelSolve(unsigned threads_count) {
    auto sol = gauss.ParallelExecute(threads_count);
    SStr::Print(sol);
  }

 private:
  SLAEGauss gauss;
};
}  // namespace s21
