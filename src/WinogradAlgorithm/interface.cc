#include "../sub/ConsoleInterface/ConsoleInterface.h"
#include "../sub/matrix/blas.h"
#include "../sub/utility/utility.h"
#include "winograd.h"
#include "winograd_parallel.h"
#include "winograd_pipeline.h"
// #include "tests/tests.h"

#include <functional>

using namespace s21;

// class TwoDouble : public Input {
//  public:
//   TwoDouble(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
//   void Func() override;
// };

class Interface : public ConsoleInterface {
 public:
  Interface();

 private:
  using type = double;
  using M = Matrix<type>;
  M A;
  M B;

  struct TwoDouble : public Input {
    TwoDouble(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };
  struct MatrixInput : public Input {
    MatrixInput(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };
  struct InputSize : public OneValueInput<unsigned int> {
    InputSize(Interface *ci, AbsMenu *next_menu);
  };
  template <class Alg>
  void Multiply() {
    M C(A.GetCols());
    Alg::Mul(A, B, C);
    std::cout << "\nResult Matrix:\n" << C << "\n\n";
  }
  void TimeTest();
};

Interface::InputSize::InputSize(Interface *ci, AbsMenu *next_menu)
    : OneValueInput(
          ci, next_menu,
          [&](unsigned int &size) {
            Interface *i = reinterpret_cast<Interface *>(ci_);
            i->A = M(size);
            i->B = M(size);
          },
          "matrix size") {}

void Interface::TwoDouble::Func() {
  Style::InputRequest({"mean", "standard deviation"});
  auto a = Read<double>();
  auto b = Read<double>();
  if (Allowed()) {
    Interface *i = reinterpret_cast<Interface *>(ci_);
    i->A.Fill([a, b] { return Random::Normal(a, b); });
    i->B.Fill([a, b] { return Random::Normal(a, b); });
  }
}

void Interface::MatrixInput::Func() {
  Interface *i = reinterpret_cast<Interface *>(ci_);
  Style::InputRequest({"matrix A"});
  i->A.Fill([&] { return Read<type>(); });
  if (Allowed()) {
    Style::InputRequest({"matrix B"});
    i->B.Fill([&] { return Read<type>(); });
  }
}

void Interface::TimeTest() {}

Interface::Interface() {
  auto home_menu = AddMenu({"Random Matrices", "Input Matrices from console",
                            "Print Matrices", "Multiply", "Time Test"});

  home_menu->Connect(
      0, AddInput(new InputSize(this, new TwoDouble(this, home_menu))));

  home_menu->Connect(
      1, AddInput(new InputSize(this, new MatrixInput(this, home_menu))));

  home_menu->Connect(2, home_menu, [&] {
    std::cout << "\nMatrix A:\n" << A << "\n\nMatrix B:\n" << B << "\n\n";
  });

  auto mul_menu = AddMenu({"Basic Winograd", "Parallel Winograd",
                           "Pipeline Winograd", "Standart", "All"});

  home_menu->Connect(3, mul_menu);

  mul_menu->Connect(0, home_menu, [&] { Multiply<Basic::Winograd<type>>(); });
  mul_menu->Connect(1, home_menu,
                    [&] { Multiply<Parallel::Winograd<type>>(); });
  mul_menu->Connect(2, home_menu, [&] { Multiply<M>(); });
  mul_menu->Connect(3, home_menu, [&] { Multiply<M>(); });
  mul_menu->Connect(4, home_menu, [&] {
    std::cout << "\nBasic Winograd:\n";
    Multiply<Basic::Winograd<type>>();
    std::cout << "\nParallel Winograd:\n";
    Multiply<Parallel::Winograd<type>>();
    std::cout << "\nPipeline Winograd:\n";
    Multiply<PipeLine::Winograd<type>>();
    std::cout << "\nStandart:\n";
    Multiply<M>();
    std::cout << '\n';
  });

  auto time_menu =
      AddMenu({"Basic Winograd", "Parallel Winograd", "Pipeline Winograd",
               "Standart", "Standart Parallel", "Blas", "All"});

  home_menu->Connect(4, new OneValueInput<unsigned int>(
                            this, home_menu,
                            [&](unsigned int repeat) {

                            },
                            "repeats"));
}

int main() {
  Interface i;
  i.Start();

  // Matrix<int> A(3);
  // std::cin >> A;
  // A.PrintFull();

  return 0;
}
