#include "../sub/ConsoleInterface/ConsoleInterface.h"
#include "../sub/matrix/blas.h"
#include "../sub/utility/utility.h"
#include "winograd.h"
#include "winograd_parallel.h"
#include "winograd_pipeline.h"
// #include "tests/tests.h"

#include <functional>

using namespace s21;

class Interface : public ConsoleInterface {
 public:
  Interface();

 private:
  using type = double;
  using M = Matrix<type>;
  M A;
  M B;

  unsigned int N_;
  unsigned int repeat_ = 1;

  struct TwoDouble : public Input {
    TwoDouble(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };
  struct MatrixInput : public Input {
    MatrixInput(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };
  struct InputSize : public OneValueInput<int> {
    InputSize(Interface *ci, AbsMenu *next_menu);
  };
  template <class Alg>
  void Multiply(const std::string &name);

  template <class Alg>
  void TimeTest(const std::string &name);

  void TimeTest(const std::function<void(const M &, const M &, M &)> &func,
                const std::string &name);

  void TimeTestPipeline();
};

template <class Alg>
void Interface::Multiply(const std::string &name) {
  M C(N_);
  Alg::Mul(A, B, C);
  std::cout << '\n' << name << " Result Matrix:\n" << C << "\n\n";
}

void Interface::TimeTestPipeline() {
  auto time_point = std::chrono::high_resolution_clock::now();
  if (repeat_ % 4 == 0 && (N_ <= 8 || N_ == 16)) {
    std::vector<M> tasks(repeat_ / 4, A);
    std::vector<M> workers(4, B);
    time_point = std::chrono::high_resolution_clock::now();
    PipeLine::Winograd<type> W(workers, N_);
    for (auto &t : tasks) {
      W.Execute(t);
    }
  } else {
    M C(N_);
    time_point = std::chrono::high_resolution_clock::now();
    for (unsigned int k = 0; k < repeat_; ++k) {
      PipeLine::Winograd<type>::Mul(A, B, C);
    }
  }
  std::cout << "Pipeline time = "
            << Time::Duration<Time::mcs>(time_point) / repeat_ << " mcs\n";
}

template <class Alg>
void Interface::TimeTest(const std::string &name) {
  M C(N_);
  auto time_point = std::chrono::high_resolution_clock::now();
  Alg W(N_);
  for (unsigned int k = 0; k < repeat_; ++k) {
    W.Execute(A, B, C);
  }
  std::cout << name
            << " time = " << Time::Duration<Time::mcs>(time_point) / repeat_
            << " mcs\n";
}

void Interface::TimeTest(
    const std::function<void(const M &, const M &, M &)> &func,
    const std::string &name) {
  M C(N_);
  auto time_point = std::chrono::high_resolution_clock::now();
  for (unsigned int k = 0; k < repeat_; ++k) {
    func(A, B, C);
  }
  std::cout << name
            << " time = " << Time::Duration<Time::mcs>(time_point) / repeat_
            << " mcs\n";
}

Interface::InputSize::InputSize(Interface *ci, AbsMenu *next_menu)
    : OneValueInput(
          ci, next_menu,
          [&](int &size) {
            Interface *i = reinterpret_cast<Interface *>(ci_);
            if (size <= 0)
              throw std::invalid_argument("matrix size must be > 0");
            i->A = M(size);
            i->B = M(size);
            i->N_ = size;
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

  mul_menu->Connect(0, home_menu,
                    [&] { Multiply<Basic::Winograd<type>>("Basic"); });
  mul_menu->Connect(1, home_menu,
                    [&] { Multiply<Parallel::Winograd<type>>("Parallel"); });
  mul_menu->Connect(2, home_menu,
                    [&] { Multiply<PipeLine::Winograd<type>>("Pipeline"); });
  mul_menu->Connect(3, home_menu, [&] { Multiply<M>("Standart"); });
  mul_menu->Connect(4, home_menu, [&] {
    Multiply<Basic::Winograd<type>>("Basic");
    Multiply<Parallel::Winograd<type>>("Parallel");
    Multiply<PipeLine::Winograd<type>>("Pipeline");
    Multiply<M>("Standart");
  });

  auto spped_menu =
      AddMenu({"Basic Winograd", "Parallel Winograd", "Pipeline Winograd",
               "Standart", "Standart Parallel", "Blas", "All"});

  home_menu->Connect(4, AddInput(new OneValueInput<int>(
                            this, spped_menu,
                            [&](int r) {
                              if (r <= 0)
                                throw std::invalid_argument(
                                    "number of iterations must be > 0");
                              repeat_ = r;
                            },
                            "number of iterations")));

  auto standart = [](const M &A, const M &B, M &C) { M::Mul(A, B, C); };
  auto standart_par = [](const M &A, const M &B, M &C) {
    M::MulParallel(A, B, C);
  };
  auto blas = [](const M &A, const M &B, M &C) { BLAS<type>::Mul(A, B, C); };

  spped_menu->Connect(0, spped_menu,
                      [&] { TimeTest<Basic::Winograd<type>>("Basic"); });
  spped_menu->Connect(1, spped_menu,
                      [&] { TimeTest<Parallel::Winograd<type>>("Parallel"); });
  spped_menu->Connect(2, spped_menu, [&] { TimeTestPipeline(); });
  spped_menu->Connect(3, spped_menu, [&] { TimeTest(standart, "Standart"); });
  spped_menu->Connect(4, spped_menu,
                      [&] { TimeTest(standart_par, "Standart Parallel"); });
  spped_menu->Connect(5, spped_menu, [&] { TimeTest(blas, "Blas"); });

  spped_menu->Connect(6, home_menu, [&] {
    TimeTest<Basic::Winograd<type>>("Basic");
    TimeTest<Parallel::Winograd<type>>("Parallel");
    TimeTestPipeline();
    TimeTest(standart, "Standart");
    TimeTest(standart_par, "Standart Parallel");
    TimeTest(blas, "Blas");
  });
}

int main() {
  Interface i;
  i.Start();
  return 0;
}
