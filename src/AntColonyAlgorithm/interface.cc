#include "../sub/ConsoleInterface/ConsoleInterface.h"
#include "AntColony_2.0.h"

using namespace s21;

class Interface : public ConsoleInterface {
  using type = unsigned int;

 public:
  Interface();

 private:
  Matrix<type> graph_;

  struct RandomGraph : public Input {
    RandomGraph(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };

  struct MatrixInput : public Input {
    MatrixInput(Interface *ci, AbsMenu *next_menu) : Input(ci, next_menu) {}
    void Func() override;
  };
};

void Interface::RandomGraph::Func() {
  Style::InputRequest({"size", "min value", "max value", "zero probability"});
  unsigned int size = Read<int>();
  auto min = Read<int>();
  auto max = Read<int>();
  auto zero = Read<double>();
  if (Allowed()) {
    if (size <= 1) throw std::invalid_argument("size must be > 1");
    if (min >= max)
      throw std::invalid_argument("min value must be < max value");
    if (min <= 0) throw std::invalid_argument("min value must be > 0");
    if (zero < 0.0 || zero > 1.0)
      throw std::invalid_argument("zero probability must be in [0.0, 1.0]");
    Interface *i = reinterpret_cast<Interface *>(ci_);
    i->graph_ = Matrix<type>(
        size, size,
        [min, max, zero](unsigned int k, unsigned int g, type &cell) {
          cell = (k != g && !Random::Bool(zero)) ? Random::Int(min, max) : 0;
        });
  }
}

void Interface::MatrixInput::Func() {
  Interface *i = reinterpret_cast<Interface *>(ci_);
  Style::InputRequest({"matrix A"});
  i->graph_.Fill([&] { return Read<type>(); });
}

Interface::Interface() {
  auto home_menu =
      AddMenu({"Random Graph", "Input Graph from console", "Print Graph",
               "Run Ant Colony", "Run Ant Colony 2.0", "Time Compare Test"});

  home_menu->Connect(0, AddInput(new RandomGraph(this, home_menu)));
  // home_menu->Connect(1, AddInput(new InputSize(this, new MatrixInput(this,
  // home_menu))));
  home_menu->Connect(
      1, new OneValueInput<unsigned int>(
             this, new MatrixInput(this, home_menu),
             [&](unsigned int &size) { graph_ = Matrix<type>(size); },
             "graph size"));
  home_menu->Connect(2, home_menu, [&] { std::cout << graph_ << '\n'; });
  home_menu->Connect(3, home_menu, [&] {
    AntColony<type> ac(graph_);
    TsmResult res = ac.Solve();
    SStr::Print(res.distance, res.vertices);
  });
  home_menu->Connect(4, home_menu, [&] {
    AntColony20<type> ac(graph_);
    TsmResult res = ac.Solve();
    SStr::Print(res.distance, res.vertices);
  });
  home_menu->Connect(
      5, new OneValueInput<unsigned int>(
             this, home_menu,
             [&](unsigned int &count) {
               auto result = Time::Compare(
                   count,
                   [&] {
                     AntColony<type> ac(graph_);
                     ac.Solve();
                   },
                   [&] {
                     AntColony20<type> ac20(graph_);
                     ac20.Solve();
                   });
               std::cout << "\tClassic Ant Colony: " << result[0] << '\n';
               std::cout << "\tAnt Colony 2.0: " << result[1] << '\n';
             },
             "number of iterations"));
}

int main() {
  Interface i;
  i.Start();
  return 0;
}
