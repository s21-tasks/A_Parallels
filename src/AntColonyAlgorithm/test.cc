#include "AntColony.h"
#include "AntColony_2.0.h"

using namespace s21;

typedef double fp_type;

int main() {
  const int size = 16;
  const float zero_probability = 0.0;
  int max_weight = 100;

  Matrix<fp_type> M1(size, size, [&](int k, int g, double &cell) {
    cell = (k != g && !Random::Bool(zero_probability))
               ? Random::Int(1, max_weight)
               : 0;
  });

  std::cout << M1;

  AntColony<fp_type> ac(M1);

  AntColony20<fp_type> ac20(M1);

  auto res = ac.Solve();
  std::cout << '\n' << res.distance << '\n';
  auto res20 = ac20.Solve();
  std::cout << '\n' << res20.distance << "\n\n";

  return 0;
}
