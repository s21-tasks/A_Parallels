#include "SLAEGauss.h"
#include "matrix.h"

int main() {
  s21::Matrix<double> coefficients({{2, 1, -1},
                                    {1, -1, 3},
                                    {3, -2, 1}});
  std::vector<double> constants = {1,
                                   -6,
                                   7};
  s21::SLAEGauss g(coefficients, constants);
  auto s = g.UsualExecute();
  for (int i = 0; i < s.size(); ++i) {
    std::cout << s[i] << " ";
  }
  std::cout << std::endl;


  return 0;
}