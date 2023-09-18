#include <list>
#include <vector>

#include "m_sstr.h"

int main() {
  s21::SStr::Print(std::vector<std::pair<int, int>>{
      std::make_pair(1, 2), std::make_pair(1, 2), std::make_pair(1, 2)});
  s21::SStr::Print(std::vector<std::vector<int>>{{1, 2, 3}, {4, 5, 6}});
  s21::SStr::Print(std::vector<int>{1, 2, 3});

  s21::SStr::Print<' '>(1, "qwe", 4.4);

  s21::SStr::Print(std::pair<std::string, double>{"qwerty", -12.2},
                   std::list<std::string>{"1", "qwe", "ccc"}, 12, 12, 12);

  return 0;
}
