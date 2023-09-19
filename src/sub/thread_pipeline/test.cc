#include <iostream>

#include "thread_pipeline.h"

void a_func(int& a) { a += 1; }

void work(int&) { std::this_thread::sleep_for(std::chrono::seconds(1)); }

void print(int& a) { std::cout << a << std::endl; }

int main() {
  // создаем экземпляр
  Pipeline<int> pipe{};

  //добавляем 1 функцию в конвейер
  std::function<void(int&)> func = a_func;
  pipe.AddStage(a_func);

  std::function<void(int&)> func3 = work;
  pipe.AddStage(func3);

  //добавляем 2 функцию в конвейер
  std::function<void(int&)> func2 = print;
  pipe.AddStage(print);

  //создаем переменную и запускаем конвейер
  int a = 89;
  pipe.Process(a);

  //  std::this_thread::sleep_for(std::chrono::seconds(2));

  return 0;
}
