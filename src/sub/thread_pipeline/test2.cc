#include <iostream>

#include "thread_pipeline.h"

void a_func(int& a) { a += 1; }

void load(int& a) { std::this_thread::sleep_for(std::chrono::seconds(1)); }

int main() {
  // создаем экземпляр
  Pipeline<int> pipe{};

  //добавляем 1 функцию в конвейер
  std::function<void(int&)> func = a_func;
  pipe.AddStage(func);

  //добавляем 2 функцию в конвейер
  std::function<void(int&)> func2 = load;
  pipe.AddStage(func2);


  pipe.AddStage([] (int &a) {
    std::cout << a << '\n';
  });

  //создаем переменную и запускаем конвейер
  int a = 89;
  int b = 12;
  int c = 13;
  pipe.Process(a);
  pipe.Process(b);
  pipe.Process(c);

  //сон чтобы главный поток не завершился раньше времени
  std::this_thread::sleep_for(std::chrono::seconds(1));

  return 0;
}
