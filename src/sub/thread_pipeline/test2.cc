#include <iostream>

#include "../utility/utility.h"
#include "thread_pipeline.h"

typedef std::vector<int> object;

void f1(object &a) {
  object b(a.size(), 5);
  for (int k = 0; k < a.size(); ++k) {
    a[k] = a[k] + b[k];
  }
}

void f2(object &a) {
  object b(a.size(), -5);
  for (int k = 0; k < a.size(); ++k) {
    a[k] = a[k] + b[k];
  }
}

void f3(object &a) {
  object b(a.size(), 10);
  for (int k = 0; k < a.size(); ++k) {
    a[k] = a[k] + b[k];
  }
}

void f4(object &a) {
  object b(a.size(), -10);
  for (int k = 0; k < a.size(); ++k) {
    a[k] = a[k] + b[k];
  }
}

int main() {
  std::vector<object> objects(10000, object(100000, 1));
  {
    // создаем экземпляр
    Pipeline<object> pipe{};

    //добавляем 1 функцию в конвейер
    std::function<void(object &)> func = f1;
    pipe.AddStage(func);

    //добавляем 2 функцию в конвейер
    std::function<void(object &)> func2 = f2;
    pipe.AddStage(func2);

    //добавляем 3 функцию в конвейер
    std::function<void(object &)> func3 = f3;
    pipe.AddStage(func3);

    //добавляем 4 функцию в конвейер
    std::function<void(object &)> func4 = f4;
    pipe.AddStage(func4);

    std::cout << "start " << '\n';

    //создаем переменную и запускаем конвейер
    for (auto &obj : objects) {
      pipe.Process(obj);
    }
    std::cout << "end 1 " << '\n';
  }
  std::cout << "end 2 " << '\n';

  //сон чтобы главный поток не завершился раньше времени
  std::this_thread::sleep_for(std::chrono::seconds(1));

  return 0;
}
