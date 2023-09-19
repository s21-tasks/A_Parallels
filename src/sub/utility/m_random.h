#pragma once

#include <iostream>
#include <map>
#include <random>
#include <thread>
#include <type_traits>

namespace s21 {

namespace Random {

template <class Type, class Distribution, class Engine, class... Args>
Type Full(Args &&...args) {
  static Engine generator{std::random_device{}()};
  static thread_local std::map<std::tuple<std::decay_t<Args>...>, Distribution>
      distributions;
  return distributions.emplace(std::tie(args...), Distribution(args...))
      .first->second(generator);
}

template <class FPType = float,
          class Distribution = std::normal_distribution<FPType>,
          class Engine = std::mt19937>
FPType Normal(FPType mean, FPType sd) {
  return Full<FPType, Distribution, Engine>(mean, sd);
}

template <class Distribution = std::bernoulli_distribution,
          class Engine = std::mt19937>
bool Bool(double chance) {
  return Full<bool, Distribution, Engine>(chance);
}

template <class FPType = float, class Engine = std::mt19937>
FPType Uniform(FPType min, FPType max = 1.0) {
  return Full<FPType, std::uniform_real_distribution<FPType>, Engine>(min, max);
}

template <class IntType = int, class Engine = std::mt19937>
IntType Int(IntType min, IntType max = std::numeric_limits<IntType>::max()) {
  return Full<IntType, std::uniform_int_distribution<IntType>, Engine>(min,
                                                                       max);
}

template <class Type, class = void>
struct Easy;

template <class Type>
struct Easy<Type, std::enable_if_t<std::is_integral<Type>::value &&
                                   !std::is_same<Type, bool>::value>> {
  static Type R(Type min, Type max) { return Int<Type>(min, max); }
};

template <typename Type>
struct Easy<Type, std::enable_if_t<std::is_floating_point<Type>::value>> {
  static Type R(Type min, Type max) { return Uniform<Type>(min, max); }
};

}  // namespace Random

}  // namespace s21
