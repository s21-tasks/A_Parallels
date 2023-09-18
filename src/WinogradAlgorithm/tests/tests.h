#pragma once

#include "../../sub/utility/utility.h"
#include "../winograd_parallel.h"
#include "../winograd_pipeline.h"

namespace s21 {

template <class T>
T RS() {
  return Random::Easy<T>::R(static_cast<T>(-10), static_cast<T>(10));
}

}  // namespace s21
