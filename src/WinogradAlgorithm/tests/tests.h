#pragma once

#include "../winograd_parallel.h"
#include "../winograd_pipeline.h"
#include "../../sub/utility/utility.h"

namespace s21 {

template<class T>
T RS() {
    return Random::Easy<T>::R(static_cast<T>(-10), static_cast<T>(10));
}

} // namespace s21
