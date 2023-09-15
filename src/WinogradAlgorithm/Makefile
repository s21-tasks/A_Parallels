.PHONY: functional_test speed_compare_test level_load_test level_load_test_optimized speed_compare_test_optimized functional_basic_test functional_parallel_test clean all

CC=g++
FLAGS=-Wall -Werror -Wextra
OPTIMIZATION=-Ofast -march=native
STANDART=-std=c++17
BLAS_FLAGS=-lpthread -lm -ldl -lblas
GTEST_FLAGS=-lgtest -lgtest_main -pthread
SRC=test.cc
EXEC=winograd

all: clean
	$(CC) $(FLAGS) $(SRC) -o winograd $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION)
	./winograd

functional_basic_test:
	$(CC) $(FLAGS) -D BASIC_WINOGRAD tests/functional.cc -o functional_basic_test $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION) $(GTEST_FLAGS)
	./functional_basic_test

functional_parallel_test:
	$(CC) $(FLAGS) -D PARALLEL_WINOGRAD tests/functional.cc -o functional_parallel_test $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION) $(GTEST_FLAGS)
	./functional_parallel_test

functional_test: functional_basic_test functional_parallel_test

level_load_test:
	$(CC) $(FLAGS) tests/level_load.cc -o level_load_test $(STANDART)
	./level_load_test

speed_compare_test:
	$(CC) $(FLAGS) tests/speed_compare.cc -o speed_compare_test $(STANDART) $(BLAS_FLAGS)
	./speed_compare_test

level_load_test_optimized:
	$(CC) $(FLAGS) tests/level_load.cc -o level_load_test_optimized $(STANDART) $(OPTIMIZATION)
	./level_load_test_optimized

speed_compare_test_optimized:
	$(CC) $(FLAGS) tests/speed_compare.cc -o speed_compare_test_optimized $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION)
	./speed_compare_test_optimized

clean:
	rm -rf winograd
	rm -rf speed_compare_test
	rm -rf functional_basic_test
	rm -rf functional_parallel_test
	rm -rf level_load_test
	rm -rf level_load_test_optimized
	rm -rf speed_compare_test_optimized
