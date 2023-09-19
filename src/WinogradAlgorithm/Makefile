.PHONY: functional_test speed_compare_test valgrind level_load_test \
		level_load_test_slow speed_compare_test_slow \
		functional_basic_test functional_parallel_test clean all \
		winograd pipeline_speed_compare pipeline_speed_compare_slow

CC=g++
FLAGS=-Wall -Werror -Wextra
OPTIMIZATION=-Ofast -march=native
STANDART=-std=c++17
BLAS_FLAGS=-lpthread -lm -ldl -lblas
GTEST_FLAGS=-lgtest -lgtest_main -pthread
INTERFACE_SRC=../sub/ConsoleInterface/ConsoleInterface.cc
THREAD_SRC=../sub/thread/m_thread.cc

all: winograd

winograd:
	$(CC) $(FLAGS) $(SRC) $(INTERFACE_SRC) interface.cc $(THREAD_SRC) -o winograd $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION)
	./winograd

fbt_build_helper:
	$(CC) $(FLAGS) -D BASIC_WINOGRAD tests/functional.cc -o functional_basic_test $(STANDART) $(GTEST_FLAGS)

functional_basic_test: fbt_build_helper
	./functional_basic_test

fpt_build_helper:
	$(CC) $(FLAGS) -D PARALLEL_WINOGRAD tests/functional.cc -o functional_parallel_test $(STANDART) $(GTEST_FLAGS)

functional_parallel_test: fpt_build_helper
	./functional_parallel_test

functional_test: functional_basic_test functional_parallel_test

speed_compare_test:
	$(CC) $(THREAD_SRC) $(FLAGS) tests/speed_compare.cc -o speed_compare_test $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION)
	./speed_compare_test

speed_compare_test_slow:
	$(CC) $(THREAD_SRC) $(FLAGS) tests/speed_compare.cc -o speed_compare_test_optimized $(STANDART) $(BLAS_FLAGS)
	./speed_compare_test_slow

level_load_test:
	$(CC) $(FLAGS) tests/level_load.cc -o level_load_test $(STANDART) $(OPTIMIZATION)
	./level_load_test

level_load_test_slow:
	$(CC) $(FLAGS) tests/level_load.cc -o level_load_test_slow $(STANDART)
	./level_load_test_slow

pipeline_speed_compare:
	$(CC) $(FLAGS) tests/pipeline.cc $(THREAD_SRC) -o pipeline_speed_compare $(STANDART) $(BLAS_FLAGS) $(OPTIMIZATION)
	./pipeline_speed_compare

pipeline_speed_compare_slow:
	$(CC) $(FLAGS) tests/pipeline.cc $(THREAD_SRC) -o pipeline_speed_compare_slow $(STANDART) $(BLAS_FLAGS)
	./pipeline_speed_compare_slow

valgrind: fbt_build_helper fpt_build_helper
	valgrind ./functional_basic_test
	valgrind ./functional_parallel_test

clean:
	rm -rf winograd
	rm -rf speed_compare_test
	rm -rf functional_basic_test
	rm -rf functional_parallel_test
	rm -rf level_load_test
	rm -rf level_load_test_slow
	rm -rf speed_compare_test_slow
	rm -rf pipeline_speed_compare
	rm -rf pipeline_speed_compare_slow
