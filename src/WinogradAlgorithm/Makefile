CC=g++
FLAGS=-Wall -Werror -Wextra
OPT_FLAGS=-Ofast -march=native
STANDART=-std=c++17
BLAS_FLAGS=-lpthread -lm -ldl -lblas
SRC=test.cc

all: clean
	$(CC) $(FLAGS) $(STANDART) $(SRC) -o test $(BLAS_FLAGS) $(OPT_FLAGS)
	./test

clean:
	rm -rf test
