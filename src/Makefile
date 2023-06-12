

all:

gauss_test:
	mkdir -p SystemsOfLinearEquationsSolve/build
	cmake -B"SystemsOfLinearEquationsSolve/build" SystemsOfLinearEquationsSolve/.
	make -C SystemsOfLinearEquationsSolve/build/.
	./SystemsOfLinearEquationsSolve/build/gauss


submodules:
	cd ..; git submodule update --init --recursive; git submodule update --recursive --remote

.PHONY: submodules



