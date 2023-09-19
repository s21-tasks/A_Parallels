SRC := $(shell find . -name "*.cc")
HEADERS := $(shell find . -name "*.h")

.PHONY: all clean ant clean_ant winograd clean_winograd gauss gauss_run clean_gauss lint get_lint echo

all: ant gauss winograd

clean: clean_gauss clean_ant clean_winograd

ant:
	@$(MAKE) -C AntColonyAlgorithm

clean_ant:
	@$(MAKE) clean -C AntColonyAlgorithm

winograd:
	@$(MAKE) -C WinogradAlgorithm

clean_winograd:
	rm -rf WinogradAlgorithm/winograd

gauss:
	mkdir -p SystemsOfLinearEquationsSolve/build
	cmake -B"SystemsOfLinearEquationsSolve/build" SystemsOfLinearEquationsSolve/.
	make -C SystemsOfLinearEquationsSolve/build/.

gauss_run: gauss
	./SystemsOfLinearEquationsSolve/build/gauss

clean_gauss:
	rm -rf SystemsOfLinearEquationsSolve/build

lint:
	@clang-format --style=Google --Werror --dry-run $(SRC) $(HEADERS) && echo "Clang-format passed"
	@cppcheck $(SRC) $(HEADERS) --language=c++

get_lint:
	@clang-format --style=Google -i $(SRC) $(HEADERS) && echo "Formatted"

echo:
	@echo "HEADERS = $(HEADERS)"
	@echo "SRC = $(SRC)"
