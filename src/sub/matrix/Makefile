all: test #speed

speed:
	g++ speed_test.cc -o speed -lblas -lOpenCL -lclBLAS
	./speed

test:
	g++ blas_cl_test.cc -o test -lblas -lOpenCL -lclBLAS -lpthread -lgtest -lgtest_main
	./test

clean:
	rm -rf test speed