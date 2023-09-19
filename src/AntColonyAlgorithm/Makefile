.PHONY: all interface clean

CC=g++
STANDART=-std=c++17
SRC=../sub/ConsoleInterface/ConsoleInterface.cc ../sub/thread/m_thread.cc interface.cc

all: interface

interface:
	$(CC) $(SRC) -o ant $(STANDART)
	./ant

clean:
	rm -rf ant
