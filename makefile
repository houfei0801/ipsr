#If use Clang, please change g++ to clang++
CC = g++
CFLAGS = -std=c++11 -O3 -fopenmp -w

SRCS =./Src/*.cpp
PROG = ./ipsr

$(PROG) : $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(PROG)
