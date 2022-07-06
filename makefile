CC = g++
CFLAGS = -std=c++11 -O3 -fopenmp -fpermissive -w

SRCS =./Src/*.cpp
PROG = ./ipsr

$(PROG) : $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) $(INCLUDE) -o $(PROG)
