#compiler and debugger
CC = g++
CFLAGS = -g -Wall

#Include directories (other header files) and libraries
#INCLUDES =
#LFLAGS =
#LIBS =

#first target entry
default: math

#create executable
math: matrix.o main.cpp
	$(CC) $(CFLAGS) -o math main.cpp matrix.o

#create object file(s)
matrix.o: matrix.h matrix.cpp
	$(CC) $(CFLAGS) -c matrix.cpp

#clean
clean:
	$(RM) math *.o *~
