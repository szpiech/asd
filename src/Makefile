# Makefile for asd

#G++FLAG = -g
G++FLAG = -O3

asd : asd.o param_t.o
	g++ -o asd asd.o param_t.o -lm -lpthread

asd.o : asd.cpp
	g++ $(G++FLAG) -c asd.cpp

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

clean:
	rm *.o
