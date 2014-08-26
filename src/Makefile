# Makefile for asd

#G++FLAG = -g
G++FLAG = -O3

asd : asd.o param_t.o gzstream.o
	g++ -o asd asd.o param_t.o gzstream.o -lm -lpthread -lz

asd.o : asd.cpp
	g++ $(G++FLAG) -c asd.cpp

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.C
	g++ $(G++FLAG) -c gzstream.C

clean:
	rm *.o
