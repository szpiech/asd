# Makefile for asd

#G++FLAG = -g

#osx
G++FLAG = -fast

#linux
#G++FLAG = -O3


asd : asd.o param_t.o gzstream.o asd-cli.o errlog.o asd-data.o asd-dist.o
	g++ -o asd asd.o param_t.o gzstream.o asd-cli.o errlog.o asd-data.o asd-dist.o -lm -lpthread -lz

asd.o : asd.cpp
	g++ $(G++FLAG) -c asd.cpp

asd-cli.o : asd-cli.cpp
	g++ $(G++FLAG) -c asd-cli.cpp

asd-data.o : asd-data.cpp
	g++ $(G++FLAG) -c asd-data.cpp

asd-dist.o : asd-dist.cpp
	g++ $(G++FLAG) -c asd-dist.cpp

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

errlog.o : errlog.cpp
	g++ $(G++FLAG) -c errlog.cpp

gzstream.o : gzstream.C
	g++ $(G++FLAG) -c gzstream.C

clean:
	rm *.o
