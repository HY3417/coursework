all: compile

default: compile

compile : c1.cpp c2.cpp c3.cpp


	#compile c1
	g++ -std=c++11 -c c1.cpp
	#compile c2
	g++ -std=c++11 -c c2.cpp
	#compile c3
	g++ -std=c++11 -c c3.cpp
	#link c1 and main
	g++ -o cal1 c1.o -L/usr/local/lib -llapack -lblas
	#link c2 and main
	g++ -o cal2 c2.o -L/usr/local/lib -llapack -lblas
	#link c3 and main
	g++ -o cal3 c3.o -L/usr/local/lib -llapack -lblas

c1 : cal1
	./cal1 0 1 1 2 0.2 250 250 0 10 5

c2 : cal2
	./cal2 0 1 1 2 0.2 250 250 0 10 5

c3 : cal3
	./cal3 0.25 1 1.3 3 0.2 250 250 0 15 8

.PHONY: clean
clean:
	rm -f *.o cal1 cal2 cal3
