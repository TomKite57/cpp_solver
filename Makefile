all:
	g++ -std=c++17 -O3 -DNDEBUG main.cpp -o cpp_solver.o

debug:
	g++ -std=c++17 -O3 main.cpp -o cpp_solver.o

clean:
	rm -rf cpp_solver.o