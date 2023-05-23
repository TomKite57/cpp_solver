all:
	g++ -std=c++20 -O3 -DNDEBUG main.cpp -o cpp_solver.o

debug:
	g++ -std=c++20 -O3 main.cpp -o cpp_solver.o

warning:
	g++ -std=c++20 -O3 -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Werror -Wfatal-errors -DNDEBUG main.cpp -o cpp_solver.o

stupid:
	g++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror -Wfatal-errors -Wcast-align -Wcast-qual -Wconversion -Wdouble-promotion -Weffc++ -Wformat=2 -Winit-self -Winline -Wlogical-op -Wmissing-include-dirs -Wnon-virtual-dtor -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-conversion -Wstrict-null-sentinel -Wstrict-overflow=5 -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=pure -Wswitch-default -Wswitch-enum -Wundef -Wuninitialized -Wuseless-cast -Wvla -Wzero-as-null-pointer-constant main.cpp -o cpp_solver.o

clean:
	rm -rf cpp_solver.o