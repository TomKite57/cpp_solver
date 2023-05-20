
#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <functional>
#include <cassert>
#include <chrono>

#include "utils.hpp"
#include "state.hpp"
#include "ODE_stepper.hpp"
#include "function.hpp"
#include "solver.hpp"

enum variable : size_t { x, xdot, y, ydot };

int main()
{
    constexpr size_t N = 5;
    State<N> state = {1.0, 2.0, 3.0, 4.0, 0.0};

    auto f1 = CircularAcceleration<y, ydot, N>(2.0);
    auto f2 = CircularAcceleration<x, xdot, N>(1.0);

    CompositeDerivative<N> dfunc{};
    dfunc.add_function(f1);
    dfunc.add_function(f2);

    RK4_stepper<N> stepper{dfunc};

    // Time loop
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        stepper(state, 0.1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state << std::endl;
    std::cout << std::endl;


    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    DynamicSolver<N> solver;
    solver.set_poststep( Incrementor<N-1, N>{});
    solver.set_ODE_step(stepper);
    // Time loop
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        solver.step(state, 0.1);
    }
    end = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state << std::endl;
    std::cout << std::endl;


    return 0;
}
