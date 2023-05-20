
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

    auto f1 = std::make_unique<CircularAcceleration<x, xdot, N>>(1.0);
    auto f2 = std::make_unique<CircularAcceleration<y, ydot, N>>(2.0);

    auto dfunc = std::make_unique<CompositeDerivative<N>>();
    dfunc->add_function(std::move(f1));
    dfunc->add_function(std::move(f2));

    //RK4_stepper<N> stepper{std::move(dfunc)};
    std::unique_ptr<ODE_stepper<N>> stepper = std::make_unique<RK4_stepper<N>>(std::move(dfunc));

    // Time loop
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        (*stepper)(state, 0.1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state << std::endl;
    std::cout << std::endl;


    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    DynamicSolver<N> solver{std::move(stepper)};
    solver.set_poststep( std::make_unique<Incrementor<N-1, N>>());
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
