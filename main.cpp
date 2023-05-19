
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

enum variable : size_t { x, xdot, y, ydot };

constexpr size_t ind(variable var) {
    return static_cast<size_t>(var);
}

int main()
{
    std::cout << x << std::endl;
    std::cout << xdot << std::endl;
    std::cout << y << std::endl;
    std::cout << ydot << std::endl;

    State<4> state = {1.0, 2.0, 3.0, 4.0};

    std::shared_ptr<Derivative<4>> f1 = std::make_shared<CircularAcceleration<x, xdot, 4>>(1.0);
    std::shared_ptr<Derivative<4>> f2 = std::make_shared<CircularAcceleration<y, ydot, 4>>(2.0);

    CompositeDerivative<4> dfunc;
    dfunc.add_function(f1);
    dfunc.add_function(f2);

    RK4_stepper<4> stepper{dfunc};

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

    return 0;
}
