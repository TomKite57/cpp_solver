
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

    // Simple Stepper
    RK4_stepper<N> stepper
    {
        [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
        [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
    };
    // Time
    std::cout << time_function([&](){stepper(state, 0.1);}, 5'000'000) << " seconds" << std::endl;


    // Dynamic Stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    DynamicSolver<N> solver;
    solver.set_poststep( Incrementor<N-1, N>{});
    solver.set_ODE_step(stepper);
    // Time
    std::cout << time_function([&](){solver.step(state, 0.1);}, 5'000'000) << " seconds" << std::endl;


    // Const Stepper with Composite Algebraic
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    ConstSolver<N> const_solver{stepper, DoNothing<N>{}, Incrementor<N-1, N>{}};
    // Time
    std::cout << time_function([&](){const_solver.step(state, 0.1);}, 5'000'000) << " seconds" << std::endl;


    // Const Stepper with Composite Algebraic
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    auto temp = CompositeAlgebraic<N>{Incrementor<N-1, N>{}, Incrementor<N-1, N>{}, Decrementor<N-1, N>{}};
    ConstSolver<N> new_const_solver{stepper, DoNothing<N>{}, temp};
    // Time
    std::cout << time_function([&](){new_const_solver.step(state, 0.1);}, 5'000'000) << " seconds" << std::endl;


    return 0;
}
