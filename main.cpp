
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

    RK4_stepper<N> stepper
    {
        [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
        [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
    };

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


    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    ConstSolver<N> const_solver{stepper, DoNothing<N>{}, Incrementor<N-1, N>{}};
    // Time loop
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        const_solver.step(state, 0.1);
    }
    end = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state << std::endl;
    std::cout << std::endl;


    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    auto temp = CompositeAlgebraic<N>{Incrementor<N-1, N>{}, Incrementor<N-1, N>{}, Decrementor<N-1, N>{}};
    //temp.add_function( Incrementor<N-1, N>{}, Incrementor<N-1, N>{}, Decrementor<N-1, N>{} );
    ConstSolver<N> new_const_solver{stepper, DoNothing<N>{}, temp};
    // Time loop
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        new_const_solver.step(state, 0.1);
    }
    end = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state << std::endl;
    std::cout << std::endl;


    return 0;
}
