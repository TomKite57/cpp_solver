
#include <iostream>
#include <iomanip>
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
    constexpr size_t W = 20;
    constexpr size_t R = 50'000'000;
    State<N> state = {1.0, 2.0, 3.0, 4.0, 0.0};

    const auto f1 = MakeWrappedDerivative<N>([](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; });
    const auto f2 = MakeWrappedDerivative<N>([](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; });

    auto TDeriv = MakeCompositeDerivative<N>(f1, f2);
    RK4_stepper<N> Tstepper(TDeriv);
    std::cout << std::setw(W) << "Template Stepper : " << timing_mean_std([&](){return time_function(Tstepper, R/10, state, 0.1);}, 10) << " seconds" << std::endl;

    // Simple Stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    RK4_stepper<N> stepper
    {
        [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
        [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
    };
    // Time
    std::cout << std::setw(W) << "Composite Stepper : " << timing_mean_std([&](){return time_function(stepper, R/10, state, 0.1);}, 10) << " seconds" << std::endl;

    // Template stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    RK4_template_stepper template_stepper
    {
        CompositeDerivative<N>
        {
            [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
            [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
        }
    };
    std::cout << std::setw(W) << "Template Stepper : " << timing_mean_std([&](){return time_function(template_stepper, R/10, state, 0.1);}, 10) << " seconds" << std::endl;

    // Dynamic Stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    DynamicSolver<N> solver;
    solver.set_poststep( Incrementor<N-1, N>{});
    solver.set_ODE_step(stepper);
    // Time
    std::cout << std::setw(W) << "Dynamic Solver : " << timing_mean_std([&](){return time_function(solver, R/10, state, 0.1);}, 10) << " seconds" << std::endl;

    // Template Solver
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    auto template_solver = MakeTemplateSolver(template_stepper, DoNothingAlgebraic<N>{}, CompositeAlgebraic<N>{Incrementor<N-1, N>{}, Incrementor<N-1, N>{}, Decrementor<N-1, N>{}});
    // Time
    std::cout << std::setw(W) << "Template Solver : " << timing_mean_std([&](){return time_function(template_solver, R/10, state, 0.1);}, 10) << " seconds" << std::endl;

    return 0;
}
