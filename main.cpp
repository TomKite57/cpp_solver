
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <valarray>
#include <string>
#include <functional>
#include <cassert>
#include <chrono>

#include "utils.hpp"
#include "state.hpp"
#include "ODE_stepper.hpp"
#include "algebraic.hpp"
#include "derivative.hpp"
#include "solver.hpp"


void timing_tests()
{
    enum variable : size_t { x, xdot, y, ydot };
    constexpr size_t N = 5;
    constexpr size_t W = 20;
    constexpr size_t REPEATS = 5'000'000;
    constexpr size_t BATCHES = 10;
    State<N> state = {1.0, 2.0, 3.0, 4.0, 0.0};

    const auto f1 = MakeWrappedDerivative<N>([](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; });
    const auto f2 = MakeWrappedDerivative<N>([](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; });

    auto TDeriv = MakeCompositeDerivative<N>(f1, f2);
    RK4Stepper<N> Tstepper(TDeriv);
    std::cout << std::setw(W) << "Template Stepper : " << timing_mean_std(BATCHES, REPEATS, Tstepper, state, 0.1) << " seconds" << std::endl;

    // Simple Stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    RK4Stepper<N> stepper
    {
        [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
        [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
    };
    // Time
    std::cout << std::setw(W) << "Composite Stepper : " << timing_mean_std(BATCHES, REPEATS, stepper, state, 0.1) << " seconds" << std::endl;

    // Template stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};

    //RK4_template_stepper<N, CompositeDerivative<N>> template_stepper
    auto template_stepper = MakeRK4Stepper<N>
    (
        CompositeDerivative<N>
        {
            [](const State<N>& s){ State<N> ds(0.0); ds[x] = s[xdot]; ds[xdot] = -s[x]; return ds; },
            [](const State<N>& s){ State<N> ds(0.0); ds[y] = s[ydot]; ds[ydot] = -2.0*s[y]; return ds; }
        }
    );

    std::cout << std::setw(W) << "Template Stepper : " << timing_mean_std(BATCHES, REPEATS, template_stepper, state, 0.1) << " seconds" << std::endl;

    // Dynamic Stepper
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    DynamicSolver<N> solver;
    solver.set_poststep( Incrementor<N, N-1>{});
    solver.set_ODE_step(stepper);
    // Time
    std::cout << std::setw(W) << "Dynamic Solver : " << timing_mean_std(BATCHES, REPEATS, solver, state, 0.1) << " seconds" << std::endl;

    // Template Solver
    state = {1.0, 2.0, 3.0, 4.0, 0.0};
    auto incrementor = MakeWrappedAlgebraic<N>([](State<N>& s, const double& dt) -> void { s[N-1] += dt; });
    auto decrementor = MakeWrappedAlgebraic<N>([](State<N>& s, const double& dt) -> void { s[N-1] -= dt; });

    auto full_func = MakeCompositeAlgebraic<N>(incrementor, template_stepper);

    auto more_stupid_algebra = MakeCompositeAlgebraic<N>(
        incrementor,
        decrementor,
        incrementor
    );
    auto template_solver = MakeSolver(template_stepper, DoNothingAlgebraic<N>{}, more_stupid_algebra);
    // Time
    std::cout << std::setw(W) << "Template Solver : " << timing_mean_std(BATCHES, REPEATS, template_solver, state, 0.1) << " seconds" << std::endl;

    auto reduced_solver = MakeSolver(template_stepper, DoNothingAlgebraic<N>{}, DoNothingAlgebraic<N>{});
    // Time
    std::cout << std::setw(W) << "Reduced Solver : " << timing_mean_std(BATCHES, REPEATS, reduced_solver, state, 0.1) << " seconds" << std::endl;
}

void bouncy_ball_test()
{
    // Constants
    enum variable : size_t { t, y, ydot };
    constexpr size_t N = 3;
    const double g = 9.81;
    const double e = 0.98;

    assert(e <= 1.0 && e >= 0.0);

    // Parameters
    const double tf = 100.0;
    const double dt = 0.01;
    const size_t steps = static_cast<size_t>(tf/dt);

    // Initial conditions
    State<N> state = {0.0, 10.0, 0.0};

    const auto bouncer = MakeWrappedAlgebraic<N>(
            [=](State<N>& s, const double& dt) -> void
            {
                if (s[y] < 0.0)
                {
                    s[y] = -s[y];
                    s[ydot] = -e*s[ydot];
                }
            }
        );

    const auto gravity = MakeWrappedDerivative<N>(
            [=](const State<N>& s) -> State<N>
            {
                State<N> ds(0.0);
                ds[y] = s[ydot];
                ds[ydot] = -g;
                return ds;
            }
        );

    const auto ode_stepper = MakeRK4Stepper<N>(gravity);
    const auto solver = MakeSolver(ode_stepper, Incrementor<N, t>{}, bouncer);

    auto out_file = std::ofstream("bouncy_ball.dat");
    out_file << "t, y, ydot" << std::endl;
    for (size_t i = 0; i < steps; ++i)
    {
        solver(state, dt);
        out_file << state << std::endl;
    }
}

int main()
{
    //timing_tests();

    bouncy_ball_test();

    return 0;
}
