
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
    PrePostSolver<N> solver;
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
    auto template_solver = MakePrePostSolver(template_stepper, DoNothingAlgebraic<N>{}, more_stupid_algebra);
    // Time
    std::cout << std::setw(W) << "Template Solver : " << timing_mean_std(BATCHES, REPEATS, template_solver, state, 0.1) << " seconds" << std::endl;

    auto reduced_solver = MakePrePostSolver(template_stepper, DoNothingAlgebraic<N>{}, DoNothingAlgebraic<N>{});
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
    const double dt = 0.001;
    const size_t steps = static_cast<size_t>(tf/dt);

    // Initial conditions
    State<N> state = {0.0, 10.0, 0.0};

    const auto bouncer = MakeWrappedAlgebraic<N>(
            [=](State<N>& s, const double&) -> void
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
    const auto alg_stpper = MakeCompositeAlgebraic<N>(Incrementor<N, t>{}, ode_stepper, bouncer);
    const auto solver = MakePrePostSolver(ode_stepper, Incrementor<N, t>{}, bouncer);

    auto out_file = std::ofstream("data/bouncy_ball.csv");
    out_file << "t, y, ydot" << std::endl;
    for (size_t i = 0; i < steps; ++i)
    {
        alg_stpper(state, dt);
        out_file << state << std::endl;
    }
}

void forced_damped_oscillator_test()
{
    // Constants
    enum variable : size_t { t, x, xdot };
    constexpr size_t N = 3;
    const double m = 1.0;
    const double k = 1.0;
    const double c = 0.1;
    const double F = 1.0;
    const double om = 1.0;

    // Parameters
    const double tf = 100.0;
    const double dt = 0.001;
    const size_t steps = static_cast<size_t>(tf/dt);

    // Initial conditions
    State<N> state = {0.0, 1.0, 0.0};

    // Derivative
    auto deriv = [=](const State<N>& s) -> State<N>
    {
        State<N> ds(0.0);
        ds[x] = s[xdot];
        ds[xdot] = (F*std::sin(om*s[t]) - c*s[xdot] - k*s[x])/m;
        return ds;
    };

    // Stepper
    auto stepper = MakeRK4Stepper<N>(deriv);
    auto solver = MakeCompositeAlgebraic<N>(Incrementor<N, t>{}, stepper);

    auto out_file = std::ofstream("data/forced_damped_oscillator.csv");
    for (size_t i = 0; i < steps; ++i)
    {
        solver(state, dt);
        out_file << state << std::endl;
    }
}

void heart_beat_test()
{
    // State
    enum variable : size_t
    {
        t,
        Qao,
        Psa
    };
    constexpr size_t N = 3;

    // Constants
    const double T = 0.0125;    // Duration of heart beat (minutes)
    const double TS = 0.0050;   // Duration of systole (minutes)
    const double TMAX = 0.002;  // Time at which flow is max (minutes)
    const double QMAX = 28.0;   // Max flow through aortic valve (liters/minute)
    const double RS = 17.86;    // Resistance of systemic circulation (mmHg/(liter/minute))
    const double CSA = 0.00175; // Capacitance of systemic arteries (liters/mmHg)

    const double Ttotal = 10.0 * T;
    const double dt = 0.01 * T;
    const size_t steps = static_cast<size_t>(Ttotal / dt);

    auto Qao_func = [=](const double &tin) -> double
    {
        const double t = tin - std::floor(tin / T) * T;
        if (t < TMAX)
            return QMAX * t / TMAX;
        else if (t < TS)
            return QMAX * (TS - t) / (TS - TMAX);
        return 0.0;
    };

    auto setter = [=](State<N>& s, const double& _dt) -> void
    {
        s[t] += _dt;
        s[Qao] = Qao_func(s[t]);
    };

    auto derivative = [=](const State<N> &s) -> State<N>
    {
        State<N> ds(0.0);
        ds[Psa] = (Qao_func(s[t]) - s[Psa] / RS) / CSA;
        return ds;
    };

    auto stepper = MakeRK4ExplicitTimeStepper<N>(derivative, setter);
    auto solver = MakeCompositeAlgebraic<N>(stepper);

    State<N> state = {0.0, 0.0, 80.0};
    auto out_file = std::ofstream("data/explicit_heart_beat.csv");
    for (size_t i = 0; i < steps; ++i)
    {
        solver(state, dt);
        out_file << state << std::endl;
    }
    out_file.close();
}

void heart_beat_two()
{
    // State
    enum variable : size_t
    {
        t,      // Set
        PLV,    // Solved
        Psa,    // Solved
        QMi,    // Set
        QAo,    // Set
        Qs,     // Set
        SMi,    // Set
        SAo     // Set
    };
    constexpr size_t N = SAo + 1;

    // Constants
    const double T = 0.0125;        // Duration of heart beat (minutes)
    const double TS = 0.0050;       // Duration of systole (minutes)
    const double RS = 17.86;        // Resistance of systemic circulation (mmHg/(liter/minute))
    const double CSA = 0.00175;     // Capacitance of systemic arteries (liters/mmHg)
    const double PLA = 5.0;         // Left atrial pressure (mmHg)
    const double RMi = 0.01;        // Mitral valve resistance (mmHg/(liter/minute))
    const double RAO = 0.01;        // Aortic valve resistance (mmHg/(liter/minute))
    const double CLVS = 0.00003;    // Compliance of left ventricle in systole (liters/mmHg)
    const double CLVD = 0.0146;     // Compliance of left ventricle in diastole (liters/mmHg)
    const double TAU_S = 0.0025;    // Time constant for ventricle pressure decay during systole (minutes)
    const double TAU_D = 0.0075;    // Time constant for ventricle pressure decay during diastole (minutes)

    // Run params
    const double Ttotal = 15.0 * T;
    const double dt = 0.01 * T;
    const size_t steps = static_cast<size_t>(Ttotal / dt);

    auto CLV_func = [=](const double &tin) -> double
    {
        const double t = tin - std::floor(tin / T) * T;
        if (t < TS)
            return CLVD * pow(CLVS/CLVD, (1.0-exp(-t/TAU_S)) / (1.0-exp(-TS/TAU_S)) );
        return CLVS * pow(CLVD/CLVS, (1.0-exp(-(t-TS)/TAU_D)) / (1.0-exp(-(T-TS)/TAU_D)) );
    };

    auto dCLV_dt = [=](const double& tin) -> double
    {
        const double t = tin - std::floor(tin / T) * T;
        if (t < TS)
            return (CLVD * pow(CLVS / CLVD, (1.0 - exp(-t / TAU_S)) * (1.0 + 1.0 / (-1.0 + pow(exp(1.0), TS / TAU_S)))) * log(CLVS / CLVD)) / (exp(t / TAU_S) * (TAU_S - TAU_S / pow(exp(1.0), TS / TAU_S)));
        return (pow(CLVD / CLVS, (1.0 - exp((-t + TS) / TAU_D)) / (1.0 - exp((-T + TS) / TAU_D))) * CLVS * exp((-t + TS) / TAU_D) * log(CLVD / CLVS)) / (TAU_D - exp((-T + TS) / TAU_D) * TAU_D);
    };

    auto derivative = [=](const State<N> &s) -> State<N>
    {
        State<N> ds(0.0);
        const double SMi = static_cast<double>(PLA <= s[PLV]);
        const double SAo = static_cast<double>(s[PLV] <= s[Psa]);
        const double CLV = CLV_func(s[t]);
   
        ds[t]   = 1.0;
        ds[PLV] = (s[SMi] * (s[PLA] - s[PLV]) / RMi - s[SAo] * (s[PLV] - s[Psa]) / RAO) / CLV;
        ds[Psa] = (s[SAo] * (s[PLV] - s[Psa]) / RAO - s[Psa] / RS) / CSA;
        return ds;
    };

    auto stepper = MakeRK4Stepper<N>(derivative);

    State<N> state = {0.0, 5.0, 80.0, 5.0, 5.0, 5.0, 0.0, 0.0};
    print(state);

    auto out_file = std::ofstream("data/heart_beat_two.csv");
    out_file << "# t, PLV, Psa, QMi, QAo, Qs, SMi, SAo" << std::endl;
    for (size_t i = 0; i < steps; ++i)
    {
        stepper(state, dt);
        out_file << state << std::endl;
        std::cout << state << std::endl;
    }
    out_file.close();
}

int main()
{
    //timing_tests();

    //bouncy_ball_test();

    //forced_damped_oscillator_test();

    //heart_beat_test();

    heart_beat_two();

    return 0;
}
