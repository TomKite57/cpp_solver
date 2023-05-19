
#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <functional>
#include <cassert>
#include <chrono>

#include "utils.hpp"
#include "state.hpp"

template<size_t N>
class DerivativeFunction
{
public:
    DerivativeFunction() = default;
    DerivativeFunction(const DerivativeFunction& d) = default;
    DerivativeFunction(DerivativeFunction&& d) = default;
    DerivativeFunction& operator=(const DerivativeFunction& d) = default;
    DerivativeFunction& operator=(DerivativeFunction&& d) = default;
    virtual ~DerivativeFunction() = default;

    State<N> operator()(const State<N>& state) const
    {
        return call(state);
    }

    virtual State<N> call(const State<N>& state) const = 0;
};

template<size_t x_ind, size_t dx_ind, size_t N>
class CircularAcceleration : public DerivativeFunction<N>
{
public:
    double _omega2;
    static constexpr size_t _x_ind = x_ind;
    static constexpr size_t _dx_ind = dx_ind;

public:
    CircularAcceleration(double omega2) : _omega2(omega2) {}

    CircularAcceleration(const CircularAcceleration& d) = default;
    CircularAcceleration(CircularAcceleration&& d) = default;
    CircularAcceleration& operator=(const CircularAcceleration& d) = default;
    CircularAcceleration& operator=(CircularAcceleration&& d) = default;
    virtual ~CircularAcceleration() = default;

    State<N> call(const State<N>& state) const override
    {
        assert(state.size() > x_ind);
        assert(state.size() > dx_ind);

        State<N> dstate(0.0, state.size());
        dstate[x_ind] = state[dx_ind];
        dstate[dx_ind] = -_omega2*state[x_ind];
        return dstate;
    }
};

template<size_t N>
class ODE_stepper
{
public:
    void operator()(State<N>& state, const double& dt) const
    {
        step(state, dt);
    }
    virtual void step(State<N>& state, const double& dt) const = 0;
};

template<size_t N>
class Euler_stepper : public ODE_stepper<N>
{
private:
    const std::function<State<N>(const State<N>&)> _dfunc;
public:
    Euler_stepper(const std::function<State<N>(const State<N>&)> dfunc) : _dfunc{dfunc} {}

    void inline step(State<N>& state, const double& dt) const override
    {
        State<N> dstate = _dfunc(state);

        assert(state.size() == dstate.size());

        state += dstate*dt;
    }
};

template<size_t N>
class RK4_stepper : public ODE_stepper<N>
{
private:
    const std::function<State<N>(const State<N>&)> _dfunc;

public:
    RK4_stepper(const std::function<State<N>(const State<N>&)>& dfunc) : _dfunc{dfunc} {}

    void inline step(State<N>& state, const double& dt) const override
    {
        State<N> k1 = _dfunc(state);
        State<N> k2 = _dfunc(state + k1*dt/2.0);
        State<N> k3 = _dfunc(state + k2*dt/2.0);
        State<N> k4 = _dfunc(state + k3*dt);

        assert(state.size() == k1.size());
        assert(state.size() == k2.size());
        assert(state.size() == k3.size());
        assert(state.size() == k4.size());

        state += (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
    }
};


int main()
{
    State<2> state = {0.0, 1.0};

    CircularAcceleration<0, 1, 2> dfunc = CircularAcceleration<0, 1, 2>(2.0);

    RK4_stepper<2> stepper{dfunc};

    // Time loop
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 5'000'000; ++i)
    {
        stepper(state, 0.1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_span.count() << " seconds" << std::endl;
    std::cout << state[0] << ", " << state[1] << std::endl;
    std::cout << std::endl;

    return 0;
}
