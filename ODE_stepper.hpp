
#ifndef ODE_STEPPER_HPP
#define ODE_STEPPER_HPP


#include <array>
#include <cassert>
#include <iostream>
#include "stddef.h"

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

#endif // ODE_STEPPER_HPP