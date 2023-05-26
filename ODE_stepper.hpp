
#ifndef ODE_STEPPER_HPP
#define ODE_STEPPER_HPP

#include <array>
#include <cassert>
#include <memory>
#include <iostream>
#include "stddef.h"

#include "state.hpp"
#include "algebraic.hpp"
#include "derivative.hpp"

// ================================================================================================= //
// ===================================== ODE Stepper Interface ===================================== //
// ================================================================================================= //
template<size_t N>
class ODEStepper: public Algebraic<N>
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    ODEStepper() = default;
    ODEStepper(const ODEStepper& d) = default;
    ODEStepper(ODEStepper&& d) = default;
    ODEStepper& operator=(const ODEStepper& d) = default;
    ODEStepper& operator=(ODEStepper&& d) = default;
    virtual ~ODEStepper() = default;

    void operator()(State<N>& state, const double& dt) const
    {
        step(state, dt);
    }

    void call(State<N>& state, const double& dt) const override
    {
        step(state, dt);
    }

    virtual void step(State<N>& state, const double& dt) const = 0;
};


// ================================================================================================= //
// ========================================= Euler Method ========================================== //
// ================================================================================================= //
template<size_t N>
class EulerStepper : public ODEStepper<N>
{
private:
    //const std::function<State<N>(const State<N>&)> _dfunc;
    const std::shared_ptr<const Derivative<N>> _dfunc;
public:
    EulerStepper(const std::shared_ptr<const Derivative<N>>& dfunc) : _dfunc{dfunc} {}

    template <typename... VS>
    requires std::is_constructible_v<CompositeDerivative<N>, VS...>
    EulerStepper(VS&&... dfuncs) : _dfunc{std::make_shared<const CompositeDerivative<N>>(std::forward<VS>(dfuncs)...)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    EulerStepper(const V& dfunc) : _dfunc{std::make_shared<const V>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    EulerStepper(V&& dfunc) : _dfunc{std::make_shared<const V>(std::move(dfunc))} {}

    void inline step(State<N>& state, const double& dt) const override
    {
        State<N> dstate = _dfunc(state);

        assert(state.size() == dstate.size());

        state += dstate*dt;
    }
};

// ================================================================================================= //
// ========================================== RK4 Method =========================================== //
// ================================================================================================= //
template<size_t N>
class RK4Stepper : public ODEStepper<N>
{
private:
    //const std::function<State<N>(const State<N>&)> _dfunc;
    const std::shared_ptr<const Derivative<N>> _dfunc;

public:
    RK4Stepper(const std::shared_ptr<const Derivative<N>>& dfunc) : _dfunc{dfunc} {}

    template <typename... VS>
    requires std::is_constructible_v<CompositeDerivative<N>, VS...>
    RK4Stepper(VS&&... dfuncs) : _dfunc{std::make_shared<const CompositeDerivative<N>>(std::forward<VS>(dfuncs)...)} {}

    template <typename V>
    requires std::is_constructible_v<WrappedDerivative<N>, V> && (!std::is_base_of_v<Derivative<N>, V>)
    RK4Stepper(const V& dfunc) : _dfunc{std::make_shared<const WrappedDerivative<N>>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    RK4Stepper(const V& dfunc) : _dfunc{std::make_shared<const V>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    RK4Stepper(V&& dfunc) : _dfunc{std::make_shared<const V>(std::move(dfunc))} {}

    void inline step(State<N>& state, const double& dt) const override
    {
        State<N> k1 = (*_dfunc)(state);
        State<N> k2 = (*_dfunc)(state + k1*dt/2.0);
        State<N> k3 = (*_dfunc)(state + k2*dt/2.0);
        State<N> k4 = (*_dfunc)(state + k3*dt);

        assert(state.size() == k1.size());
        assert(state.size() == k2.size());
        assert(state.size() == k3.size());
        assert(state.size() == k4.size());

        state += (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
    }
};

// ================================================================================================= //
// ========================================= Template test ========================================= //
// ================================================================================================= //

template <size_t N, typename T>
requires std::invocable<T, State<N>>
class TRK4Stepper : public ODEStepper<N>
{
private:
    const T _dfunc;

public:
    TRK4Stepper(const T& dfunc) : _dfunc{dfunc} {}

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

template <size_t N, typename T>
auto MakeRK4Stepper(const T& dfunc) -> TRK4Stepper<N, T>
{
    return TRK4Stepper<N, T>(dfunc);
}

#endif // ODE_STEPPER_HPP