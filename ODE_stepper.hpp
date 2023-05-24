
#ifndef ODE_STEPPER_HPP
#define ODE_STEPPER_HPP


#include <array>
#include <cassert>
#include <memory>
#include <iostream>
#include "stddef.h"
#include "function.hpp"

// ================================================================================================= //
// ===================================== ODE Stepper Interface ===================================== //
// ================================================================================================= //
template<size_t N>
class ODE_stepper
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    ODE_stepper() = default;
    ODE_stepper(const ODE_stepper& d) = default;
    ODE_stepper(ODE_stepper&& d) = default;
    ODE_stepper& operator=(const ODE_stepper& d) = default;
    ODE_stepper& operator=(ODE_stepper&& d) = default;
    virtual ~ODE_stepper() = default;

    void operator()(State<N>& state, const double& dt) const
    {
        step(state, dt);
    }
    virtual void step(State<N>& state, const double& dt) const = 0;
};


// ================================================================================================= //
// ========================================= Euler Method ========================================== //
// ================================================================================================= //
template<size_t N>
class Euler_stepper : public ODE_stepper<N>
{
private:
    //const std::function<State<N>(const State<N>&)> _dfunc;
    const std::shared_ptr<const Derivative<N>> _dfunc;
public:
    Euler_stepper(const std::shared_ptr<const Derivative<N>>& dfunc) : _dfunc{dfunc} {}

    template <typename... VS>
    requires std::is_constructible_v<CompositeDerivative<N>, VS...>
    Euler_stepper(VS&&... dfuncs) : _dfunc{std::make_shared<const CompositeDerivative<N>>(std::forward<VS>(dfuncs)...)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    Euler_stepper(const V& dfunc) : _dfunc{std::make_shared<const V>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    Euler_stepper(V&& dfunc) : _dfunc{std::make_shared<const V>(std::move(dfunc))} {}

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
class RK4_stepper : public ODE_stepper<N>
{
private:
    //const std::function<State<N>(const State<N>&)> _dfunc;
    const std::shared_ptr<const Derivative<N>> _dfunc;

public:
    RK4_stepper(const std::shared_ptr<const Derivative<N>>& dfunc) : _dfunc{dfunc} {}

    template <typename... VS>
    requires std::is_constructible_v<CompositeDerivative<N>, VS...>
    RK4_stepper(VS&&... dfuncs) : _dfunc{std::make_shared<const CompositeDerivative<N>>(std::forward<VS>(dfuncs)...)} {}

    template <typename V>
    requires std::is_constructible_v<WrappedDerivative<N>, V> && (!std::is_base_of_v<Derivative<N>, V>)
    RK4_stepper(const V& dfunc) : _dfunc{std::make_shared<const WrappedDerivative<N>>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    RK4_stepper(const V& dfunc) : _dfunc{std::make_shared<const V>(dfunc)} {}

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    RK4_stepper(V&& dfunc) : _dfunc{std::make_shared<const V>(std::move(dfunc))} {}

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

template <typename T>
requires std::is_base_of_v<Derivative<T::_N>, T>
class RK4_template_stepper : public ODE_stepper<T::_N>
{
private:
    const T _dfunc;

public:
    RK4_template_stepper(const T& dfunc) : _dfunc{dfunc} {}

    void inline step(State<T::_N>& state, const double& dt) const override
    {
        State<T::_N> k1 = _dfunc(state);
        State<T::_N> k2 = _dfunc(state + k1*dt/2.0);
        State<T::_N> k3 = _dfunc(state + k2*dt/2.0);
        State<T::_N> k4 = _dfunc(state + k3*dt);

        assert(state.size() == k1.size());
        assert(state.size() == k2.size());
        assert(state.size() == k3.size());
        assert(state.size() == k4.size());

        state += (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
    }
};

#endif // ODE_STEPPER_HPP