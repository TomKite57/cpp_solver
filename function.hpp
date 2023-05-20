
#ifndef FUNCTION_HPP
#define FUNCTION_HPP


#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <iostream>
#include "stddef.h"

// ================================================================================================= //
// ===================================== Derivative functions ====================================== //
// ================================================================================================= //
template<size_t N>
class Derivative
{
public:
    Derivative() = default;
    Derivative(const Derivative& d) = default;
    Derivative(Derivative&& d) = default;
    Derivative& operator=(const Derivative& d) = default;
    Derivative& operator=(Derivative&& d) = default;
    virtual ~Derivative() = default;

    virtual void add_function(std::unique_ptr<Derivative<N>>&& function) { throw std::runtime_error("Class Not Composite"); }

    State<N> operator()(const State<N>& state) const
    {
        return call(state);
    }

    virtual State<N> call(const State<N>& state) const = 0;
};

template<size_t N>
class CompositeDerivative : public Derivative<N>
{
private:
    std::vector<std::unique_ptr<Derivative<N>>> _functions;

public:
    CompositeDerivative() = default;
    CompositeDerivative(const CompositeDerivative& d) = default;
    CompositeDerivative(CompositeDerivative&& d) = default;
    CompositeDerivative& operator=(const CompositeDerivative& d) = default;
    CompositeDerivative& operator=(CompositeDerivative&& d) = default;
    virtual ~CompositeDerivative() = default;

    void add_function(std::unique_ptr<Derivative<N>>&& function) override
    {
        _functions.push_back(std::move(function));
    }

    virtual State<N> call(const State<N>& state) const override
    {
        State<N> result = State<N>(0.0, N);

        for (const auto& f : _functions)
        {
            result += (*f)(state);
        }

        return result;
    }
};


template<size_t x_ind, size_t dx_ind, size_t N>
class CircularAcceleration : public Derivative<N>
{
public:
    const double _omega2;
    static constexpr size_t _x_ind = x_ind;
    static constexpr size_t _dx_ind = dx_ind;

public:
    CircularAcceleration(double omega2) : _omega2{omega2} {}

    CircularAcceleration(const CircularAcceleration& d): _omega2(d._omega2) {}
    CircularAcceleration(CircularAcceleration&& d): _omega2(d._omega2) {}
    CircularAcceleration& operator=(const CircularAcceleration& d){ _omega2 = d._omega2; return *this; }
    CircularAcceleration& operator=(CircularAcceleration&& d){ _omega2 = d._omega2; return *this; }
    virtual ~CircularAcceleration() {}

    State<N> call(const State<N>& state) const override
    {
        assert(state.size() > _x_ind);
        assert(state.size() > _dx_ind);

        State<N> dstate(0.0, N);
        dstate[_x_ind] = state[_dx_ind];
        dstate[_dx_ind] = -_omega2*state[_x_ind];
        return dstate;
    }
};


// ================================================================================================= //
// ====================================== Algebraic Functions ====================================== //
// ================================================================================================= //
template<size_t N>
class Algebraic
{
public:
    Algebraic() = default;
    Algebraic(const Algebraic& d) = default;
    Algebraic(Algebraic&& d) = default;
    Algebraic& operator=(const Algebraic& d) = default;
    Algebraic& operator=(Algebraic&& d) = default;
    virtual ~Algebraic() = default;

    void operator()(State<N>& state, const double& dt) const
    {
        return call(state, dt);
    }

    virtual void call(State<N>& state, const double& dt) const = 0;
};

template<size_t ind, size_t N>
class Incrementor : public Algebraic<N>
{
public:
    Incrementor() = default;
    Incrementor(const Incrementor& d) = default;
    Incrementor(Incrementor&& d) = default;
    Incrementor& operator=(const Incrementor& d) = default;
    Incrementor& operator=(Incrementor&& d) = default;
    virtual ~Incrementor() = default;

    void call(State<N>& state, const double& dt) const override
    {
        state[ind] += dt;
    }
};



#endif // FUNCTION_HPP