
#ifndef FUNCTION_HPP
#define FUNCTION_HPP


#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <iostream>
#include "stddef.h"

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
    std::vector<std::shared_ptr<Derivative<N>>> _functions;

public:
    CompositeDerivative() = default;
    CompositeDerivative(const CompositeDerivative& d): _functions(d._functions) {};
    CompositeDerivative(CompositeDerivative&& d): _functions(std::move(d._functions)) {};
    CompositeDerivative& operator=(const CompositeDerivative& d) { _functions = d._functions; return *this;}
    CompositeDerivative& operator=(CompositeDerivative&& d) { _functions = std::move(d._functions); return *this; };
    virtual ~CompositeDerivative() = default;

    void add_function(std::shared_ptr<Derivative<N>> function)
    {
        _functions.push_back(function);
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

#endif // FUNCTION_HPP