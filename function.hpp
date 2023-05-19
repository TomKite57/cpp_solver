
#ifndef FUNCTION_HPP
#define FUNCTION_HPP


#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <iostream>
#include "stddef.h"

template<size_t N>
class Function
{
public:
    Function() = default;
    Function(const Function& d) = default;
    Function(Function&& d) = default;
    Function& operator=(const Function& d) = default;
    Function& operator=(Function&& d) = default;
    virtual ~Function() = default;

    State<N> operator()(const State<N>& state) const
    {
        return call(state);
    }

    virtual State<N> call(const State<N>& state) const = 0;
};


template<size_t N>
class CompositeFunction : public Function<N>
{
private:
    std::vector<std::shared_ptr<Function<N>>> _functions;

public:
    CompositeFunction() = default;
    CompositeFunction(const CompositeFunction& d): _functions(d._functions) {};
    CompositeFunction(CompositeFunction&& d): _functions(std::move(d._functions)) {};
    CompositeFunction& operator=(const CompositeFunction& d) { _functions = d._functions; return *this;}
    CompositeFunction& operator=(CompositeFunction&& d) { _functions = std::move(d._functions); return *this; };
    virtual ~CompositeFunction() = default;

    //CompositeFunction(std::initializer_list<std::unique_ptr<Function<N>>> functions): _functions{functions} {}

    void add_function(std::shared_ptr<Function<N>> function)
    {
        _functions.push_back(function);
    }

    virtual State<N> call(const State<N>& state) const override
    {
        State<N> result = State<N>(0.0, N);

        for (size_t i = 0; i < _functions.size(); ++i)
        {
            result += (*_functions[i])(state);
        }

        return result;
    }
};


template<size_t x_ind, size_t dx_ind, size_t N>
class CircularAcceleration : public Function<N>
{
public:
    double _omega2;
    static constexpr size_t _x_ind = x_ind;
    static constexpr size_t _dx_ind = dx_ind;

public:
    CircularAcceleration(double omega2) : _omega2(omega2) {}

    CircularAcceleration(const CircularAcceleration& d): _omega2(d._omega2) {}
    CircularAcceleration(CircularAcceleration&& d): _omega2(d._omega2) {}
    CircularAcceleration& operator=(const CircularAcceleration& d){ _omega2 = d._omega2; return *this; }
    CircularAcceleration& operator=(CircularAcceleration&& d){ _omega2 = d._omega2; return *this; }
    virtual ~CircularAcceleration() {}

    State<N> call(const State<N>& state) const override
    {
        assert(state.size() > x_ind);
        assert(state.size() > dx_ind);

        State<N> dstate(0.0, N);
        dstate[x_ind] = state[dx_ind];
        dstate[dx_ind] = -_omega2*state[x_ind];
        return dstate;
    }
};

#endif // FUNCTION_HPP