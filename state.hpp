
#ifndef STATE_HPP
#define STATE_HPP

#include <array>
#include <cassert>
#include <iostream>
#include "stddef.h"

template<size_t NN>
class State
{
public:
    static constexpr size_t N = NN;

    friend State<N> operator*(double scalar, const State<N>& state)
    {
        return state*scalar;
    }

    friend std::ostream& operator<<(std::ostream& os, const State<N>& state)
    {
        os << "State(";
        for (size_t i = 0; i < N; ++i)
        {
            os << state[i];
            if (i < N-1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }

    // Constructor that initializes all elements of the array to the given value
    State(double initial_value) : _data{ }
    {
        _data.fill(initial_value);
    }

    State(double initial_value, size_t size) : _data{ }
    {
        assert(size == N);
        _data.fill(initial_value);
    }

    // Constructor that initializes the array with the given values
    State(std::initializer_list<double> values) : _data{ }
    {
        assert(values.size() == N);

        size_t i = 0;
        for (auto value : values)
        {
            _data[i] = value;
            ++i;
        }
    }

    const double& operator[](size_t index) const
    {
        assert(index < N);
        assert(index >= 0);
        return _data[index];
    }

    double& operator[](size_t index)
    {
        assert(index < N);
        assert(index >= 0);
        return _data[index];
    }

    // Piecewise addition operator
    State<N> operator+(const State<N>& other) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] + other._data[i];
        }
        return result;
    }

    // Piecewise subtraction operator
    State<N> operator-(const State<N>& other) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] - other._data[i];
        }
        return result;
    }

    // Piecewise multiplication operator
    State<N> operator*(const State<N>& other) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] * other._data[i];
        }
        return result;
    }

    // Scaling operator
    State<N> operator*(double scalar) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] * scalar;
        }
        return result;
    }

    State<N> operator/(const State<N>& other) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] / other._data[i];
        }
    return result;
    }

    // Scalar division operator
    State<N> operator/(double scalar) const
    {
        State<N> result(0.0);
        for (size_t i = 0; i < N; ++i)
        {
            result._data[i] = _data[i] / scalar;
        }
        return result;
    }

    void operator+=(const State<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] += other._data[i];
        }
    }

    // Accessor functions
    size_t size() const
    {
        return N;
    }

private:
    std::array<double, N> _data;
};

# endif // STATE_HPP