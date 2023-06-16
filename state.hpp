
#ifndef STATE_HPP
#define STATE_HPP

#include <array>
#include <cassert>
#include <iostream>
#include "stddef.h"


// ================================================================================================= //
// ========================================== State Class ========================================== //
// ================================================================================================= //
template<size_t N>
class State
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    bool has_inf() const
    {
        for (size_t i = 0; i < N; ++i)
            if (std::isinf(_data[i]) || std::isinf(-_data[i]))
                return true;
        return false;
    }

    bool has_nan() const
    {
        for (size_t i = 0; i < N; ++i)
            if (std::isnan(_data[i]) || std::isnan(-_data[i]))
                return true;
        return false;
    }

    friend State<N> operator*(double scalar, const State<N>& state)
    {
        return state*scalar;
    }

    friend std::ostream& operator<<(std::ostream& os, const State<N>& state)
    {
        for (size_t i = 0; i < N; ++i)
        {
            os << state[i];
            if (i < N-1)
            {
                os << ", ";
            }
        }
        return os;
    }

    friend void print(const State<N>& state)
    {
        std::cout << "State(" << state << ")\n";
    }

    // Constructor that initializes all elements of the array to the given value
    State(double initial_value) : _data{ }
    {
        _data.fill(initial_value);
    }

    // Constructor that initializes the array with the given values
    State(const std::initializer_list<double>& values) : _data{ }
    {
        assert(values.size() == N);

        size_t i = 0;
        for (const auto& value : values)
        {
            _data[i] = value;
            ++i;
        }
    }

    const double& operator[](size_t index) const
    {
        assert(index < N);
        return _data[index];
    }

    double& operator[](size_t index)
    {
        assert(index < N);
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

    void operator-=(const State<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] -= other._data[i];
        }
    }

    void operator*=(const State<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] *= other._data[i];
        }
    }

    void operator*=(double scalar)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] *= scalar;
        }
    }

    void operator/=(const State<N>& other)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] /= other._data[i];
        }
    }

    void operator/=(double scalar)
    {
        for (size_t i = 0; i < N; ++i)
        {
            _data[i] /= scalar;
        }
    }

private:
    std::array<double, N> _data;
};

# endif // STATE_HPP