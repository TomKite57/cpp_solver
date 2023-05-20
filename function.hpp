
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

    State<N> operator()(const State<N>& state) const { return call(state); }

    virtual State<N> call(const State<N>& state) const = 0;
};

template<size_t N>
class CompositeDerivative : public Derivative<N>
{
private:
    std::vector<std::shared_ptr<Derivative<N>>> _functions;

public:
    CompositeDerivative() = default;
    CompositeDerivative(const CompositeDerivative& d) = default;
    CompositeDerivative(CompositeDerivative&& d) = default;
    CompositeDerivative& operator=(const CompositeDerivative& d) = default;
    CompositeDerivative& operator=(CompositeDerivative&& d) = default;
    virtual ~CompositeDerivative() = default;

    // Variadic constructor
    template <typename V, typename... VS>
    requires std::is_base_of_v<Derivative<N>, V>
    CompositeDerivative(const V& dfunc, const VS&... dfuncs)
    {
        add_function(dfunc, dfuncs...);
    }

    // Base cases
    void add_function(const std::shared_ptr<Derivative<N>>& function)
    {
        _functions.push_back(function);
    }

    void add_function() { return; }

    // Recursive case
    template <typename V, typename... Args>
    requires std::is_base_of_v<Derivative<N>, V>
    void add_function(const V& dfunc, Args&&... args)
    {
        _functions.push_back(std::make_shared<V>(dfunc));
        add_function(std::forward<Args>(args)...);
    }

    template <typename V, typename... Args>
    requires std::is_base_of_v<Derivative<N>, V>
    void add_function(V&& dfunc, Args&&... args)
    {
        _functions.push_back(std::make_shared<V>(std::move(dfunc)));
        add_function(std::forward<Args>(args)...);
    }

    virtual State<N> call(const State<N>& state) const override
    {
        State<N> result = State<N>(0.0);

        for (const auto& f : _functions)
            result += (*f)(state);

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
        static_assert(state.size() > _x_ind);
        static_assert(state.size() > _dx_ind);

        State<N> dstate(0.0);
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

template<size_t ind, size_t N>
class Decrementor : public Algebraic<N>
{
public:
    Decrementor() = default;
    Decrementor(const Decrementor& d) = default;
    Decrementor(Decrementor&& d) = default;
    Decrementor& operator=(const Decrementor& d) = default;
    Decrementor& operator=(Decrementor&& d) = default;
    virtual ~Decrementor() = default;

    void call(State<N>& state, const double& dt) const override
    {
        state[ind] -= dt;
    }
};

template<size_t N>
class CompositeAlgebraic : public Algebraic<N>
{
private:
    std::vector<std::shared_ptr<Algebraic<N>>> _functions;

public:
    CompositeAlgebraic() = default;
    CompositeAlgebraic(const CompositeAlgebraic& d) = default;
    CompositeAlgebraic(CompositeAlgebraic&& d) = default;
    CompositeAlgebraic& operator=(const CompositeAlgebraic& d) = default;
    CompositeAlgebraic& operator=(CompositeAlgebraic&& d) = default;
    virtual ~CompositeAlgebraic() = default;

    // Variadic constructor
    template <typename V, typename... VS>
    requires std::is_base_of_v<Algebraic<N>, V>
    CompositeAlgebraic(const V& dfunc, const VS&... dfuncs)
    {
        add_function(dfunc, dfuncs...);
    }

    // Base case for adding functions
    void add_function(const std::shared_ptr<Algebraic<N>>& function)
    {
        _functions.push_back(function);
    }

    void add_function() { return; }

    // Recursive case for adding functions
    template <typename V, typename... VS>
    requires std::is_base_of_v<Algebraic<N>, V>
    void add_function(const V& dfunc, VS&&... dfuncs)
    {
        _functions.push_back(std::make_shared<V>(dfunc));
        add_function(std::forward<VS>(dfuncs)...);
    }

    template <typename V, typename... VS>
    requires std::is_base_of_v<Algebraic<N>, V>
    void add_function(V&& dfunc, VS&&... dfuncs)
    {
        _functions.push_back(std::make_shared<V>(std::move(dfunc)));
        add_function(std::forward<VS>(dfuncs)...);
    }

    virtual void call(State<N>& state, const double& dt) const override
    {
        for (const auto& f : _functions)
        {
            (*f)(state, dt);
        }
    }
};

template<size_t ind, size_t N>
class CircularPosition : public Algebraic<N>
{
public:
    const double _omega;
    static constexpr size_t _ind = ind;

public:
    CircularPosition(double omega) : _omega{omega} {}

    CircularPosition(const CircularPosition& d): _omega(d._omega) {}
    CircularPosition(CircularPosition&& d): _omega(d._omega) {}
    CircularPosition& operator=(const CircularPosition& d){ _omega = d._omega; return *this; }
    CircularPosition& operator=(CircularPosition&& d){ _omega = d._omega; return *this; }
    virtual ~CircularPosition() {}

    void call(State<N>& state, const double& dt) const override
    {
        static_assert(state.size() > _ind);

        state[_ind] = std::fmod(state[_ind] + _omega*dt, 2.0*M_PI);
    }
};


template<size_t N>
class DoNothing : public Algebraic<N>, public Derivative<N>
{
public:
    DoNothing() = default;
    DoNothing(const DoNothing& d) = default;
    DoNothing(DoNothing&& d) = default;
    DoNothing& operator=(const DoNothing& d) = default;
    DoNothing& operator=(DoNothing&& d) = default;
    virtual ~DoNothing() = default;

    void call(State<N>& state, const double& dt) const override { return; }

    State<N> call(const State<N>& state) const override { return State<N>(0.0); }
};

#endif // FUNCTION_HPP