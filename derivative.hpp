
#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP


#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <iostream>
#include "stddef.h"

#include "state.hpp"

// ================================================================================================= //
// ===================================== Derivative functions ====================================== //
// ================================================================================================= //
template<size_t N>
class Derivative
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    Derivative() = default;
    Derivative(const Derivative& d) = default;
    Derivative(Derivative&& d) = default;
    Derivative& operator=(const Derivative& d) = default;
    Derivative& operator=(Derivative&& d) = default;
    virtual ~Derivative() = default;

    State<N> operator()(const State<N>& state) const { return call(state); }

    virtual State<N> call(const State<N>& state) const = 0;
};


// ====================================== Dynamic derivatives ====================================== //
// Wrapper
template<size_t N>
class WrappedDerivative : public Derivative<N>
{
private:
    std::function<State<N>(const State<N>&)> _dfunc;

public:
    WrappedDerivative() = default;
    WrappedDerivative(const WrappedDerivative& d) = default;
    WrappedDerivative(WrappedDerivative&& d) = default;
    WrappedDerivative& operator=(const WrappedDerivative& d) = default;
    WrappedDerivative& operator=(WrappedDerivative&& d) = default;
    virtual ~WrappedDerivative() = default;

    WrappedDerivative(const std::function<State<N>(const State<N>&)>& dfunc) : _dfunc{dfunc} {}

    void set_dfunc(const std::function<State<N>(const State<N>&)>& dfunc) { _dfunc = dfunc; }

    virtual State<N> call(const State<N>& state) const override { return _dfunc(state); }
};

// Composite
template<size_t N>
class CompositeDerivative : public Derivative<N>
{
private:
    std::vector<std::shared_ptr<Derivative<N>>> _functions;

    // add_function_helpers
    template <typename V>
    requires std::constructible_from<WrappedDerivative<N>, V> && (!std::is_base_of_v<Derivative<N>, V>)
    void add_function_helper(const V& function)
    {
        _functions.push_back(std::make_shared<WrappedDerivative<N>>(function));
    }

    template <typename V>
    requires std::is_base_of_v<Derivative<N>, V>
    void add_function_helper(const V& function)
    {
        _functions.push_back(std::make_shared<V>(function));
    }

    void add_function_helper(const std::shared_ptr<Derivative<N>>& function)
    {
        _functions.push_back(function);
    }

    void add_function_helper() { return; }

public:
    CompositeDerivative() = default;
    CompositeDerivative(const CompositeDerivative& d) = default;
    CompositeDerivative(CompositeDerivative&& d) = default;
    CompositeDerivative& operator=(const CompositeDerivative& d) = default;
    CompositeDerivative& operator=(CompositeDerivative&& d) = default;
    virtual ~CompositeDerivative() = default;

    // Variadic constructor
    template <typename V, typename... VS>
    CompositeDerivative(const V& dfunc, const VS&... dfuncs): _functions{}
    {
        add_function(dfunc, dfuncs...);
    }

    // Variadic add_function
    template <typename V, typename... Args>
    void add_function(const V& dfunc, Args&&... args)
    {
        add_function_helper(dfunc);
        add_function(std::forward<Args>(args)...);
    }

    // add_function base case
    void add_function() { return; }

    virtual State<N> call(const State<N>& state) const override
    {
        State<N> result = State<N>(0.0);

        for (const auto& f : _functions)
            result += (*f)(state);

        return result;
    }
};

// ====================================== Template derivatives ===================================== //
// Wrapper
template<size_t N, typename V>
requires std::invocable<const V&, const State<N>&> && std::same_as<std::invoke_result_t<const V&, const State<N>&>, State<N>>
class TWrappedDerivative : public Derivative<N>
{
private:
    const V _dfunc;

public:
    TWrappedDerivative() = delete;
    TWrappedDerivative(const TWrappedDerivative& d) = default;
    TWrappedDerivative(TWrappedDerivative&& d) = default;
    TWrappedDerivative& operator=(const TWrappedDerivative& d) = default;
    TWrappedDerivative& operator=(TWrappedDerivative&& d) = default;
    virtual ~TWrappedDerivative() = default;

    TWrappedDerivative(const V& dfunc) : _dfunc{dfunc} {}

    virtual State<N> call(const State<N>& state) const override { return _dfunc(state); }
};

template<size_t N, typename V>
requires std::invocable<const V&, const State<N>&> && std::same_as<std::invoke_result_t<const V&, const State<N>&>, State<N>>
auto MakeWrappedDerivative(const V& dfunc) -> TWrappedDerivative<N, V>
{
    return TWrappedDerivative<N, V>(dfunc);
}

// Composite
template<size_t N, typename... Args>
class TCompositeDerivative : public Derivative<N>
{
private:
    const std::tuple<const Args...> _functions;

public:
    TCompositeDerivative() = delete;
    TCompositeDerivative(const TCompositeDerivative& d) = default;
    TCompositeDerivative(TCompositeDerivative&& d) = default;
    TCompositeDerivative& operator=(const TCompositeDerivative& d) = default;
    TCompositeDerivative& operator=(TCompositeDerivative&& d) = default;
    virtual ~TCompositeDerivative() = default;

    // Variadic constructor
    template <typename... VS>
    TCompositeDerivative(const VS&... dfuncs): _functions{dfuncs...} {}

    virtual State<N> call(const State<N>& state) const override
    {
        State<N> result = State<N>(0.0);

        std::apply([&](auto&&... args) { ((result += args(state)), ...); }, _functions);

        return result;
    }
};

template<size_t N, typename... Args>
auto MakeCompositeDerivative(const Args&... args) -> TCompositeDerivative<N, Args...>
{
    return TCompositeDerivative<N, Args...>(args...);
}

// ====================================== Concrete derivatives ===================================== //
// SHM
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
// ====================================== Do Nothing Function ====================================== //
// ================================================================================================= //
template<size_t N>
class DoNothingDerivative : public Derivative<N>
{
public:
    DoNothingDerivative() = default;
    DoNothingDerivative(const DoNothingDerivative& d) = default;
    DoNothingDerivative(DoNothingDerivative&& d) = default;
    DoNothingDerivative& operator=(const DoNothingDerivative& d) = default;
    DoNothingDerivative& operator=(DoNothingDerivative&& d) = default;
    virtual ~DoNothingDerivative() = default;

    State<N> call(const State<N>&) const override { return State<N>(0.0); }
};



#endif // FUNCTION_HPP