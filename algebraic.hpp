
#ifndef ALGEBRAIC_HPP
#define ALGEBRAIC_HPP


#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <iostream>
#include "stddef.h"


// ================================================================================================= //
// ====================================== Algebraic Functions ====================================== //
// ================================================================================================= //
template<size_t N>
class Algebraic
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    Algebraic() = default;
    Algebraic(const Algebraic& d) = default;
    Algebraic(Algebraic&& d) = default;
    Algebraic& operator=(const Algebraic& d) = default;
    Algebraic& operator=(Algebraic&& d) = default;
    virtual ~Algebraic() = default;

    void operator()(State<N>& state, const double& dt) const { return call(state, dt); }

    virtual void call(State<N>& state, const double& dt) const = 0;
};

// ======================================= Dynamic algebraic ======================================= //
// Wrapper
template<size_t N>
class WrappedAlgebraic : public Algebraic<N>
{
private:
    const std::function<void(State<N>&, const double&)> _func;

public:
    WrappedAlgebraic() = default;
    WrappedAlgebraic(const WrappedAlgebraic& d) = default;
    WrappedAlgebraic(WrappedAlgebraic&& d) = default;
    WrappedAlgebraic& operator=(const WrappedAlgebraic& d) = default;
    WrappedAlgebraic& operator=(WrappedAlgebraic&& d) = default;
    virtual ~WrappedAlgebraic() = default;

    WrappedAlgebraic(const std::function<void(State<N>&, const double&)>& func) : _func{func} {}

    void set_func(const std::function<void(State<N>&, const double&)>& func) { _func = func; }

    virtual void call(State<N>& state, const double& dt) const override { _func(state, dt); }
};

// Composite
template<size_t N>
class CompositeAlgebraic : public Algebraic<N>
{
private:
    std::vector<std::shared_ptr<const Algebraic<N>>> _functions;

    // add_function_helpers
    template <typename V>
    requires std::constructible_from<WrappedAlgebraic<N>, V> && (!std::is_base_of_v<Algebraic<N>, V>)
    void add_function_helper(const V& function)
    {
        _functions.push_back(std::make_shared<const WrappedAlgebraic<N>>(function));
    }

    template <typename V>
    requires std::is_base_of_v<Algebraic<N>, V>
    void add_function_helper(const V& function)
    {
        _functions.push_back(std::make_shared<const V>(function));
    }

    void add_function_helper(const std::shared_ptr<Algebraic<N>>& function)
    {
        _functions.push_back(function);
    }

    void add_function_helper() { return; }

public:
    CompositeAlgebraic() = default;
    CompositeAlgebraic(const CompositeAlgebraic& d) = default;
    CompositeAlgebraic(CompositeAlgebraic&& d) = default;
    CompositeAlgebraic& operator=(const CompositeAlgebraic& d) = default;
    CompositeAlgebraic& operator=(CompositeAlgebraic&& d) = default;
    virtual ~CompositeAlgebraic() = default;

    // Variadic constructor
    template <typename V, typename... VS>
    CompositeAlgebraic(const V& dfunc, const VS&... dfuncs): _functions{}
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

    virtual void call(State<N>& state, const double& dt) const override
    {
        for (const auto& f : _functions)
        {
            (*f)(state, dt);
        }
    }
};

// ======================================= Template algebraic ====================================== //
// Wrapper
template<size_t N, typename V>
class TWrappedAlgebraic : public Algebraic<N>
{
private:
    const V _func;

public:
    TWrappedAlgebraic() = delete;
    TWrappedAlgebraic(const TWrappedAlgebraic& d) = default;
    TWrappedAlgebraic(TWrappedAlgebraic&& d) = default;
    TWrappedAlgebraic& operator=(const TWrappedAlgebraic& d) = default;
    TWrappedAlgebraic& operator=(TWrappedAlgebraic&& d) = default;
    ~TWrappedAlgebraic() = default;

    TWrappedAlgebraic(const V& func) : _func{func} {}

    void call(State<N>& state, const double& dt) const override { _func(state, dt); }
};

template<size_t N, typename V>
requires std::invocable<V, State<N>&, const double&>
auto MakeWrappedAlgebraic(const V& func) -> TWrappedAlgebraic<N, V>
{
    return TWrappedAlgebraic<N, V>(func);
}

// Composite
template<size_t N, typename... Args>
class TCompositeAlgebraic : public Algebraic<N>
{
private:
    const std::tuple<const Args...> _functions;

public:
    TCompositeAlgebraic() = delete;
    TCompositeAlgebraic(const TCompositeAlgebraic& d) = default;
    TCompositeAlgebraic(TCompositeAlgebraic&& d) = default;
    TCompositeAlgebraic& operator=(const TCompositeAlgebraic& d) = default;
    TCompositeAlgebraic& operator=(TCompositeAlgebraic&& d) = default;
    ~TCompositeAlgebraic() = default;

    TCompositeAlgebraic(const Args&... dfuncs) : _functions{dfuncs...} {}

    void call(State<N>& state, const double& dt) const override
    {
        std::apply([&](const auto&... f) { (f(state, dt), ...); }, _functions);
    }
};

template<size_t N, typename... Args>
auto MakeCompositeAlgebraic(const Args&... dfuncs) -> TCompositeAlgebraic<N, Args...>
{
    return TCompositeAlgebraic<N, Args...>(dfuncs...);
}

// ======================================= Concrete algebraic ====================================== //
// Increment time
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

// Decrement time
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


// ================================================================================================= //
// ====================================== Do Nothing Function ====================================== //
// ================================================================================================= //
template<size_t N>
class DoNothingAlgebraic : public Algebraic<N>
{
public:
    DoNothingAlgebraic() = default;
    DoNothingAlgebraic(const DoNothingAlgebraic& d) = default;
    DoNothingAlgebraic(DoNothingAlgebraic&& d) = default;
    DoNothingAlgebraic& operator=(const DoNothingAlgebraic& d) = default;
    DoNothingAlgebraic& operator=(DoNothingAlgebraic&& d) = default;
    virtual ~DoNothingAlgebraic() = default;

    void call(State<N>&, const double&) const override { return; }
};


#endif // FUNCTION_HPP