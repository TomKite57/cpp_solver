
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <memory>

#include "state.hpp"
#include "algebraic.hpp"
#include "derivative.hpp"

// ================================================================================================= //
// =================================== General Solver Interface ==================================== //
// ================================================================================================= //
template <size_t N>
class Solver
{
public:
    static constexpr size_t size() { return N; }
    static constexpr size_t _N = N;

    Solver() = default;
    Solver(const Solver& s) = default;
    Solver(Solver&& s) = default;
    Solver& operator=(const Solver& s) = default;
    Solver& operator=(Solver&& s) = default;
    virtual ~Solver() = default;

    void operator()(State<N>& state, const double& dt) const
    {
        step(state, dt);
    }
    virtual void step(State<N>& state, const double& dt) const = 0;
};

// PrePostSolver
template <size_t N>
class PrePostSolver : public Solver<N>
{
private:
    std::shared_ptr<ODEStepper<N>> _ODE_stepper{nullptr};
    std::shared_ptr<Algebraic<N>> _prestep{nullptr};
    std::shared_ptr<Algebraic<N>> _poststep{nullptr};


public:
    PrePostSolver() = default;
    PrePostSolver(const std::shared_ptr<ODEStepper<N>>& stepper) : _ODE_stepper{stepper} {}

    template <typename V>
    requires std::is_base_of_v<ODEStepper<N>, V>
    PrePostSolver(const V& stepper): _ODE_stepper{std::make_shared<V>(stepper)} {}

    template <typename V>
    requires std::is_base_of_v<ODEStepper<N>, V>
    PrePostSolver(V&& stepper): _ODE_stepper{std::make_shared<V>(std::move(stepper))} {}

    PrePostSolver(const PrePostSolver& d) = default;
    PrePostSolver(PrePostSolver&& d) = default;
    PrePostSolver& operator=(const PrePostSolver& d) = default;
    PrePostSolver& operator=(PrePostSolver&& d) = default;
    virtual ~PrePostSolver() = default;

    // Setters for the steppers
    void set_ODE_step(const std::shared_ptr<ODEStepper<N>>& function) { _ODE_stepper = function; }
    void set_prestep(const std::shared_ptr<Algebraic<N>>& function) { _prestep = function; }
    void set_poststep(const std::shared_ptr<Algebraic<N>>& function) { _poststep = function; }

    template <typename V>
    requires std::is_base_of_v<ODEStepper<N>, V>
    void set_ODE_step(const V& function) { _ODE_stepper = std::make_shared<V>(function); }

    template <typename V>
    requires std::is_base_of_v<ODEStepper<N>, V>
    void set_ODE_step(V&& function) { _ODE_stepper = std::make_shared<V>(std::move(function)); }

    // Template versions of the prestep
    template <typename V>
    requires std::is_base_of_v<Algebraic<N>, V>
    void set_prestep(const V& function) { _prestep = std::make_shared<V>(function); }

    template <typename V>
    requires std::is_base_of_v<Algebraic<N>, V>
    void set_prestep(V&& function) { _prestep = std::make_shared<V>(std::move(function)); }

    // Template versions of the poststep
    template <typename V>
    requires std::is_base_of_v<Algebraic<N>, V>
    void set_poststep(const V& function) { _poststep = std::make_shared<V>(function); }

    template <typename V>
    requires std::is_base_of_v<Algebraic<N>, V>
    void set_poststep(V&& function) { _poststep = std::make_shared<V>(std::move(function)); }

    virtual void step(State<N>& state, const double& dt) const override
    {
        if (_prestep) (*_prestep)(state, dt);
        if (_ODE_stepper) (*_ODE_stepper)(state, dt);
        if (_poststep) (*_poststep)(state, dt);
    }
};

// ================================================================================================= //
// ======================================== Template Solver ======================================== //
// ================================================================================================= //

template <size_t N, typename ODE, typename PRE, typename POST>
requires std::is_base_of_v<ODEStepper<N>, ODE> && std::is_base_of_v<Algebraic<N>, PRE> && std::is_base_of_v<Algebraic<N>, POST> \
         && (N == ODE::_N) && (N == PRE::_N) && (N == POST::_N)
class TPrePostSolver : public Solver<N>
{
private:
    const ODE _ODE_stepper;
    const PRE _prestep;
    const POST _poststep;

public:
    TPrePostSolver(const ODE& ode, const PRE& pre, const POST& post):
    _ODE_stepper{ode}, _prestep{pre}, _poststep{post}
    {}

    TPrePostSolver(const TPrePostSolver& d) = default;
    TPrePostSolver(TPrePostSolver&& d) = default;
    TPrePostSolver& operator=(const TPrePostSolver& d) = default;
    TPrePostSolver& operator=(TPrePostSolver&& d) = default;
    virtual ~TPrePostSolver() = default;

    virtual void step(State<N>& state, const double& dt) const override
    {
        _prestep(state, dt);
        _ODE_stepper(state, dt);
        _poststep(state, dt);
    }
};

template <typename ODE, typename PRE, typename POST>
requires std::is_base_of_v<ODEStepper<ODE::_N>, ODE> && std::is_base_of_v<Algebraic<ODE::_N>, PRE> && std::is_base_of_v<Algebraic<ODE::_N>, POST>
auto MakePrePostSolver(const ODE& ode, const PRE& pre, const POST& post) -> TPrePostSolver<ODE::_N, std::decay_t<ODE>, std::decay_t<PRE>, std::decay_t<POST>>
{
    return TPrePostSolver<ODE::_N, std::decay_t<ODE>, std::decay_t<PRE>, std::decay_t<POST>>(ode, pre, post);
}

#endif // SOLVER_HPP