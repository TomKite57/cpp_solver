
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

template <size_t N>
class DynamicSolver : public Solver<N>
{
private:
    std::shared_ptr<ODE_stepper<N>> _ODE_stepper{nullptr};
    std::shared_ptr<Algebraic<N>> _prestep{nullptr};
    std::shared_ptr<Algebraic<N>> _poststep{nullptr};


public:
    DynamicSolver() = default;
    DynamicSolver(const std::shared_ptr<ODE_stepper<N>>& stepper) : _ODE_stepper{stepper} {}

    template <typename V>
    requires std::is_base_of_v<ODE_stepper<N>, V>
    DynamicSolver(const V& stepper): _ODE_stepper{std::make_shared<V>(stepper)} {}

    template <typename V>
    requires std::is_base_of_v<ODE_stepper<N>, V>
    DynamicSolver(V&& stepper): _ODE_stepper{std::make_shared<V>(std::move(stepper))} {}

    DynamicSolver(const DynamicSolver& d) = default;
    DynamicSolver(DynamicSolver&& d) = default;
    DynamicSolver& operator=(const DynamicSolver& d) = default;
    DynamicSolver& operator=(DynamicSolver&& d) = default;
    virtual ~DynamicSolver() = default;

    // Setters for the steppers
    void set_ODE_step(const std::shared_ptr<ODE_stepper<N>>& function) { _ODE_stepper = function; }
    void set_prestep(const std::shared_ptr<Algebraic<N>>& function) { _prestep = function; }
    void set_poststep(const std::shared_ptr<Algebraic<N>>& function) { _poststep = function; }

    template <typename V>
    requires std::is_base_of_v<ODE_stepper<N>, V>
    void set_ODE_step(const V& function) { _ODE_stepper = std::make_shared<V>(function); }

    template <typename V>
    requires std::is_base_of_v<ODE_stepper<N>, V>
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
// ========================================= Const solver ========================================== //
// ================================================================================================= //

template <size_t N>
class ConstSolver : public Solver<N>
{
private:
    const std::shared_ptr<const ODE_stepper<N>> _ODE_stepper;
    const std::shared_ptr<const Algebraic<N>> _prestep;
    const std::shared_ptr<const Algebraic<N>> _poststep;

public:
    template <typename ODE, typename PRE, typename POST>
    requires std::is_base_of_v<ODE_stepper<N>, ODE> && std::is_base_of_v<Algebraic<N>, PRE> && std::is_base_of_v<Algebraic<N>, POST>
    ConstSolver(const ODE& ode, const PRE& pre, const POST& post):
    _ODE_stepper{std::make_shared<const ODE>(ode)}, _prestep{std::make_shared<const PRE>(pre)}, _poststep{std::make_shared<const POST>(post)}
    {}


    ConstSolver(const ConstSolver& d) = default;
    ConstSolver(ConstSolver&& d) = default;
    ConstSolver& operator=(const ConstSolver& d) = default;
    ConstSolver& operator=(ConstSolver&& d) = default;
    virtual ~ConstSolver() = default;

    virtual void step(State<N>& state, const double& dt) const override
    {
        (*_prestep)(state, dt);
        (*_ODE_stepper)(state, dt);
        (*_poststep)(state, dt);
    }
};

// ================================================================================================= //
// ======================================== Template Solver ======================================== //
// ================================================================================================= //

template <size_t N, typename ODE, typename PRE, typename POST>
requires std::is_base_of_v<ODE_stepper<N>, ODE> && std::is_base_of_v<Algebraic<N>, PRE> && std::is_base_of_v<Algebraic<N>, POST> \
         && (N == ODE::_N) && (N == PRE::_N) && (N == POST::_N)
class TemplateSolver : public Solver<N>
{
private:
    const ODE _ODE_stepper;
    const PRE _prestep;
    const POST _poststep;

public:
    TemplateSolver(const ODE& ode, const PRE& pre, const POST& post):
    _ODE_stepper{ode}, _prestep{pre}, _poststep{post}
    {}

    TemplateSolver(const TemplateSolver& d) = default;
    TemplateSolver(TemplateSolver&& d) = default;
    TemplateSolver& operator=(const TemplateSolver& d) = default;
    TemplateSolver& operator=(TemplateSolver&& d) = default;
    virtual ~TemplateSolver() = default;

    virtual void step(State<N>& state, const double& dt) const override
    {
        _prestep(state, dt);
        _ODE_stepper(state, dt);
        _poststep(state, dt);
    }
};

template <typename ODE, typename PRE, typename POST>
requires std::is_base_of_v<ODE_stepper<ODE::_N>, ODE> && std::is_base_of_v<Algebraic<ODE::_N>, PRE> && std::is_base_of_v<Algebraic<ODE::_N>, POST>
auto MakeTemplateSolver(const ODE& ode, const PRE& pre, const POST& post) -> TemplateSolver<ODE::_N, std::decay_t<ODE>, std::decay_t<PRE>, std::decay_t<POST>>
{
    return TemplateSolver<ODE::_N, std::decay_t<ODE>, std::decay_t<PRE>, std::decay_t<POST>>(ode, pre, post);
}

#endif // SOLVER_HPP