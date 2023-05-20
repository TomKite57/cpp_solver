
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <memory>

#include "state.hpp"
#include "function.hpp"

// ================================================================================================= //
// =================================== General Solver Interface ==================================== //
// ================================================================================================= //
template <size_t N>
class Solver
{
public:
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
    std::unique_ptr<ODE_stepper<N>> _stepper;
    std::unique_ptr<Algebraic<N>> _prestep{nullptr};
    std::unique_ptr<Algebraic<N>> _poststep{nullptr};


public:
    DynamicSolver() = delete;
    DynamicSolver(std::unique_ptr<ODE_stepper<N>>&& stepper) : _stepper{std::move(stepper)} {}
    DynamicSolver(const DynamicSolver& d) = default;
    DynamicSolver(DynamicSolver&& d) = default;
    DynamicSolver& operator=(const DynamicSolver& d) = default;
    DynamicSolver& operator=(DynamicSolver&& d) = default;
    virtual ~DynamicSolver() = default;

    //void set_derivative(std::unique_ptr<Derivative<N>>&& function) { _derivative = function; }
    void set_prestep(std::unique_ptr<Algebraic<N>>&& function) { _prestep = std::move(function); }
    void set_poststep(std::unique_ptr<Algebraic<N>>&& function) { _poststep = std::move(function); }

    virtual void step(State<N>& state, const double& dt) const override
    {
        if (_prestep) (*_prestep)(state, dt);
        if (_stepper) (*_stepper)(state, dt);
        if (_poststep) (*_poststep)(state, dt);
    }

};

#endif // SOLVER_HPP