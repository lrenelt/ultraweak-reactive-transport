// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_TEST_PROBLEMS_HH
#define DUNE_ULTRAWEAK_TEST_PROBLEMS_HH

#include <cmath>

#include "dune/pdelab.hh"

namespace BoundaryProfiles {
  using RF = double;

  RF C2bump(const RF& x) {
    return (x>=0.25)*(x<=0.75)*(256*pow(x,4) - 512*pow(x,3) + 352*pow(x,2) - 96*x + 9);
  }

  RF sineBump(const RF& x, const unsigned int n=3) {
    return pow(sin(n*M_PI*x), 2);
  }

  RF L2bump(const RF& x, const unsigned int n=3) {
    return pow(sin(n*M_PI*x),2) >= 0.35;
  }
}

// TODO: solve in here?
template<typename BaseProblem, typename DGFType>
class DarcyVelocityAdapter : public BaseProblem
{
public:
  using Traits = typename BaseProblem::Traits;

  DarcyVelocityAdapter(const BaseProblem& baseProblem, const DGFType& velocityDgf)
    : BaseProblem(baseProblem), baseProblem_(baseProblem), velocityDgf_(velocityDgf) {}

  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    const auto& insidePos = el.geometryInInside().global(x);

    const auto& vel = b(el.inside(), insidePos);
    const auto val = vel.dot(el.unitOuterNormal(x));
    const double tol = 1e-10;
    if (val > tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
    else if (val < -tol)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::None;
  }

  template<typename Element, typename X>
  auto b (const Element& el, const X& x) const
  {
    typename Traits::RangeType ret(0.0);
    velocityDgf_.evaluate(el, x, ret);
    return ret;
  }

private:
  const BaseProblem& baseProblem_;
  const DGFType& velocityDgf_;
};

template <typename GV, typename RF_>
class PureTransportProblem : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF_>
{
public:
  using RF = RF_;
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;
  using Traits = typename Base::Traits;

  PureTransportProblem(Dune::ParameterTree pTree) :
    Base(), pTree_(pTree),
    reaction_(pTree_.get<RF>("problem.reaction")) {}

  // no diffusion
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    return typename Traits::PermTensorType(0.0);
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    using std::min, std::max;
    if (std::min(global[0],global[1]) < 1e-5)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }

  // curved velocity field
  template<typename Element, typename X>
  auto b (const Element& el, const X& x) const
  {
    using std::sqrt;
    return typename Traits::RangeType({sqrt(0.5), sqrt(0.5)});
  }

  // reaction coefficient
  template<typename Element, typename X>
  RF c (const Element& el, const X& x) const
  {
    const double halfReactionBlockHeight = 0.1;
    const auto& global = el.geometry().global(x);
    using std::abs;
    return reaction_ * (abs(global[1]-0.5) < halfReactionBlockHeight);
  }

  // Dirichlet condition
  template<typename Element, typename X>
  RF g (const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    return BoundaryProfiles::sineBump(global[1], pTree_.get<int>("problem.nInflowBumps"));
  }

private:
  Dune::ParameterTree pTree_;
  RF reaction_;
};

#endif // DUNE_ULTRAWEAK_TEST_PROBLEMS_HH
