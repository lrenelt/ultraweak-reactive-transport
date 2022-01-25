// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_CATALYSATOR_PROBLEMS_HH
#define DUNE_ULTRAWEAK_CATALYSATOR_PROBLEMS_HH

#include <dune/pdelab.hh>

/**
   Pressure gradient from left to right with low permeability block at
   0.4 < x < 0.6
*/
template <typename GV, typename RF>
class CatalysatorProblem1 : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>
{
public:
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>;
  using Traits = typename Base::Traits;

private:
  const RF halfReactionBlockWidth_ = 0.1;

public:

  CatalysatorProblem1(Dune::ParameterTree& pTree) : Base(),
    pTree_(pTree), I(0.0), minPerm_(pTree_.get<double>("darcy.min_permeability")) {
    // precompute unity tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      I[i][i] = 1.0;
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    if (global[0] < 1e-10 or global[0] > 1-1e-10)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  // Permeability tensor (isotropic)
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    const auto& global = el.geometry().global(x);

    using std::abs;
    if (abs(global[0] - 0.5) < halfReactionBlockWidth_)
      return minPerm_ * I;
    else
      return I;
  }

  // Dirichlet condition
  template<typename Element, typename X>
  auto g (const Element& el, const X& x) const
  {
    const auto gval = 1.0 - el.geometry().global(x)[0];
    return -gval;
  }

  //! Neumann boundary condition
  template<typename Intersection, typename X>
  auto j (const Intersection& is, const X& x) const
  {
    return 0.0;
  }

private:
  Dune::ParameterTree& pTree_;
  typename Traits::PermTensorType I;
  const RF minPerm_;
};


/**
   Pressure gradient from left to right with low permeability block at
   the bottom
*/
template <typename GV, typename RF>
class CatalysatorProblem2 : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>
{
public:
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>;
  using Traits = typename Base::Traits;

private:
  const RF reactionBlockHeight_ = 0.25;

public:

  CatalysatorProblem2(Dune::ParameterTree& pTree) : Base(),
    pTree_(pTree), I(0.0), minPerm_(pTree_.get<RF>("darcy.min_permeability")) {
    // precompute unity tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      I[i][i] = 1.0;
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    if ((global[0] < 1e-10 or global[0] > 1-1e-10) and
        global[1] > reactionBlockHeight_)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  // Permeability tensor (isotropic)
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    const auto& global = el.geometry().global(x);

    using std::abs;
    if (global[1] < reactionBlockHeight_)
      return minPerm_ * I;
    else
      return I;
  }

  // Dirichlet condition
  template<typename Element, typename X>
  auto g (const Element& el, const X& x) const
  {
    const auto gval = 1.0 - el.geometry().global(x)[0];
    return -gval;
  }

  //! Neumann boundary condition
  template<typename Intersection, typename X>
  auto j (const Intersection& is, const X& x) const
  {
    return 0.0;
  }

private:
  Dune::ParameterTree& pTree_;
  typename Traits::PermTensorType I;
  const RF minPerm_;
};


/**
   Pressure gradient from top left to bottom right with a horizontal
   low permeability block in the middle
*/
template <typename GV, typename RF>
class CatalysatorProblem3 : public Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>
{
public:
  using Base = Dune::PDELab::ConvectionDiffusionModelProblem<GV, RF>;
  using Traits = typename Base::Traits;

private:
  const RF halfReactionBlockHeight_ = 0.1;
  const RF openingHeight_ = 1./3;

public:

  CatalysatorProblem3(Dune::ParameterTree& pTree) : Base(),
    pTree_(pTree), I(0.0), minPerm_(pTree_.get<RF>("darcy.min_permeability")) {
    // precompute unity tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      I[i][i] = 1.0;
  }

  // Boundary condition type
  template<typename Element, typename X>
  auto bctype(const Element& el, const X& x) const
  {
    auto global = el.geometry().global(x);
    if (((global[0] < 1e-10) and (global[1] > 1-openingHeight_) ) or
        ((global[0] > 1 - 1e-10) and (global[1] < openingHeight_)))
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  // Permeability tensor (isotropic)
  template<typename Element, typename X>
  auto A (const Element& el, const X& x) const
  {
    const auto& global = el.geometry().global(x);

    using std::abs;
    if (abs(global[1] - 0.5) < halfReactionBlockHeight_)
      return minPerm_ * I;
    else
      return I;
  }

  // Dirichlet condition
  template<typename Element, typename X>
  auto g (const Element& el, const X& x) const
  {
    const auto gval = 1.0 - el.geometry().global(x)[0];
    return -gval;
  }

  //! Neumann boundary condition
  template<typename Intersection, typename X>
  auto j (const Intersection& is, const X& x) const
  {
    return 0.0;
  }

private:
  Dune::ParameterTree& pTree_;
  typename Traits::PermTensorType I;
  const RF minPerm_;
};

#endif  // DUNE_ULTRAWEAK_CATALYSATOR_PROBLEMS_HH
