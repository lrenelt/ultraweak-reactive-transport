#ifndef DUNE_ULTRAWEAK_TRANSPORT_DISCRETEGRIDFUNCTIONTRANSPORTDIFFOP_HH
#define DUNE_ULTRAWEAK_TRANSPORT_DISCRETEGRIDFUNCTIONTRANSPORTDIFFOP_HH

#include <vector>

#include <dune/pdelab.hh>

// TODO: remove
using namespace Dune;
using namespace PDELab;

template<typename GFS, typename X, typename Problem>
class DiscreteGridFunctionTransportDiffOp : public GridFunctionBase<
  GridFunctionTraits<typename GFS::Traits::GridViewType,
                     typename GFS::Traits::FiniteElementType::Traits::LocalBasisType
                       ::Traits::RangeFieldType,
                     GFS::Traits::FiniteElementType::Traits::LocalBasisType
                       ::Traits::dimRange,
                     typename GFS::Traits::FiniteElementType::Traits::LocalBasisType
                       ::Traits::RangeType>,
  DiscreteGridFunctionTransportDiffOp<GFS,X, Problem>>
{
private:
  using GV = typename GFS::Traits::GridViewType;
  using LBTraits = typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits;
  using RF = typename LBTraits::RangeFieldType;
  using RangeType = typename LBTraits::RangeType;
;
public:
  using Traits = GridFunctionTraits<GV,RF,LBTraits::dimRange,RangeType>;

private:
  using BaseT = GridFunctionBase<Traits,
                                 DiscreteGridFunctionTransportDiffOp<GFS,X,Problem>>;

  using DGFScalar = DiscreteGridFunction<GFS,X>;
  using DGFGradient = DiscreteGridFunctionGradient<GFS,X>;

public:
  DiscreteGridFunctionTransportDiffOp(const GFS& gfs, const X& x_, const Problem& problem, const bool rescaling = false) :
    problem_(problem),
    rescaling_(rescaling),
    dgfScalar_(gfs, x_),
    dgfGradient_(gfs, x_) {}

  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& x,
                       typename Traits::RangeType& y) const
  {
    dgfScalar_.evaluate(e,x,phi);
    dgfGradient_.evaluate(e,x,gradphi);

    //evaluate data functions
    const auto velocity = problem_.b(e, x);
    const auto c = problem_.c(e, x);

    y = - velocity.dot(gradphi) + c * phi;
    if (rescaling_)
      y /= velocity.two_norm();
  }

  //! get a reference to the GridView
  const typename Traits::GridViewType& getGridView() const
  { return dgfScalar_.getGridView(); }


private:
  const Problem& problem_;
  const bool rescaling_;

  mutable typename DGFScalar::Traits::RangeType phi;
  mutable typename DGFGradient::Traits::RangeType gradphi;
  DGFScalar dgfScalar_;
  DGFGradient dgfGradient_;
};

#endif  // DUNE_ULTRAWEAK_TRANSPORT_DISCRETEGRIDFUNCTIONTRANSPORTDIFFOP_HH
