#ifndef DUNE_ULTRAWEAK_TRANSPORT_NORMALEQLOCALOPERATOR_HH
#define DUNE_ULTRAWEAK_TRANSPORT_NORMALEQLOCALOPERATOR_HH

#include <dune/pdelab.hh>


template<typename Problem, typename FEM>
class SymmetricTransportOperator :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // define flags controlling global assembler
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };

  using LocalBasisType = typename
    FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using CacheType = Dune::PDELab::LocalBasisCache<LocalBasisType>;

SymmetricTransportOperator (const Problem& problem_, const int intorder_ = 2,
                            const bool rescaling=false, const double rescalingOut = 1.0) :
   problem(problem_), intorder(intorder_), rescaling_(rescaling), rescalingOut_(rescalingOut) {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu,
                     const X& x, const LFSV& lfsv,
                     R& r) const {
    using RF = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::Traits::SizeType;

    const auto& cell = eg.entity();
    const auto geo = eg.geometry();

    auto gradphi = makeJacobianContainer(lfsu);
    std::vector<RF> Bstarv(lfsu.size());

    for (const auto& qp : quadratureRule(geo, intorder)) {
      auto& phi = cache.evaluateFunction(qp.position(), lfsu.finiteElement().localBasis());
      auto& js = cache.evaluateJacobian(qp.position(), lfsu.finiteElement().localBasis());
      const auto S = geo.jacobianInverseTransposed(qp.position());

      for (size_type i=0; i<lfsu.size(); ++i)
        S.mv(js[i][0], gradphi[i][0]);

      // evaluate data functions
      const auto velocity = problem.b(cell, qp.position());
      const auto c = problem.c(cell, qp.position());

      // evaluate adjoint operator
      for (size_type i=0; i<lfsu.size(); i++)
        Bstarv[i] = -velocity.dot(gradphi[i][0]) + c*phi[i];

      RF Bstaru = 0.0;
      for (size_type i=0; i<lfsu.size(); i++)
        Bstaru += x(lfsu,i)*Bstarv[i];

      RF factor = qp.weight() * geo.integrationElement(qp.position());
      if (rescaling_)
        factor /= velocity.two_norm();
      for (size_type i=0; i<lfsu.size(); i++)
        r.accumulate(lfsu, i, factor * Bstaru * Bstarv[i]);
    }
  }

  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu,
                        const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    using RF = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::Traits::SizeType;

    // get Jacobian and determinant
    // assume the transformation is linear
    const auto& cell = eg.entity();
    const auto geo = eg.geometry();

    auto gradphi = makeJacobianContainer(lfsu);
    std::vector<RF> Bstarv(lfsu.size());

    for (const auto& qp : quadratureRule(geo, intorder)) {
      auto& phi = cache.evaluateFunction(qp.position(), lfsu.finiteElement().localBasis());
      auto& js = cache.evaluateJacobian(qp.position(), lfsu.finiteElement().localBasis());
      const auto S = geo.jacobianInverseTransposed(qp.position());

      for (size_type i=0; i<lfsu.size(); ++i)
        S.mv(js[i][0], gradphi[i][0]);

      // evaluate data functions
      const auto velocity = problem.b(cell, qp.position());
      const auto c = problem.c(cell, qp.position());

      // evaluate adjoint operator
      for (size_type i=0; i<lfsu.size(); i++)
        Bstarv[i] = -velocity.dot(gradphi[i][0]) + c*phi[i];

      RF factor = qp.weight() * geo.integrationElement(qp.position());
      if (rescaling_)
        factor /= velocity.two_norm();
      for (size_type j=0; j<lfsu.size(); j++)
        for (size_type i=0; i<lfsv.size(); i++)
          mat.accumulate(lfsv, i, lfsu, j, factor * Bstarv[i] * Bstarv[j]);
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s,
                      const LFSV& lfsv_s, R& r_s) const {
    using RF = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using RangeType = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
    using size_type = typename LFSU::Traits::SizeType;

    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = ig.geometry();
    const auto geoInInside = ig.geometryInInside();
    const auto insideCell = ig.inside();

    // do inflow boundary check at the midpoint
    const auto velocity = problem.b(insideCell, insideCell.geometry().local(geo.center()));

    auto normal = ig.centerUnitOuterNormal();

    RangeType u;
    std::vector<RangeType> phi(lfsu_s.size());

    // assemble u v on outflow boundary
    if (ig.boundary() and velocity.dot(normal) > -1e-12) {
      for (const auto& qp : quadratureRule(geo, intorder)) {
        const auto& insidepos = geoInInside.global(qp.position());

        // scalar functions
        phi = cache.evaluateFunction(insidepos, lfsu_s.finiteElement().localBasis());

        // evaluate u
        u = 0.0;
        for (size_type i=0; i<lfsu_s.size(); i++)
          u += x_s(lfsu_s,i) * phi[i];

        const RF factor = rescalingOut_ * qp.weight() * geo.integrationElement(qp.position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s, i, factor * u * phi[i]);
      }
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s,
                      const LFSV& lfsv_s, M& mat_ss) const {
    using RF = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using RangeType = typename
      LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
    using size_type = typename LFSU::Traits::SizeType;

    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = ig.geometry();
    const auto geoInInside = ig.geometryInInside();
    const auto insideCell = ig.inside();

    // do inflow boundary check at the midpoint
    const auto velocity = problem.b(insideCell, insideCell.geometry().local(geo.center()));

    std::vector<RangeType> phi(lfsu_s.size());

    // assemble u v on outflow boundary
    if (ig.boundary() and velocity.dot(ig.centerUnitOuterNormal()) > -1e-12) {
      for (const auto& qp : quadratureRule(geo, intorder)) {
        const auto& insidepos = geoInInside.global(qp.position());

        // scalar functions
        phi = cache.evaluateFunction(insidepos, lfsu_s.finiteElement().localBasis());

        const RF factor = rescalingOut_ * qp.weight() * geo.integrationElement(qp.position());
        for (size_type i=0; i<lfsu_s.size(); i++)
          for (size_type j=0; j<lfsv_s.size(); j++)
            mat_ss.accumulate(lfsu_s, i, lfsv_s, j, factor * phi[i] * phi[j]);
      }
    }
  }

  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const {
    using RF = typename
      LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSV::Traits::SizeType;

    const auto& cell = eg.entity();
    const auto geo = eg.geometry();

    for (const auto& qp : quadratureRule(geo, intorder)) {
      const auto fval = problem.f(cell, qp.position());
      auto& psi = cache.evaluateFunction(qp.position(), lfsv.finiteElement().localBasis());
      const RF integrationFactor = qp.weight() * geo.integrationElement(qp.position());
      for (size_type i=0; i<lfsv.size(); i++)
        r.accumulate(lfsv, i, integrationFactor * -fval * psi[i]);
    }
  }

  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv_s, R& r_s) const {
    using RF = typename
      LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSV::Traits::SizeType;

    const auto& interface = ig.intersection();
    const auto geo = ig.geometry();

    for (const auto& qp : quadratureRule(geo, intorder)) {
      auto bctype = problem.bctype(ig.intersection(), qp.position());
      if (ig.boundary() and bctype == Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet) {
        const auto& insidepos = ig.geometryInInside().global(qp.position());
        auto& phi = cache.evaluateFunction(insidepos, lfsv_s.finiteElement().localBasis());

        const auto gval = problem.g(interface, qp.position());
        auto velocity = problem.b(ig.inside(), insidepos);

        using std::abs;
        const auto scalingFactor = abs(velocity * ig.centerUnitOuterNormal());

        RF integrationFactor = qp.weight() * geo.integrationElement(qp.position());
        for (size_type i=0; i<lfsv_s.size(); i++)
          r_s.accumulate(lfsv_s, i, integrationFactor * -gval * phi[i] * scalingFactor);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& z, const LFSV& lfsv,
                              R& r) const {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
  void jacobian_apply_boundary(const IG& ig, const LFSU& lfsu_s, const X& z_s,
                               const LFSV& lfsv_s, Y& y_s) const {
    alpha_boundary(ig,lfsu_s,z_s,lfsv_s,y_s);
  }

private:
  const Problem& problem;
  int intorder;
  const bool rescaling_;
  const double rescalingOut_;

  CacheType cache;
};

#endif // DUNE_ULTRAWEAK_TRANSPORT_NORMALEQLOCALOPERATOR_HH
