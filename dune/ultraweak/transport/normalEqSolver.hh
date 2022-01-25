// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_TRANSPORT_NORMAL_EQ_SOLVER_HH
#define DUNE_ULTRAWEAK_TRANSPORT_NORMAL_EQ_SOLVER_HH

#include <dune/pdelab.hh>

#include <dune/istl/solverfactory.hh>

#include "discreteGridFunctionTransportDiffOp.hh"
#include "normalEqLocalOperator.hh"

#include "../constraints.hh"
#include "../fullBlockInvert.hh"
#include "../reconstruction.hh"
#include "../schurcomplement.hh"
#include "../solvingManager.hh"

template<typename GV, std::size_t order_>
struct NormalEqTraits {
  using GridView = GV;
  static const int dim = GV::dimension;
  static const int order = order_;
  using DF = typename GV::Grid::ctype; // type for ccordinates
  using RF = double; //Dune::Float128;           // type for computations

  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,order>;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, RF>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
};

// defines reconstruction traits based on normal eq traits
template<typename T>
struct DefaultReconstructionTraits {
  using GridView = typename T::GridView;
  static const int dim = T::dim;
  static const int order = T::order;
  using DF = typename T::DF;
  using RF = typename T::RF;

  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,order, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange>;
  using NOCON = Dune::PDELab::NoConstraints;
  using BlockedVBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
  using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GridView,FEM,NOCON,BlockedVBE>;
  using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, RF>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
};

// Wrapper class for all solving aspects
template< std::size_t order, typename GV, typename Problem>
class NormalEqSolver : public SolvingManager<GV, Problem, NormalEqTraits<GV,order>> {
public:
  using Traits = NormalEqTraits<GV,order>;
  using ReconstructionTraits = DefaultReconstructionTraits<Traits>;

private:
  using Base = SolvingManager<GV, Problem, Traits>;

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using NormalEqOperatorType = SymmetricTransportOperator<Problem, typename Traits::FEM>;
  using NormalEqGO = Dune::PDELab::GridOperator<
    typename Traits::GridFunctionSpace, typename Traits::GridFunctionSpace,
    NormalEqOperatorType, MBE,
    typename Traits::RF, typename Traits::RF, typename Traits::RF,
    typename Base::CC, typename Base::CC>;

private:
  using SchurMatType = Dune::PDELab::Backend::Native<typename NormalEqGO::Jacobian>;
  using Y = Dune::PDELab::Backend::Native<typename Traits::CoefficientVector>;
  using SchurOperatorType = Dune::MatrixAdapter<SchurMatType, Y,Y>;

public:
  NormalEqSolver(const GV& gv, Problem& problem, Dune::ParameterTree& pTree) :
    Base(gv, typename Traits::FEM(gv), problem, pTree),
    reconstructionFem(), reconstructionGfs(gv, reconstructionFem),
    reconstructionCoeffs(reconstructionGfs)
  {
    gfs.name("Numerical solution of the normal equations");
    reconstructionGfs.name("DG");

    // Make a global operator
    const int dim = Traits::dim;
    MBE mbe(1<<(dim+1)); // guess nonzeros per row

    const int intorder = Traits::order + 2;
    lop_ = std::make_shared<NormalEqOperatorType>(problem, intorder, pTree.template get<bool>("rescaling"));

    normaleqgo = std::make_shared<NormalEqGO>(gfs, cc, gfs, cc, *lop_, mbe);

    // initialize solvers
    Dune::initSolverFactories<SchurOperatorType>();
  }

  Dune::InverseOperatorResult solvingStats;

  /**
     Solves linear transport using the 'normal equations'
  */
  void solve(Dune::ParameterTree& solverConfig) {
    // assemble the system matrix
    typename NormalEqGO::Jacobian wrappedSchurMat(*normaleqgo, 0.0);
    typename NormalEqGO::Domain tempEvalPoint(gfs, 0.0);
    normaleqgo->jacobian(tempEvalPoint, wrappedSchurMat);

    // assemble the rhs
    typename NormalEqGO::Range rhs(gfs, 0.0);
    normaleqgo->residual(tempEvalPoint, rhs);
    rhs *= -1;
    auto natRHS = Dune::PDELab::Backend::native(rhs);

    // create linear operator
    SchurMatType schurMat = Dune::PDELab::Backend::native(wrappedSchurMat);
    const auto schurop = std::make_shared<SchurOperatorType>(schurMat);

    // create solver
    const auto solver = getSolverFromFactory(schurop, solverConfig);

    // solve
    auto solution = std::make_shared<Y>(schurMat.N());
    *solution = 0.0;
    solver->apply(*solution, natRHS, solvingStats);
    coeffs.attach(solution); // attach the underlying container

    if (solvingStats.condition_estimate > 0)
      std::cout << "Convergence rate: " << solvingStats.conv_rate << ", condition (estimate): " << solvingStats.condition_estimate << std::endl;

    writeToVTK();
  }

  void writeToVTK() const {
    // write output files
    using GFS = typename Traits::GridFunctionSpace;
    using CoefficientVector = typename Traits::CoefficientVector;

    using DGF = typename Traits::DiscreteGridFunction;
    using VTKF = Dune::PDELab::VTKGridFunctionAdapter<DGF>;

    const int refinements = pTree.template get<int>("visualization.refinements");
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals(refinements), false);
    DGF dgf = this->getDiscreteGridFunction();
    vtkwriter.addVertexData(std::make_shared<VTKF>(dgf, "test_sol"));

    using RecDGF = DiscreteGridFunctionTransportDiffOp<GFS, CoefficientVector, Problem>;
    using VTKF_REC = Dune::PDELab::VTKGridFunctionAdapter<RecDGF>;

    auto recdgf = getDiscreteGridFunctionReconstruction();
    vtkwriter.addCellData(std::make_shared<VTKF_REC>(recdgf, "reconstruction_subsampling"));

    std::string filename = "normal_eq_sol";
    if (Traits::order == 2)
      filename += "_so";
    vtkwriter.write(filename, Dune::VTK::appendedraw);
  }

  auto getDiscreteGridFunctionReconstruction() const {
    using GFS = typename Traits::GridFunctionSpace;
    using Coeffs = typename Traits::CoefficientVector;
    return DiscreteGridFunctionTransportDiffOp<GFS, Coeffs, Problem>(
             gfs, coeffs, problem, pTree.template get<bool>("rescaling"));
  }

  const auto& getReconstructionGfs() const {
    return reconstructionGfs;
  }

public:
  using Base::problem;
  using Base::gv;
  using Base::gfs;
  using Base::cc;
  using Base::coeffs;
  using Base::pTree;

  typename ReconstructionTraits::FEM reconstructionFem;
  typename ReconstructionTraits::GridFunctionSpace reconstructionGfs;
  typename ReconstructionTraits::CoefficientVector reconstructionCoeffs;

private:
  std::shared_ptr<NormalEqOperatorType> lop_;
  std::shared_ptr<NormalEqGO> normaleqgo;
};

#endif  // DUNE_ULTRAWEAK_TRANSPORT_NORMAL_EQ_SOLVER_HH
