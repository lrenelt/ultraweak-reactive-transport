// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_DGSOLVER_HH
#define DUNE_ULTRAWEAK_DGSOLVER_HH

#include "dune/pdelab.hh"

#include "solvingManager.hh"

template<class GV, std::size_t order_>
struct DGSolverTraits {
  static const int dim = GV::dimension;
  using DF = typename GV::Grid::ctype;
  using RF = double;
  static const int order = order_;

  using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<DF, RF, order, dim>;
  using Constraints = Dune::PDELab::NoConstraints;
  using VectorBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed, Dune::QkStuff::QkSize<order, dim>::value>;
  using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<GV, FEM, Constraints, VectorBackend>;
  using CoefficientVector = Dune::PDELab::Backend::Vector<GridFunctionSpace, DF>;
  using DiscreteGridFunction = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
};

// Wrapper class for all solving aspects
template<std::size_t order, typename GV, typename Problem, typename Traits = DGSolverTraits<GV,order>>
class DGSolver : public SolvingManager<GV, Problem, Traits> {
private:
  using Base = SolvingManager<GV, Problem, Traits>;

public:
  DGSolver(const GV& gv_, Problem& problem_, Dune::ParameterTree& pTree_) : Base(gv_, typename Traits::FEM(), problem_, pTree_) {
    gfs.name("dg_sol");
  }

  void solve() {
    const int dim = Traits::dim;
    using DF = typename Traits::DF;
    using RF = typename Traits::RF;

    // Local operator
    using LocalOperator = Dune::PDELab::ConvectionDiffusionDG<Problem, typename Traits::FEM>;
    LocalOperator localOperator(problem);

    // constraints
    using GridFunctionSpace = typename Traits::GridFunctionSpace;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<RF>::Type;
    ConstraintsContainer constraintsContainer;
    constraintsContainer.clear();
    Dune::PDELab::constraints(gfs, constraintsContainer);

    // Grid operator
    gfs.update();
    using std::pow;
    const int dofestimate = pow(2, dim) * gfs.maxLocalSize();
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MatrixBackend matrixBackend(dofestimate);
    using GridOperator = Dune::PDELab::GridOperator<GridFunctionSpace,
                                                    GridFunctionSpace,
                                                    LocalOperator,
                                                    MatrixBackend,
                                                    DF,
                                                    RF,
                                                    RF,
                                                    ConstraintsContainer,
                                                    ConstraintsContainer>;
    GridOperator gridOperator(gfs, constraintsContainer,
                              gfs, constraintsContainer,
                              localOperator, matrixBackend);

    std::cout << "gfs with " << gfs.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << constraintsContainer.size() << " dofs generated  "<< std::endl;

    // Solution vector
    coeffs = 0.0;

    // Solve matrix-based with PDELab preset
    const int verbosity = 0;
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver linearSolver(verbosity);
    using CoefficientVector = typename Traits::CoefficientVector;
    using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, CoefficientVector>;

    // create grid function
    using DGF = Dune::PDELab::DiscreteGridFunction<GridFunctionSpace, CoefficientVector>;
    DGF dgf(gfs, coeffs);

    const double reduction = pTree.template get<double>("reduction");
    Solver solver(gridOperator, linearSolver, coeffs, reduction);
    solver.apply();
    writeVTK();
    writeVelocityFieldVTK();
  }

  // Visualization
  void writeVTK() {
    using VTKWriter = Dune::SubsamplingVTKWriter<GV>;
    Dune::RefinementIntervals subint(1);
    VTKWriter vtkwriter(gv, subint);
    std::string vtkfile("dg_ref_sol");
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, coeffs,
                                         Dune::PDELab::vtk::defaultNameScheme());
    vtkwriter.write(vtkfile, Dune::VTK::ascii);
  }

  void writeVelocityFieldVTK() {
    auto blambda = [this](const auto& el, const auto& x){ return this->problem.b(el, x); };
    auto bGridFunction = Dune::PDELab::makeGridFunctionFromCallable(gv,blambda);

    using DF = typename GV::Grid::ctype;
    using RF = double;
    using FEM = Dune::PDELab::RaviartThomasLocalFiniteElementMap<GV,DF,RF,0>;
    using CON = Dune::PDELab::NoConstraints;
    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    using VEC = Dune::PDELab::Backend::Vector<GFS,DF>;

    FEM fem(gv);
    GFS gfs(gv, fem);
    gfs.name("interpolation");
    VEC coefficientVector(gfs);
    Dune::PDELab::interpolate(bGridFunction, gfs, coefficientVector);

    // set up vtk writer
    using VTKWriter = Dune::SubsamplingVTKWriter<GV>;
    Dune::RefinementIntervals subint(1);
    VTKWriter vtkwriter(gv, subint);
    std::string vtkfile("velocity_field");
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, coefficientVector, Dune::PDELab::vtk::defaultNameScheme());

    // write output
    vtkwriter.write(vtkfile, Dune::VTK::ascii);
  }

protected:
  using Base::problem;
  using Base::gv;
  using Base::gfs;
  using Base::coeffs;
  using Base::pTree;
};

#endif // DUNE_ULTRAWEAK_DGSOLVER_HH
