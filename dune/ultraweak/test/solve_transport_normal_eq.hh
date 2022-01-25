// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ULTRAWEAK_TEST_SOLVE_TRANSPORT_NORMAL_EQ_HH
#define DUNE_ULTRAWEAK_TEST_SOLVE_TRANSPORT_NORMAL_EQ_HH

#include "problems.hh"
#include "../darcy-velocity/catalysatorProblems.hh"
#include "../darcy-velocity/darcyMixedSolver.hh"
#include "../darcy-velocity/darcyProblem.hh"
#include "../transport/normalEqSolver.hh"

void solveTransportNormalEq(Dune::ParameterTree& pTree) {
    const int dim = 2;
    using RF = double;
    static const std::size_t order = 1;

    // make grid
    const int refinement = pTree.get<int>("grid.refinement");

    typedef Dune::YaspGrid<dim> Grid;
    Dune::FieldVector<double, dim> domain(1.0);
    std::array<int, dim> domainDims;

    domain[0] = 1.0;
    domain[1] = 1.0;
    domainDims[0] = pTree.get<int>("grid.yasp_x");
    domainDims[1] = pTree.get<int>("grid.yasp_y");

    auto gridp = std::make_shared<Grid>(domain, domainDims);
    gridp->globalRefine(refinement);
    auto gv = gridp->leafGridView();

    Dune::ParameterTree solverpTree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("solver_config.ini", solverpTree);

    using GV = typename Grid::LeafGridView;

    PureTransportProblem<GV,RF> baseProblem(pTree);
    if (pTree.get<bool>("darcy.useDarcy")) {
      using DarcyProblem = CatalysatorProblem3<GV,RF>;
      DarcyProblem darcyProblem(pTree);
      DarcyMixedSolver darcysolver(gv, darcyProblem, pTree);
      darcysolver.solve();
      auto darcy_velocity = darcysolver.getDiscreteGridFunction();

      DarcyVelocityAdapter problem(baseProblem, darcy_velocity);
      using Solver = NormalEqSolver<GV,decltype(problem),order>;
      Solver solver(gv, problem, pTree);
      solver.solve(solverpTree);
    }
    else {
      using Solver = NormalEqSolver<GV,decltype(baseProblem),order>;
      Solver solver(gv, baseProblem, pTree);
      solver.solve(solverpTree);
    }
}

void solveTransportNormalEqSecondOrder(Dune::ParameterTree& pTree) {
    const int dim = 2;
    using RF = double;
    static const std::size_t order = 2;

    // make grid
    const int refinement = pTree.get<int>("grid.refinement");

    typedef Dune::YaspGrid<dim> Grid;
    Dune::FieldVector<double, dim> domain(1.0);
    std::array<int, dim> domainDims;

    domain[0] = 1.0;
    domain[1] = 1.0;
    domainDims[0] = pTree.get<int>("grid.yasp_x");
    domainDims[1] = pTree.get<int>("grid.yasp_y");

    auto gridp = std::make_shared<Grid>(domain, domainDims);
    gridp->globalRefine(refinement);
    auto gv = gridp->leafGridView();

    Dune::ParameterTree solverpTree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("solver_config.ini", solverpTree);

    using GV = typename Grid::LeafGridView;

    PureTransportProblem<GV,RF> baseProblem(pTree);
    if (pTree.get<bool>("darcy.useDarcy")) {
      using DarcyProblem = CatalysatorProblem3<GV,RF>;
      DarcyProblem darcyProblem(pTree);
      DarcyMixedSolver darcysolver(gv, darcyProblem, pTree);
      darcysolver.solve();
      auto darcy_velocity = darcysolver.getDiscreteGridFunction();

      DarcyVelocityAdapter problem(baseProblem, darcy_velocity);
      using Solver = NormalEqSolver<GV,decltype(problem),order>;
      Solver solver(gv, problem, pTree);
      solver.solve(solverpTree);
    }
    else {
      using Solver = NormalEqSolver<GV,decltype(baseProblem),order>;
      Solver solver(gv, baseProblem, pTree);
      solver.solve(solverpTree);
    }
}

#endif  // DUNE_ULTRAWEAK_TEST_SOLVE_TRANSPORT_NORMAL_EQ_HH
