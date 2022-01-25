// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ULTRAWEAK_TESTUTIL_HH
#define DUNE_ULTRAWEAK_TESTUTIL_HH

#include <iostream>
#include <fstream>

#include <dune/common/parametertreeparser.hh>

#include <dune/pdelab/test/l2norm.hh>

#include "dgsolver.hh"
#include "darcy-velocity/catalysatorProblems.hh"
#include "darcy-velocity/darcyProblem.hh"
#include "darcy-velocity/darcyMixedSolver.hh"
#include "refinementadapter.hh"
#include "test/problems.hh"
#include "transport/normalEqSolver.hh"

/**
   Write computed data to a csv file
*/
// TODO: pass gridwidths (how to calculate?)
void writeToCSV(std::string filename, const std::vector<double>& data, const std::vector<typename Dune::InverseOperatorResult> res, Dune::ParameterTree& pTree, double h0, double refMin) {
  // open file
  std::ofstream file;
  filename += ".csv";
  file.open(filename);

  // write header
  file << "gridwidth, l2error, condition_estimate, iterations, time \n";

  // write data
  for (size_t i=0; i<data.size(); i++) {
    file << std::to_string(pow(0.5, refMin+i) * h0) << ", " << std::to_string(data[i]) << ", ";
    file << res[i].condition_estimate << ", " << res[i].iterations << ", " << res[i].elapsed <<  "\n";
  }

  // close the file
  file.close();
}


template<typename Grid>
void doTest(std::shared_ptr<Grid> gridp, Dune::ParameterTree pTree) {
  std::size_t ref_min, ref_max;
  ref_min = pTree.get<int>("grid.refMin");
  ref_max = pTree.get<int>("grid.refMax");

  gridp->globalRefine(ref_max);

  for (std::size_t ref = ref_min; ref<=ref_max; ++ref)
    solve_normal_eq(gridp->levelGridView(ref), pTree);
}

// TODO: pass reference function and numerical solve function
template<std::size_t order, typename Grid, typename Problem>
void doTestWithRefSol(std::shared_ptr<Grid> gridp, Dune::ParameterTree pTree,
                      Dune::ParameterTree pTreeConvTest,
                      Problem& problem) {
  using GV = typename Grid::LevelGridView;
  using Solver = NormalEqSolver<order, GV, Problem>;
  using RF = typename Solver::Traits::RF;

  // read convergence test specific parameters
  std::size_t refMin = pTreeConvTest.get<int>("grid.refMin");
  std::size_t refMax = pTreeConvTest.get<int>("grid.refMax");
  std::string filename = pTreeConvTest.get<std::string>("filename");

  gridp->globalRefine(refMax);

  // get finest refinement level grid view
  auto gv = gridp->levelGridView(refMax);

  // calculate DG-reference solution
  std::cout << "\n\nSolving for DG-reference solution..." << std::endl;;
  using DGSolver = DGSolver<order,GV,Problem>;
  DGSolver dgsolver(gv, problem, pTree);
  dgsolver.solve();
  auto refSol = dgsolver.getDiscreteGridFunction();
  std::cout << "Solving for DG-reference solution...done!" << std::endl;

  // initialize QoI vectors
  std::vector<double> l2errors;
  std::vector<typename Dune::InverseOperatorResult> solvingStats;

  Dune::ParameterTree solverpTree;
  Dune::ParameterTreeParser ptreeparser;
  ptreeparser.readINITree("solver_config.ini", solverpTree);

  // loop over all refined grids
  for (std::size_t ref = refMin; ref<=refMax; ++ref) {
    std::cout << "\n\nRefinement stage " << ref << "/" << refMax << " ..." << std::endl;
    auto gvRef = gridp->levelGridView(ref);

    // solve the normal equations
    Solver normaleqsolver(gvRef, problem, pTree);
    normaleqsolver.solve(solverpTree);
    auto numSol = normaleqsolver.getDiscreteGridFunctionReconstruction();
    solvingStats.push_back(normaleqsolver.solvingStats);

    // adapt solution to finest grid level and compute the error
    auto numSolRef = RefinementAdapter(numSol, gv, refMax-ref);
    auto diff = Dune::PDELab::DifferenceAdapter(refSol, numSolRef);
    const std::size_t intorder = order+1;
    l2errors.push_back(l2norm(diff,intorder));
  }

  // write output data
  double h0 = 1. / pTreeConvTest.get<int>("grid.yasp_x"); // TODO: remove
  if (order == 2)
    filename += "_so";
  writeToCSV(filename, l2errors, solvingStats, pTree, h0, refMin);
}

template<std::size_t order, typename Grid, typename Problem>
void doTestWithRefSolDarcy(std::shared_ptr<Grid> gridp, Dune::ParameterTree pTree,
                      Dune::ParameterTree pTreeConvTest,
                      Problem& baseProblem) {
  using GV = typename Grid::LevelGridView;
  using RF = double;

  // read convergence test specific parameters
  std::size_t refMin = pTreeConvTest.get<int>("grid.refMin");
  std::size_t refMax = pTreeConvTest.get<int>("grid.refMax");
  std::string filename = pTreeConvTest.get<std::string>("filename");

  gridp->globalRefine(refMax);

  // get finest refinement level grid view
  auto gv = gridp->levelGridView(refMax);

  // calculate Darcy-velocity field
  std::cout << "\n\nSolving for fine Darcy field..." << std::endl;;
  using DarcyProblem = CatalysatorProblem3<GV, RF>;

  DarcyProblem darcyProblem(pTree);
  DarcyMixedSolver darcysolver(gv, darcyProblem, pTree);
  darcysolver.solve();
  auto darcyVelocity = darcysolver.getDiscreteGridFunction();
  std::cout << "Solving for fine Darcy field...done." << std::endl;;

  DarcyVelocityAdapter problem(baseProblem, darcyVelocity);

  // calculate DG-reference solution
  std::cout << "\n\nSolving for DG-reference solution..." << std::endl;;
  static const std::size_t orderdg = 1;
  using DGSolver = DGSolver<orderdg,GV,decltype(problem)>;
  DGSolver dgsolver(gv, problem, pTree);
  dgsolver.solve();
  auto refSol = dgsolver.getDiscreteGridFunction();
  std::cout << "Solving for DG-reference solution...done!" << std::endl;

  // initialize QoI vectors
  std::vector<double> l2errors;
  std::vector<typename Dune::InverseOperatorResult> solvingStats;

  Dune::ParameterTree solverpTree;
  Dune::ParameterTreeParser ptreeparser;
  ptreeparser.readINITree("solver_config.ini", solverpTree);

  // loop over all refined grids
  using DGFType = DarcyMixedSolver<GV, DarcyProblem>::DiscreteGridFunction;
  using Solver = NormalEqSolver<order,GV, DarcyVelocityAdapter<Problem,DGFType>>;
  for (std::size_t ref = refMin; ref<=refMax; ++ref) {
    std::cout << "\n\nRefinement stage " << ref << "/" << refMax << " ..." << std::endl;
    auto gvRef = gridp->levelGridView(ref);

    // solve Darcy again on the coarse grid
    // TODO: reuse fine solution?
    DarcyMixedSolver darcysolver(gvRef, darcyProblem, pTree);
    darcysolver.solve();
    const auto darcyVelocityRef = darcysolver.getDiscreteGridFunction();
    std::cout << "Solved Darcy" << std::endl;
    DarcyVelocityAdapter problemRef(baseProblem, darcyVelocityRef);

    // solve the normal equations
    Solver normaleqsolver(gvRef, problemRef, pTree);
    normaleqsolver.solve(solverpTree);
    auto numSol = normaleqsolver.getDiscreteGridFunctionReconstruction();
    solvingStats.push_back(normaleqsolver.solvingStats);

    // adapt solution to finest grid level and compute the error
    auto numSolRef = RefinementAdapter(numSol, gv, refMax-ref);
    auto diff = Dune::PDELab::DifferenceAdapter(refSol, numSolRef);
    const std::size_t intorder = order+1;
    l2errors.push_back(l2norm(diff,intorder));
  }

  // write output data
  double h0 = 1. / pTreeConvTest.get<int>("grid.yasp_x"); // TODO: remove
  if (order == 2)
    filename += "_so";
  writeToCSV(filename, l2errors, solvingStats, pTree, h0, refMin);
}

#endif  // DUNE_ULTRAWEAK_TESTUTIL_HH
