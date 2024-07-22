//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXAMPLE_BELOS_SOLVE
#define KOKKOS_EXAMPLE_BELOS_SOLVE

// Tpetra
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "CGSolve.hpp" // for result_struct

#include "TrilinosCouplings_config.h"
#if defined(HAVE_TRILINOSCOUPLINGS_BELOS) && defined(HAVE_TRILINOSCOUPLINGS_MUELU)

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "SGPreconditioner.hpp"
#include "MeanBasedPreconditioner.hpp"
#include "MueLuPreconditioner.hpp"

namespace Kokkos {
namespace Example {

// Functor to copy the nodal coordinates out of the mesh fixture into
// a multi-vector.  This functor is necessary because the respective views
// have potentially different layouts.
template <typename coords_type, typename coords_vec_type>
struct FillCoords {
  typedef typename coords_type::execution_space execution_space;
  typedef typename coords_type::size_type size_type;

  const coords_type m_coords;
  const coords_vec_type m_coords_vec;
  const size_type m_dim;

  FillCoords( const coords_type& coords, const coords_vec_type& coords_vec )
    : m_coords(coords), m_coords_vec(coords_vec), m_dim(coords.extent(1))
  {
    // Note:  coords contains off-processor halo nodes and thus is longer
    // than coords_vec, which is the same length as the solution vector.
    // These extra halo nodes are stored at the end.
    Kokkos::parallel_for( m_coords_vec.extent(0), *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const
  {
    for (size_type j=0; j<m_dim; ++j)
      m_coords_vec(i,j) = m_coords(i,j);
  }
};

template <typename coords_type, typename coords_vec_type>
void fill_coords( const coords_type& coords,
                  const coords_vec_type& coords_vec ) {
  FillCoords<coords_type,coords_vec_type>(coords, coords_vec);
}

template <typename S, typename V, typename O>
struct ExtractEnsembleIts {
  static std::vector<int>
  apply(const Belos::SolverManager<S,V,O>& solver) {
    return std::vector<int>();
  }
};

template <class SM, class SV, class LO, class GO, class N, class Mesh>
result_struct
belos_solve(
  Tpetra::CrsMatrix<SM,LO,GO,N>& A,
  const Tpetra::MultiVector<SV,LO,GO,N>& b,
  Tpetra::MultiVector<SV,LO,GO,N>& x,
  Teuchos::RCP< Tpetra::Operator<SM,LO,GO,N> >& precOp,
  const Mesh& mesh,
  const int use_muelu,
  const int use_mean_based,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const unsigned max_iter = 200,
  const typename Kokkos::ArithTraits<SV>::mag_type tolerance =
    Kokkos::ArithTraits<SV>::epsilon())
{
  typedef Tpetra::Operator<SM,LO,GO,N> OperatorType;
  typedef Tpetra::MultiVector<SV,LO,GO,N> VectorType;
  typedef typename VectorType::dot_type BelosScalarType;
  typedef Belos::LinearProblem<BelosScalarType, VectorType, OperatorType> ProblemType;
  typedef Belos::PseudoBlockCGSolMgr<BelosScalarType, VectorType, OperatorType> CGSolverType;
  typedef Belos::PseudoBlockGmresSolMgr<BelosScalarType, VectorType, OperatorType> GmresSolverType;
  typedef Belos::SolverManager<BelosScalarType, VectorType, OperatorType> SolverType;
  typedef SGPreconditioner<SM,LO,GO,N> PreconditionerType;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::ParameterList;

  // Create some timers used by Belos so we can access them
  Teuchos::RCP<Teuchos::Time> time_mat_vec =
    Teuchos::TimeMonitor::getNewTimer("Belos: Operation Op*x");
  Teuchos::RCP<Teuchos::Time> time_prec_apply =
    Teuchos::TimeMonitor::getNewTimer("Belos: Operation Prec*x");
  Teuchos::RCP<Teuchos::Time> time_total;
  Teuchos::RCP<Teuchos::Time> time_prec_setup =
    Teuchos::TimeMonitor::getNewTimer("Total MueLu setup time");

  //--------------------------------
  // Create preconditioner if requested and we weren't given one
  if (use_muelu && precOp == Teuchos::null) {
    Teuchos::TimeMonitor timeMon(*time_prec_setup);
    RCP<PreconditionerType> preconditioner;

    // Create tpetra-vector storing coordinates for repartitioning
    typename Mesh::node_coord_type node_coords = mesh.node_coord();
    //Teuchos::RCP<VectorType> coords =
    Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> > coords =
      Teuchos::rcp(new Tpetra::MultiVector<double,LO,GO,N>(x.getMap(), node_coords.extent(1)));
    fill_coords(node_coords, coords->getLocalViewDevice(Tpetra::Access::ReadWrite));

    RCP<ParameterList> mueluParams = Teuchos::sublist(fenlParams, "MueLu");
    if (use_mean_based) {
      precOp = build_mean_based_muelu_preconditioner(rcpFromRef(A), mueluParams,
                                                     coords);
    }
    else {
      precOp = build_muelu_preconditioner(rcpFromRef(A), mueluParams, coords);
    }
  }

  //--------------------------------
  // Set up linear solver
  RCP<ParameterList> belosParams = Teuchos::sublist(fenlParams, "Belos");

  if (!(belosParams->isParameter("Convergence Tolerance")))
    belosParams->set("Convergence Tolerance", tolerance);
  if (!(belosParams->isParameter("Maximum Iterations")))
    belosParams->set("Maximum Iterations", Teuchos::as<int>(max_iter));
  if (!(belosParams->isParameter("Output Frequency")))
    belosParams->set("Output Frequency", 1);

  RCP<ProblemType> problem =
    rcp(new ProblemType(rcpFromRef(A), rcpFromRef(x), rcpFromRef(b)));

  std::string belos_solver = belosParams->get("Belos Solver", "CG");
  RCP<SolverType> solver;
  if (belos_solver == "CG" ||
      belos_solver == "cg") {
    time_total =
      Teuchos::TimeMonitor::getNewTimer("Belos: PseudoBlockCGSolMgr total solve time");
    solver = rcp(new CGSolverType(problem, belosParams));
  }
  else if (belos_solver == "GMRES" ||
           belos_solver == "Gmres" ||
           belos_solver == "gmres") {
    time_total =
      Teuchos::TimeMonitor::getNewTimer("Belos: PseudoBlockGmresSolMgr total solve time");
    solver = rcp(new GmresSolverType(problem, belosParams));
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "Invalid solver " << belos_solver);

  if (use_muelu || use_mean_based){
     problem->setRightPrec(precOp);
  }
  const bool isSet = problem->setProblem();

  TEUCHOS_TEST_FOR_EXCEPTION(!isSet, std::runtime_error,
                             "Belos failed to set problem correctly.");

  //--------------------------------
  // Solve for nonlinear update
  Belos::ReturnType result = solver->solve();
  bool converged = (result == Belos::Converged);
  TEUCHOS_TEST_FOR_EXCEPTION(!converged, std::runtime_error,
                             "Belos solver did not converge!");

  result_struct cgsolve;
  cgsolve.iteration = solver->getNumIters();
  //cgsolve.norm_res = solver->achievedTol();

  cgsolve.iter_time = time_total->totalElapsedTime();
  cgsolve.total_time = time_total->totalElapsedTime();
  cgsolve.matvec_time = time_mat_vec->totalElapsedTime();
  cgsolve.prec_apply_time = time_prec_apply->totalElapsedTime();
  cgsolve.prec_setup_time = time_prec_setup->totalElapsedTime();

  // Extract the number of iterations for ensembles
  cgsolve.ensemble_its =
    ExtractEnsembleIts<BelosScalarType, VectorType, OperatorType>::apply(*solver);
  // if (cgsolve.ensemble_its.size() > 0) {
  //   std::cout << std::endl << "ensemble iterations = ";
  //   for (std::size_t i=0; i<cgsolve.ensemble_its.size(); ++i)
  //     std::cout << cgsolve.ensemble_its[i] << " ";
  //   std::cout << std::endl;
  // }

  return cgsolve;
}


} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else

namespace Kokkos {
namespace Example {

template <class SM, class SV, class LO, class GO, class N, class Mesh>
result_struct
belos_solve(
  Tpetra::CrsMatrix<SM,LO,GO,N>& A,
  const Tpetra::MultiVector<SV,LO,GO,N>& b,
  Tpetra::MultiVector<SV,LO,GO,N>& x,
  Teuchos::RCP< Tpetra::Operator<SM,LO,GO,N> >& precOp,
  const Mesh& mesh,
  const int use_muelu,
  const int use_mean_based,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const unsigned max_iter = 200,
  const typename Kokkos::ArithTraits<SV>::mag_type tolerance =
    Kokkos::ArithTraits<SV>::epsilon())
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "Belos solver requested but not compiled!");
}

} // namespace Example
} // namespace Kokkos

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_BELOS_SOLVE */
