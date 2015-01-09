/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
  typedef typename coords_type::device_type device_type;
  typedef typename coords_type::size_type size_type;

  const coords_type m_coords;
  const coords_vec_type m_coords_vec;
  const size_type m_dim;

  FillCoords( const coords_type& coords, const coords_vec_type& coords_vec )
    : m_coords(coords), m_coords_vec(coords_vec), m_dim(coords.dimension_1())
  {
    // Note:  coords contains off-processor halo nodes and thus is longer
    // than coords_vec, which is the same length as the solution vector.
    // These extra halo nodes are stored at the end.
    Kokkos::parallel_for( m_coords_vec.dimension_0(), *this );
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

template <class SM, class SV, class LO, class GO, class N, class Mesh>
result_struct
belos_solve(
  const Teuchos::RCP<Tpetra::CrsMatrix<SM,LO,GO,N> >& A,
  const Teuchos::RCP<Tpetra::Vector<SV,LO,GO,N> >& b,
  const Teuchos::RCP<Tpetra::Vector<SV,LO,GO,N> >& x,
  const Mesh& mesh,
  const int use_muelu,
  const int use_mean_based,
  const unsigned max_iter = 200,
  const typename Kokkos::Details::ArithTraits<SV>::mag_type tolerance =
    Kokkos::Details::ArithTraits<SV>::epsilon())
{
  typedef Tpetra::Operator<SM,LO,GO,N> OperatorType;
  typedef Tpetra::MultiVector<SV,LO,GO,N> VectorType;
  typedef typename VectorType::dot_type BelosScalarType;
  typedef Belos::LinearProblem<BelosScalarType, VectorType, OperatorType> ProblemType;
  typedef Belos::PseudoBlockCGSolMgr<BelosScalarType, VectorType, OperatorType> SolverType;
  typedef SGPreconditioner<SM,LO,GO,N> PreconditionerType;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Create some timers used by Belos so we can access them
  Teuchos::RCP<Teuchos::Time> time_mat_vec =
    Teuchos::TimeMonitor::getNewTimer("Belos: Operation Op*x");
  Teuchos::RCP<Teuchos::Time> time_prec_apply =
    Teuchos::TimeMonitor::getNewTimer("Belos: Operation Prec*x");
  Teuchos::RCP<Teuchos::Time> time_total =
    Teuchos::TimeMonitor::getNewTimer("Belos: PseudoBlockCGSolMgr total solve time");
  Teuchos::RCP<Teuchos::Time> time_prec_setup =
    Teuchos::TimeMonitor::getNewTimer("Total MueLu setup time");

  //--------------------------------
  // Create preconditioner
  RCP<PreconditionerType> preconditioner;
  RCP<OperatorType> precOp;

  if (use_muelu) {
    Teuchos::TimeMonitor timeMon(*time_prec_setup);

    // Create tpetra-vector storing coordinates for repartitioning
    typename Mesh::node_coord_type node_coords = mesh.node_coord();
    Teuchos::RCP<VectorType> coords =
      Teuchos::rcp(new VectorType(x->getMap(), node_coords.dimension_1()));
    fill_coords(node_coords, coords->getDualView().d_view);

    std::string xmlFileName="muelu.xml";
    if (use_mean_based) {
      preconditioner =
        Teuchos::rcp(new MeanBasedPreconditioner<SM, LO, GO, N>());
      precOp = preconditioner->setupPreconditioner(A, xmlFileName, coords);
    }
    else {
      preconditioner =
        Teuchos::rcp(new MueLuPreconditioner<SM, LO, GO, N>());
      precOp = preconditioner->setupPreconditioner(A, xmlFileName, coords);
    }
  }

  //--------------------------------
  // Set up linear solver
  RCP<ParameterList> belosParams = Teuchos::parameterList();
  // Read in any params from xml file
  Teuchos::updateParametersFromXmlFileAndBroadcast(
    "belos.xml", belosParams.ptr(),*A->getComm());

  if (!(belosParams->isParameter("Convergence Tolerance")))
    belosParams->set("Convergence Tolerance", tolerance);
  if (!(belosParams->isParameter("Maximum Iterations")))
    belosParams->set("Maximum Iterations", Teuchos::as<int>(max_iter));
  if (!(belosParams->isParameter("Output Frequency")))
    belosParams->set("Output Frequency", 1);

  RCP<ProblemType> problem = rcp(new ProblemType(A, x, b));
  RCP<SolverType> solver = rcp(new SolverType(problem, belosParams));

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
  const Teuchos::RCP<Tpetra::CrsMatrix<SM,LO,GO,N> >& A,
  const Teuchos::RCP<Tpetra::Vector<SV,LO,GO,N> >& b,
  const Teuchos::RCP<Tpetra::Vector<SV,LO,GO,N> >& x,
  const Mesh& mesh,
  const int use_muelu,
  const int use_mean_based,
  const unsigned max_iter = 200,
  const typename Kokkos::Details::ArithTraits<SV>::mag_type tolerance =
    Kokkos::Details::ArithTraits<SV>::epsilon())
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
