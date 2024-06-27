// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//-----------------------------------------------------
// Ifpack2::SupportGraph is an implementation
// of Vaidya's maximum weight spanning tree preconditioner
//------------------------------------------------------

#ifndef IFPACK2_SUPPORTGRAPH_DEF_HPP
#define IFPACK2_SUPPORTGRAPH_DEF_HPP

// Ifpack2's CMake system should (and does) prevent Trilinos from
// attempting to build or install this class, if Lemon is not enabled.
// We check for this case regardless, in order to catch any bugs that
// future development might introduce in the CMake scripts.

#ifdef HAVE_IFPACK2_LEMON
#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#else
#  error "Ifpack2::SupportGraph requires that Trilinos be built with Lemon support."
#endif // HAVE_IFPACK2_LEMON

#include "Ifpack2_Heap.hpp"
#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>

namespace Ifpack2 {

template <class MatrixType>
SupportGraph<MatrixType>::
SupportGraph (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero()),
  Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one()),
  Randomize_ (1),
  NumForests_ (1),
  KeepDiag_ (Teuchos::ScalarTraits<magnitude_type>::one()),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{}


template <class MatrixType>
SupportGraph<MatrixType>::~SupportGraph () {}


template <class MatrixType>
void SupportGraph<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  magnitude_type absThresh = STM::zero();
  magnitude_type relThresh = STM::one();

  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double>
                                    ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double>
                                    ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try{
    Randomize_ = params.get<int> ("MST: randomize");
  }
  catch (InvalidParameterName&) {
  }

  if (absThresh != -STM::one()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one()) {
    Rthresh_ = relThresh;
  }

  try{
    NumForests_ = params.get<int> ("MST: forest number");
  }
  catch (InvalidParameterName&) {
  }

  try{
    KeepDiag_ = params.get<magnitude_type> ("MST: keep diagonal");
  }
  catch (InvalidParameterName&) {
  }
}



template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
SupportGraph<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null(), std::runtime_error, "Ifpack2::SupportGraph::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
SupportGraph<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
SupportGraph<MatrixType>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null(), std::runtime_error, "Ifpack2::ILUT::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
SupportGraph<MatrixType>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap();
}


template <class MatrixType>
bool SupportGraph<MatrixType>::hasTransposeApply() const {
  return true;
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumInitialize() const
{
  return(NumInitialize_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumCompute () const {
  return(NumCompute_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumApply () const {
  return(NumApply_);
}


template <class MatrixType>
double SupportGraph<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template<class MatrixType>
double SupportGraph<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template<class MatrixType>
double SupportGraph<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}


template<class MatrixType>
void SupportGraph<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null() && A->getComm()->getSize() == 1 &&
    A->getLocalNumRows() != A->getLocalNumCols(),
    std::runtime_error, "Ifpack2::ILUT::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getLocalNumRows() << " rows and "
    << A->getLocalNumCols() << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous support
  // graph computation.
  IsInitialized_ = false;
  IsComputed_ = false;
  A_local_ = Teuchos::null;
  Support_ = Teuchos::null;
  solver_ = Teuchos::null;
  A_ = A;
}



template<class MatrixType>
void
SupportGraph<MatrixType>::findSupport ()
{
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vec_type;



  const scalar_type zero = STS::zero();
  const scalar_type one = STS::one();

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));

  size_t num_verts = A_local_->getLocalNumRows();
  size_t num_edges
    = (A_local_->getLocalNumEntries() - A_local_->getLocalNumDiags())/2;


  // Create data structures for the BGL code
  // and temp data structures for extraction
  lemon::ListGraph graph;
  for (size_t row = 0; row < num_verts; ++row) {
    graph.addNode();
  }
  lemon::ListGraph::EdgeMap<magnitude_type> edgeWeights(graph);

  size_t num_entries;
  size_t max_num_entries = A_local_->getLocalMaxNumRowEntries();

  std::vector<scalar_type> valuestemp (max_num_entries);
  std::vector<local_ordinal_type> indicestemp (max_num_entries);

  std::vector<magnitude_type> diagonal (num_verts);

  Tpetra::ArrayView<scalar_type> values (valuestemp);
  Tpetra::ArrayView<local_ordinal_type> indices (indicestemp);

  // Extract from the tpetra matrix keeping only one edge per pair
  // (assume symmetric)
  size_t offDiagCount = 0;
  for (size_t row = 0; row < num_verts; ++row) {
    A_local_->getLocalRowCopy (row, indices, values, num_entries);
    for (size_t colIndex = 0; colIndex < num_entries; ++colIndex) {
      if(row == Teuchos::as<size_t>(indices[colIndex])) {
        diagonal[row] = values[colIndex];
      }

      if((row < Teuchos::as<size_t>(indices[colIndex]))
         && (values[colIndex] < zero)) {
        lemon::ListGraph::Edge edge
          = graph.addEdge(graph.nodeFromId(row),
                          graph.nodeFromId(Teuchos::as<size_t>
                                           (indices[colIndex])));
        edgeWeights[edge] = values[colIndex];

        if (Randomize_) {
          // Add small random pertubation.
          edgeWeights[edge] *= one +
            STS::magnitude(STS::rmin() * STS::random());
        }

        offDiagCount++;
      }
    }
  }


  // Run Kruskal, actually maximal weight ST since edges are negative
  std::vector<lemon::ListGraph::Edge> spanningTree;
  lemon::kruskal(graph, edgeWeights, std::back_inserter(spanningTree));

  // Create array to store the exact number of non-zeros per row
  Teuchos::ArrayRCP<size_t> NumNz (num_verts, 1);

  // Find the degree of all the vertices
  for (size_t i = 0; i != spanningTree.size(); ++i) {
    lemon::ListGraph::Edge e = spanningTree[i];

    local_ordinal_type localsource = graph.id(graph.u(e));
    local_ordinal_type localtarget = graph.id(graph.v(e));

    // We only want upper triangular entries, might need to swap
    if (localsource > localtarget) {
      localsource = localtarget;
      localtarget = graph.id(graph.u(e));
    }

    NumNz[localsource] += 1;
  }


  // Create an stl vector of stl vectors to hold indices and values
  std::vector<std::vector<local_ordinal_type> > Indices (num_verts);
  std::vector<std::vector<magnitude_type> > Values (num_verts);

  for (size_t i = 0; i < num_verts; ++i) {
    Indices[i].resize(NumNz[i]);
    Values[i].resize(NumNz[i]);
  }

  // The local ordering might be different from the global
  // ordering and we need the local number of non-zeros per row
  // to correctly allocate the preconditioner matrix memory
  Teuchos::ArrayRCP<size_t> localnumnz (num_verts, 1);

  for (size_t i = 0; i < num_verts; ++i) {
     Indices[i][0] = i;
  }


  // Add each spanning forest (tree) to the support graph and
  // remove it from original graph
  for (int i = 0; i < NumForests_; ++i) {
    // If a tree has already been added then we need to rerun Kruskall and
    // update the arrays containing size information
    if (i > 0) {
      spanningTree.clear();
      lemon::kruskal(graph, edgeWeights, std::back_inserter(spanningTree));

      for (size_t i = 0; i != spanningTree.size(); ++i) {
        NumNz[graph.id(graph.u(spanningTree[i]))] += 1;
      }

      // FIXME (mfh 14 Nov 2013) Are you sure that all this resizing
      // is a good idea?
      for (size_t i = 0; i < num_verts; ++i) {
        Indices[i].resize(NumNz[i]);
        Values[i].resize(NumNz[i]);
      }
    }

    for (size_t i = 0; i != spanningTree.size(); ++i) {
      lemon::ListGraph::Edge e = spanningTree[i];

      local_ordinal_type localsource = graph.id(graph.u(e));
      local_ordinal_type localtarget = graph.id(graph.v(e));

      if (localsource > localtarget) {
        localsource = localtarget;
        localtarget = graph.id(graph.u(e));
      }

      // Assume standard Laplacian with constant row-sum.
      // Edge weights are negative, so subtract to make diagonal positive
      Values[localtarget][0] -= edgeWeights[e];
      Values[localsource][0] -= edgeWeights[e];

      Indices[localsource][localnumnz[localsource]] = localtarget;
      Values[localsource][localnumnz[localsource]] = edgeWeights[e];
      localnumnz[localsource] += 1;

      graph.erase(e);
    }
  }

  // Set diagonal to weighted average of Laplacian preconditioner
  // and the original matrix

  // First compute the "diagonal surplus" (in the original input matrix)
  // If input is a (pure, Dirichlet) graph Laplacian , this will be 0
  vec_type ones (A_local_->getDomainMap());
  vec_type surplus (A_local_->getRangeMap());

  ones.putScalar(one);
  A_local_->apply(ones, surplus);

  Teuchos::ArrayRCP<const scalar_type> surplusaccess = surplus.getData(0);

  for (size_t i = 0; i < num_verts; ++i) {
    if (surplusaccess[i] > zero) {
      Values[i][0] += surplusaccess[i];
    }

    // If the original diagonal is less than the row sum then we aren't going to
    // use it regardless of the diagonal option, shouldn't happen for proper
    // Laplacian
    if (diagonal[i] < Values[i][0]) {
      diagonal[i] = Values[i][0];
    }

    Values[i][0] = KeepDiag_*diagonal[i] + (one-KeepDiag_) * Values[i][0];

    // Modify the diagonal with user specified scaling
    if (Rthresh_) {
      Values[i][0] *= Rthresh_;
    }
    if (Athresh_) {
      Values[i][0] += Athresh_;
    }
  }

  // Create the CrsMatrix for the support graph
  Support_ = rcp (new crs_matrix_type (A_local_->getRowMap(),
                                       A_local_->getColMap(),
                                       localnumnz));

  // Fill in the matrix with the stl vectors for each row
  for (size_t row = 0; row < num_verts; ++row) {
    Teuchos::ArrayView<local_ordinal_type>
      IndicesInsert (Indices[Teuchos::as<local_ordinal_type> (row)]);
    Teuchos::ArrayView<scalar_type>
      ValuesInsert (Values[Teuchos::as<local_ordinal_type> (row)]);
    Support_->insertLocalValues (row, IndicesInsert, ValuesInsert);
  }

  Support_->fillComplete();

}

template<class MatrixType>
Teuchos::RCP<const typename SupportGraph<MatrixType>::row_matrix_type>
SupportGraph<MatrixType>::
makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A->getComm()->getSize() > 1) {
    return Teuchos::rcp (new LocalFilter<row_matrix_type> (A));
  } else {
    return A;
  }
}

template<class MatrixType>
void SupportGraph<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null()) {
    timer = TimeMonitor::getNewCounter(timerName);
  }
  double startTime = timer->wallTime();
  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null(), std::runtime_error, "Ifpack2::SupportGraph::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;
    Support_ = Teuchos::null;
    solver_ = Teuchos::null;

    A_local_ = makeLocalFilter(A_); // Compute the local filter.
    findSupport(); // Compute the support.

    // Set up the solver and compute the symbolic factorization.
    solver_ = Amesos2::create<crs_matrix_type, MV> ("amesos2_cholmod", Support_);
    solver_->symbolicFactorization();

    IsInitialized_ = true;
    ++NumInitialize_;
  } // Stop timing here.

  InitializeTime_ += (timer->wallTime() - startTime);
}



template<class MatrixType>
void SupportGraph<MatrixType>::compute () {
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Don't count initialization in the compute() time.
  if (! isInitialized()) {
    initialize();
  }

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter(timerName);
  if (timer.is_null()) {
    timer = TimeMonitor::getNewCounter(timerName);
  }
  double startTime = timer->wallTime();
  { // Start timing here.
    Teuchos::TimeMonitor timeMon (*timer);
    solver_->numericFactorization();
    IsComputed_ = true;
    ++NumCompute_;
  } // Stop timing here.

  ComputeTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
void
SupportGraph<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,
                                 local_ordinal_type,
                                 global_ordinal_type,
                                 node_type>& X,
       Tpetra::MultiVector<scalar_type,
                           local_ordinal_type,
                           global_ordinal_type,
                           node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  typedef scalar_type DomainScalar;
  typedef scalar_type RangeScalar;
  typedef Tpetra::MultiVector<DomainScalar, local_ordinal_type,
    global_ordinal_type, node_type> MV;

  RCP<FancyOStream> out = getFancyOStream(rcpFromRef(std::cout));

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter(timerName);
  if (timer.is_null()) {
    timer = TimeMonitor::getNewCounter(timerName);
  }
  double startTime = timer->wallTime();
  { // Start timing here.
    Teuchos::TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed(), std::runtime_error,
      "Ifpack2::SupportGraph::apply: You must call compute() to compute the "
      "incomplete factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::SupportGraph::apply: X and Y must have the same number of "
      "columns.  X has " << X.getNumVectors() << " columns, but Y has "
      << Y.getNumVectors() << " columns.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      beta != STS::zero(), std::logic_error,
      "Ifpack2::SupportGraph::apply: This method does not currently work when "
      "beta != 0.");

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    {
      if (X.aliases(Y)) {
        Xcopy = rcp (new MV (X, Teuchos::Copy));
      } else {
        Xcopy = rcpFromRef (X);
      }
    }

    if (alpha != STS::one ()) {
      Y.scale (alpha);
    }

    RCP<MV> Ycopy = rcpFromRef (Y);

    solver_->setB (Xcopy);
    solver_->setX (Ycopy);
    solver_->solve ();
  } // Stop timing here.

  ++NumApply_;

  ApplyTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
std::string SupportGraph<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::SupportGraph\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Matrix: not null"
       << ", Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]";
  }

  os << "}";
  return os.str ();
}


template <class MatrixType>
void SupportGraph<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl
    = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
  Teuchos::OSTab tab (out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if(vl != VERB_NONE && getComm()->getRank() == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===================================================================="
      "===========" << endl;
    out << "Absolute threshold: " << getAbsoluteThreshold() << endl;
    out << "Relative threshold: " << getRelativeThreshold() << endl;

    if (isComputed()) {
      out << "Number of nonzeros in A: " << A_->getGlobalNumEntries() << endl;
      out << "Number of nonzeros in A_local: "
          << A_local_->getGlobalNumEntries() << endl;
      out << "Number of edges in support graph: "
          << Support_->getGlobalNumEntries() - Support_->getGlobalNumDiags()
          << endl;

      const double popFrac =
        static_cast<double> (Support_->getGlobalNumEntries() -
                             Support_->getGlobalNumDiags()) /
        (A_local_->getGlobalNumEntries() - A_local_->getGlobalNumDiags());

      out << "Fraction of off diagonals of supportgraph/off diagonals of "
        "original: " << popFrac << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s) " << endl;
    out << "------------    -------    ---------------" << endl;
    out << "initialize()    " << setw(7) << getNumInitialize() << "    "
        << setw(15) << getInitializeTime() << endl;
    out << "compute()       " << setw(7) << getNumCompute()    << "    "
        << setw(15) << getComputeTime()    << endl;
    out << "apply()         " << setw(7) << getNumApply()      << "    "
        << setw(15) << getApplyTime()      << endl;
    out << "===================================================================="
      "===========" << endl;
    out << endl;

    solver_->printTiming(out, verbLevel);
  }
}


}//namespace Ifpack2

#endif /* IFPACK2_SUPPORTGRAPH_DEF_HPP */

