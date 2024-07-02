// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_ComputeGatherMap.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_computeRowAndColumnOneNorms.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_leftAndOrRightScaleCrsMatrix.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_LAPACK.hpp"
#include "KokkosBlas1_abs.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas3_gemm.hpp"

#if defined(HAVE_IFPACK2_AZTECOO)
#  include "AztecOO.h"
#endif // defined(HAVE_IFPACK2_AZTECOO)

#if defined(HAVE_IFPACK2_EPETRA)
#  include "Epetra_Comm.h"
#  if defined(HAVE_MPI)
#    include "mpi.h"
#    include "Epetra_MpiComm.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif
#  include "Epetra_CrsMatrix.h"
#  include "Epetra_Map.h"
#  include "Epetra_Vector.h"
#  if defined(HAVE_TPETRACORE_MPI)
#    include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  endif // HAVE_TPETRACORE_MPI
#endif // defined(HAVE_IFPACK2_EPETRA)

#include <algorithm> // std::transform
#include <cctype> // std::toupper
#include <functional>
#include <memory> // std::unique_ptr
#include <sstream>
#include <tuple>

namespace { // (anonymous)

// See example here:
//
// http://en.cppreference.com/w/cpp/string/byte/toupper
std::string stringToUpper (std::string s)
{
  std::transform (s.begin (), s.end (), s.begin (),
                  [] (unsigned char c) { return std::toupper (c); });
  return s;
}

#if defined(HAVE_IFPACK2_EPETRA)

// I know it looks weird to use a unique_ptr here, but Epetra_Comm
// actually has its own shallow-copy semantics, so we don't need
// shared_ptr here.
std::unique_ptr<Epetra_Comm>
tpetraToEpetraComm (const Teuchos::Comm<int>& tpetraComm)
{
#if defined(HAVE_TPETRACORE_MPI) && ! defined(HAVE_MPI)
#  error "MPI is enabled in Tpetra, but is not enabled in Epetra."
#elif ! defined(HAVE_TPETRACORE_MPI) && defined(HAVE_MPI)
#  error "MPI is not enabled in Tpetra, but is enabled in Epetra."
#endif

#if defined(HAVE_TPETRACORE_MPI) && defined(HAVE_MPI)
  MPI_Comm rawMpiComm = Tpetra::Details::extractMpiCommFromTeuchos (tpetraComm);
  std::unique_ptr<Epetra_MpiComm> retComm (new Epetra_MpiComm (rawMpiComm));
#else
  std::unique_ptr<Epetra_SerialComm> retComm (new Epetra_SerialComm);
#endif

  return std::move (retComm);
}

// Given a Tpetra::Map map_t, and comm_e, the result of
// tpetraToEpetraComm(*(map_t.getComm())), return the corresponding
// Epetra_Map.
template<class LO, class GO, class NT>
Epetra_Map
tpetraToEpetraMap (const Tpetra::Map<LO, GO, NT>& map_t,
                   const Epetra_Comm& comm_e)
{
  const int gblNumInds = static_cast<int> (map_t.getGlobalNumElements ());
  const int lclNumInds = static_cast<int> (map_t.getLocalNumElements ());
  const int indexBase = static_cast<int> (map_t.getIndexBase ());

  if (map_t.isContiguous ()) {
    if (map_t.isUniform ()) {
      return Epetra_Map (gblNumInds, indexBase, comm_e);
    }
    else {
      return Epetra_Map (gblNumInds, lclNumInds, indexBase, comm_e);
    }
  }
  else {
    auto myGblInds = map_t.getMyGlobalIndices (); // returns host data
    return Epetra_Map (gblNumInds, lclNumInds, myGblInds.data (), indexBase, comm_e);
  }
}

template<class TpetraDistObjectType>
int
getRankSafelyFromTpetraDistObject (const TpetraDistObjectType& X)
{
  const auto map = X.getMap ();
  if (map.is_null ()) {
    return 0;
  }
  else {
    const auto comm = map->getComm ();
    if (comm.is_null ()) {
      return 0;
    }
    else {
      return comm->getRank ();
    }
  }
}

template<class LO, class GO, class NT>
std::pair<int, std::string>
deep_copy (Epetra_Vector& X_e,
           const Tpetra::Vector<double, LO, GO, NT>& X_t)
{
  using host_view_type = Kokkos::View<double*, Kokkos::LayoutLeft,
    Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  int lclErrCode = 0;

  const int lclNumRows = X_e.Map ().NumMyElements ();
  if (lclNumRows != static_cast<int> (X_t.getLocalLength ())) {
    lclErrCode = -1;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): X_t.getLocalLength()"
      " = " << X_t.getLocalLength () << " != map_e.NumMyElements() = "
            << lclNumRows << "." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }

  double* X_e_lcl_raw = nullptr;
  const int curErrCode = X_e.ExtractView (&X_e_lcl_raw);
  if (curErrCode != 0) {
    lclErrCode = -2;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): "
            << "Epetra_Vector::ExtractView failed with nonzero error code "
            << curErrCode << std::endl;
    return {lclErrCode, errStrm.str ()};
  }
  if (X_e_lcl_raw == nullptr && lclNumRows != 0) {
    lclErrCode = -3;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): "
      "Epetra_Vector::ExtractView returns success, but output a null pointer, "
      "even though the object's Map reports nonzero (" << lclNumRows << ") "
      "rows on the calling process." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }

  host_view_type X_e_lcl (X_e_lcl_raw, lclNumRows);
  if (X_t.need_sync_device ()) {
    auto X_t_lcl_2d = X_t.getLocalViewHost (Tpetra::Access::ReadOnly);
    auto X_t_lcl = Kokkos::subview (X_t_lcl_2d, Kokkos::ALL (), 0);
    Kokkos::deep_copy (X_e_lcl, X_t_lcl);
  }
  else {
    auto X_t_lcl_2d = X_t.getLocalViewDevice (Tpetra::Access::ReadOnly);
    auto X_t_lcl = Kokkos::subview (X_t_lcl_2d, Kokkos::ALL (), 0);
    Kokkos::deep_copy (X_e_lcl, X_t_lcl);
  }

  return {lclErrCode, ""};
}

template<class LO, class GO, class NT>
std::pair<int, std::string>
deep_copy (Tpetra::Vector<double, LO, GO, NT>& X_t,
           const Epetra_Vector& X_e)
{
  using host_view_type = Kokkos::View<const double*, Kokkos::LayoutLeft,
    Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  int lclErrCode = 0;

  const int lclNumRows = X_e.Map ().NumMyElements ();
  if (lclNumRows != static_cast<int> (X_t.getLocalLength ())) {
    lclErrCode = -1;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): X_t.getLocalLength()"
      " = " << X_t.getLocalLength () << " != map_e.NumMyElements() = "
            << lclNumRows << "." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }

  double* X_e_lcl_raw = nullptr;
  const int curErrCode = X_e.ExtractView (&X_e_lcl_raw);
  if (curErrCode != 0) {
    lclErrCode = -2;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): "
            << "Epetra_Vector::ExtractView failed with nonzero error code "
            << curErrCode << std::endl;
    return {lclErrCode, errStrm.str ()};
  }
  if (X_e_lcl_raw == nullptr && lclNumRows != 0) {
    lclErrCode = -3;
    const int myRank = getRankSafelyFromTpetraDistObject (X_t);
    std::ostringstream errStrm;
    errStrm << "Proc " << myRank << ": deep_copy(Vector): "
      "Epetra_Vector::ExtractView returns success, but output a null pointer, "
      "even though the object's Map reports nonzero (" << lclNumRows << ") "
      "rows on the calling process." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }

  host_view_type X_e_lcl (X_e_lcl_raw, lclNumRows);
  if (X_t.need_sync_device ()) {
    auto X_t_lcl_2d = X_t.getLocalViewHost (Tpetra::Access::OverwriteAll);
    auto X_t_lcl = Kokkos::subview (X_t_lcl_2d, Kokkos::ALL (), 0);
    Kokkos::deep_copy (X_t_lcl, X_e_lcl);
  }
  else {
    auto X_t_lcl_2d = X_t.getLocalViewDevice (Tpetra::Access::OverwriteAll);
    auto X_t_lcl = Kokkos::subview (X_t_lcl_2d, Kokkos::ALL (), 0);
    Kokkos::deep_copy (X_t_lcl, X_e_lcl);
  }

  return {lclErrCode, ""};
}

template<class LO, class GO, class NT>
std::pair<int, std::string>
deep_copy (Epetra_MultiVector& X_e,
           const Tpetra::MultiVector<double, LO, GO, NT>& X_t)
{
  const int myRank = getRankSafelyFromTpetraDistObject (X_t);
  std::ostringstream errStrm;
  int lclErrCode = 0;

  // This is supposed to be globally consistent, so it should throw.
  const int numCols = X_e.NumVectors ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (X_t.getNumVectors () != std::size_t (numCols), std::invalid_argument,
     "deep_copy(MultiVector): X_t.getNumVectors() = " << X_t.getNumVectors ()
     << " != X_e.NumVectors() = " << numCols << ".");

  const int lclNumRows = X_e.Map ().NumMyElements ();
  if (lclNumRows != static_cast<int> (X_t.getLocalLength ())) {
    lclErrCode = 0;
    errStrm << "Proc " << myRank << ": deep_copy(MultiVector): "
      "X_t.getLocalLength() = " << X_t.getLocalLength ()
      << " != map_e.NumMyElements() = " << lclNumRows << "." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }
  else {
    for (int col = 0; col < numCols; ++col) {
      Epetra_Vector* X_e_cur = X_e(col);
      auto X_t_cur = X_t.getVector (col);
      const std::pair<int, std::string> result = deep_copy (*X_e_cur, *X_t_cur);
      if (result.first != 0) {
        errStrm << "Proc " << myRank << ": deep_copy(MultiVector): For column "
          << col << ", deep_copy failed with local error code " << lclErrCode
          << " and the following message: " << result.second;
        return {lclErrCode, errStrm.str ()};
      }
    }
    return {lclErrCode, ""};
  }
}

template<class LO, class GO, class NT>
std::pair<int, std::string>
deep_copy (Tpetra::MultiVector<double, LO, GO, NT>& X_t,
           const Epetra_MultiVector& X_e)
{
  const int myRank = getRankSafelyFromTpetraDistObject (X_t);
  std::ostringstream errStrm;
  int lclErrCode = 0;

  // This is supposed to be globally consistent, so it should throw.
  const int numCols = X_e.NumVectors ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (X_t.getNumVectors () != std::size_t (numCols), std::invalid_argument,
     "deep_copy(MultiVector): X_t.getNumVectors() = " << X_t.getNumVectors ()
     << " != X_e.NumVectors() = " << numCols << ".");

  const int lclNumRows = X_e.Map ().NumMyElements ();
  if (lclNumRows != static_cast<int> (X_t.getLocalLength ())) {
    lclErrCode = 0;
    errStrm << "Proc " << myRank << ": deep_copy(MultiVector): "
      "X_t.getLocalLength() = " << X_t.getLocalLength ()
      << " != map_e.NumMyElements() = " << lclNumRows << "." << std::endl;
    return {lclErrCode, errStrm.str ()};
  }
  else {
    for (int col = 0; col < numCols; ++col) {
      const Epetra_Vector* X_e_cur = X_e(col);
      auto X_t_cur = X_t.getVectorNonConst (col);
      const std::pair<int, std::string> result = deep_copy (*X_t_cur, *X_e_cur);
      if (result.first != 0) {
        errStrm << "Proc " << myRank << ": deep_copy(MultiVector): For column "
          << col << ", deep_copy failed with local error code " << lclErrCode
          << " and the following message: " << result.second;
        return {lclErrCode, errStrm.str ()};
      }
    }
    return {lclErrCode, ""};
  }
}

template<class LO, class GO, class NT>
Epetra_Vector
tpetraToEpetraVector (const Tpetra::MultiVector<double, LO, GO, NT>& X_t,
                      const Epetra_Map& map_e)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (X_t.getNumVectors () != std::size_t (1), std::invalid_argument,
     "This function requires that X_t have only one column.");

  const int lclNumRows = map_e.NumMyElements ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (lclNumRows != static_cast<int> (X_t.getLocalLength ()),
     std::runtime_error, "X_t.getLocalLength() = " << X_t.getLocalLength ()
     << " != map_e.NumMyElements() = " << lclNumRows << ".");

  Epetra_Vector X_e (map_e);
  const std::pair<int, std::string> result = deep_copy (X_e, X_t);
  using Teuchos::outArg;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MAX;

  // This function is a kind of constructor for an Epetra_Vector, so
  // we shouldn't try to return an error code.  This is why we
  // synchronize on the error code here.
  int lclErrCode = result.first > 0 ? -result.first : result.first;
  int gblErrCode = 0;
  if (! X_t.getMap ().is_null () && ! X_t.getMap ()->getComm ().is_null ()) {
    const auto& comm = * (X_t.getMap ()->getComm ());
    reduceAll<int, int> (comm, REDUCE_MAX, lclErrCode, outArg (gblErrCode));
    if (gblErrCode != 0) {
      std::ostringstream err;
      Tpetra::Details::gathervPrint (err, result.second, comm);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err.str ());
    }
  }

  return X_e;
}

template<class LO, class GO, class NT>
Epetra_CrsMatrix
tpetraToEpetraCrsMatrix (const Tpetra::CrsMatrix<double, LO, GO, NT>& A_t,
                         const Epetra_Map& rowMap,
                         const Epetra_Map& colMap,
                         const Epetra_Map& domMap,
                         const Epetra_Map& ranMap)
{
  static_assert (std::is_same<LO, int>::value, "This code only works when LO = int.");
  const LO lclNumRows = rowMap.NumMyElements ();

  std::vector<int> numEntPerRow (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const size_t numEnt = A_t.getNumEntriesInLocalRow (lclRow);
    numEntPerRow[lclRow] = static_cast<int> (numEnt);
  }
  // We can use static profile, since we know the structure in advance.
  Epetra_CrsMatrix A_e (Copy, rowMap, colMap, numEntPerRow.data (), true);

  using tmatrix_t = Tpetra::CrsMatrix<double, LO, GO, NT>;

  typename tmatrix_t::nonconst_local_inds_host_view_type 
           lclColInds ("ifpack2::lclColInds", A_t.getLocalMaxNumRowEntries());
  typename tmatrix_t::nonconst_values_host_view_type 
           vals ("ifpack2::vals", A_t.getLocalMaxNumRowEntries());

  int lclErrCode = 0;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    size_t numEnt = A_t.getNumEntriesInLocalRow (lclRow);

    A_t.getLocalRowCopy (static_cast<LO> (lclRow), lclColInds, vals, numEnt);
    lclErrCode = A_e.InsertMyValues (lclRow, static_cast<int> (numEnt),
      vals.data(), lclColInds.data());
    if (lclErrCode != 0) {
      break;
    }
  }

  if (lclErrCode > 0) {
    lclErrCode = -lclErrCode; // make sure that REDUCE_MIN works
  }
  using Teuchos::outArg;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  int gblErrCode = 0; // output argument
  reduceAll<int, int> (* (A_t.getMap ()->getComm ()), REDUCE_MIN,
                       lclErrCode, outArg (gblErrCode));
  TEUCHOS_TEST_FOR_EXCEPTION
    (gblErrCode != 0, std::runtime_error, "tpetraToEpetraCrsMatrix: "
     "InsertMyValues failed on some process!");

  try {
    lclErrCode = A_e.FillComplete (domMap, ranMap);
  }
  catch (...) {
    // Epetra likes to throw integers, not std::exception subclasses.
    lclErrCode = -1;
  }
  if (lclErrCode > 0) {
    lclErrCode = -lclErrCode; // make sure that REDUCE_MIN works
  }
  reduceAll<int, int> (* (A_t.getMap ()->getComm ()), REDUCE_MIN,
                       lclErrCode, outArg (gblErrCode));
  TEUCHOS_TEST_FOR_EXCEPTION
    (gblErrCode != 0, std::runtime_error, "tpetraToEpetraCrsMatrix: "
     "Epetra_CrsMatrix::FillComplete failed on some process!");

  return A_e;
}

// Return Epetra versions of (A, x, b).
template<class LO, class GO, class NT>
std::tuple<Epetra_CrsMatrix, Epetra_Vector, Epetra_Vector>
tpetraToEpetraLinearSystem (const Tpetra::CrsMatrix<double, LO, GO, NT>& A_t,
                            const Tpetra::MultiVector<double, LO, GO, NT>& x_t,
                            const Tpetra::MultiVector<double, LO, GO, NT>& b_t)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A_t.isFillComplete (), std::invalid_argument,
     "The input matrix A must be fill complete.");

  Teuchos::RCP<const Teuchos::Comm<int>> comm_t =
    A_t.getMap ().is_null () ? Teuchos::null : A_t.getMap ()->getComm ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (comm_t.is_null (), std::runtime_error, "Tpetra's Comm is null.");

  const bool domMapsSame = A_t.getDomainMap ()->isSameAs (* (x_t.getMap ()));
  const bool ranMapsSame = A_t.getRangeMap ()->isSameAs (* (b_t.getMap ()));
  TEUCHOS_TEST_FOR_EXCEPTION
    (! domMapsSame, std::invalid_argument, "Domain Map of input matrix A "
     "is not the same as the Map of input vector x.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! ranMapsSame, std::invalid_argument, "Range Map of input matrix A "
     "is not the same as the Map of input vector b.");

  std::unique_ptr<Epetra_Comm> comm_e = tpetraToEpetraComm (*comm_t);

  Epetra_Map rowMap_e = tpetraToEpetraMap (* (A_t.getRowMap ()), *comm_e);
  Epetra_Map colMap_e = tpetraToEpetraMap (* (A_t.getColMap ()), *comm_e);
  Epetra_Map domMap_e = tpetraToEpetraMap (* (A_t.getDomainMap ()), *comm_e);
  Epetra_Map ranMap_e = tpetraToEpetraMap (* (A_t.getRangeMap ()), *comm_e);

  Epetra_CrsMatrix A_e = tpetraToEpetraCrsMatrix (A_t, rowMap_e, colMap_e,
                                                  domMap_e, ranMap_e);
  Epetra_Vector x_e = tpetraToEpetraVector (x_t, domMap_e);
  Epetra_Vector b_e = tpetraToEpetraVector (b_t, ranMap_e);

  // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
  //return {A_e, x_e, b_e};
  return std::make_tuple (A_e, x_e, b_e);
}

#endif // defined(HAVE_IFPACK2_EPETRA)

#if defined(HAVE_IFPACK2_AZTECOO) && defined(HAVE_IFPACK2_EPETRA)

std::tuple<double, int, bool>
solveEpetraLinearSystemWithAztecOO (std::tuple<Epetra_CrsMatrix, Epetra_Vector, Epetra_Vector>& sys,
                                    const std::string& solverType,
                                    const double convTol,
                                    const int maxNumIters)
{
  using std::get;
  Epetra_CrsMatrix& A = get<0> (sys);
  Epetra_Vector& x = get<1> (sys);
  Epetra_Vector& b = get<2> (sys);
  Epetra_LinearProblem problem (&A, &x, &b);
  AztecOO solver (problem);
  //(void) solver.SetAztecOption (AZ_precond, AZ_sym_GS);
  (void) solver.SetAztecOption (AZ_precond, AZ_none);

  const std::string upperSolverType = stringToUpper (solverType);
  int errCode = 0;
  int aztecOO_solverType = AZ_gmres;
  if (upperSolverType == "BICGSTAB") {
    aztecOO_solverType = AZ_bicgstab;
  }
  else if (upperSolverType == "CG") {
    aztecOO_solverType = AZ_cg;
  }
  else if (upperSolverType == "CGS") {
    aztecOO_solverType = AZ_cgs;
  }
  else if (upperSolverType == "GMRES") {
    aztecOO_solverType = AZ_gmres;
  }
  else if (upperSolverType == "TFQMR") {
    aztecOO_solverType = AZ_tfqmr;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "Linear solver type \"" << solverType
       << "\" not currently supported for AztecOO.");
  }
  errCode = solver.SetAztecOption (AZ_solver, aztecOO_solverType);
  if (errCode != 0) {
    // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
    //return {-1.0, 0, false};
    return std::make_tuple (-1.0, 0, false);
  }

  errCode = solver.Iterate (maxNumIters, convTol);
  bool success = true;
  if (errCode == 0) {
    success = true;
  }
  else if (errCode == 1) {
    success = false; // reached max number of iterations
  }
  else { // negative means error code
    success = false;
  }
  // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
  //return {solver.ScaledResidual (), solver.NumIters (), success};
  return std::make_tuple (solver.ScaledResidual (), solver.NumIters (), success);
}

#endif // defined(HAVE_IFPACK2_AZTECOO) && defined(HAVE_IFPACK2_EPETRA)

template<class MV>
typename MV::dot_type accurate_dot (const MV& X, const MV& Y)
{
  using LO = typename MV::local_ordinal_type;
  using dot_type = typename MV::dot_type;

  const LO lclNumRows = X.getLocalLength ();
  auto X_lcl_2d = X.getLocalViewHost(Tpetra::Access::ReadOnly);
  auto X_lcl = Kokkos::subview (X_lcl_2d, Kokkos::ALL (), 0);
  auto Y_lcl_2d = Y.getLocalViewHost(Tpetra::Access::ReadOnly);
  auto Y_lcl = Kokkos::subview (Y_lcl_2d, Kokkos::ALL (), 0);

  long double sum = 0.0;
  for (LO i = 0; i < lclNumRows; ++i) {
    const long double x_i = X_lcl (i);
    const long double y_i = Y_lcl (i);
    sum = std::fma (x_i, y_i, sum);
  }

  return dot_type (sum);
}

template<class MV>
auto norm (const MV& X) -> decltype (X.getVector (0)->norm2 ())
{
  return X.getVector (0)->norm2 ();
}

template<class MV, class OP>
void residual (MV& R, const MV& B, const OP& A, const MV& X)
{
  using STS = Teuchos::ScalarTraits<typename MV::scalar_type>;

  Tpetra::deep_copy (R, B);
  A.apply (X, R, Teuchos::NO_TRANS, -STS::one (), STS::one ());
}

template<class MV>
typename MV::dot_type dot (const MV& X, const MV& Y)
{
  // Teuchos::Array<typename MV::dot_type> dots (X.getNumVectors ());
  // Y.dot (X, dots ());
  // return dots[0];

  return accurate_dot (X, Y);
}

template<class MV>
bool
AZ_breakdown_f (MV& v, MV& w, const typename MV::dot_type v_dot_w)
{
  using STS = Teuchos::ScalarTraits<typename MV::scalar_type>;

  const auto v_norm = norm (v);
  const auto w_norm = norm (w);
  return STS::magnitude (v_dot_w) <= 100.0 * v_norm * w_norm * STS::eps ();
}

template<class MV, class OP>
std::tuple<typename MV::mag_type, int, bool>
bicgstab_aztecoo (MV& x,
                  const OP& A,
                  const OP* const M,
                  const MV& b,
                  const int max_it,
                  const typename MV::mag_type tol)
{
  using STS = Teuchos::ScalarTraits<typename MV::scalar_type>;
  using STM = Teuchos::ScalarTraits<typename MV::mag_type>;
  using dot_type = typename MV::dot_type;
  using mag_type = typename MV::mag_type;

  int iter = 0;

  bool brkdown_will_occur = false;
  dot_type alpha = STS::one ();
  dot_type beta = STS::zero ();
  dot_type omega = STS::one ();
  dot_type rhonm1 = STS::one ();
  dot_type rhon = STS::zero ();
  dot_type sigma = STS::zero ();
  mag_type brkdown_tol = STS::eps ();
  mag_type scaled_r_norm = -STM::one ();
  mag_type actual_residual = -STM::one ();
  mag_type rec_residual = -STM::one ();
  dot_type dtemp = STS::zero ();

  MV phat (b.getMap (), b.getNumVectors (), false);
  MV p (b.getMap (), b.getNumVectors (), false);
  MV shat (b.getMap (), b.getNumVectors (), false);
  MV s (b.getMap (), b.getNumVectors (), false);
  MV r (b.getMap (), b.getNumVectors (), false);
  MV r_tld (b.getMap (), b.getNumVectors (), false);
  MV v (b.getMap (), b.getNumVectors (), false);

  residual (r, b, A, x); // r = b - A*x;

  // "v, p <- 0"
  v.putScalar (STS::zero ());
  p.putScalar (STS::zero ());

  // "set rtilda" [sic]
  constexpr int my_AZ_aux_vec = 0; // AZ_resid = 0; AZ_rand (?) = 1
  constexpr int my_AZ_resid = 0;
  if (my_AZ_aux_vec == my_AZ_resid) {
    Tpetra::deep_copy (r_tld, r);
  }
  else {
    r_tld.randomize ();
  }

  // AZ_compute_global_scalars does all this, neatly bundled into a
  // single all-reduce.
  const mag_type b_norm = norm (b);
  actual_residual = norm (r);
  rec_residual = actual_residual;
  scaled_r_norm = rec_residual / b_norm;
  rhon = dot (r_tld, r);
  if (scaled_r_norm <= tol) {
    // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
    //return {scaled_r_norm, 0, true};
    return std::make_tuple (scaled_r_norm, 0, true);
  }

  for (iter = 1; iter <= max_it; ++iter) {
    //std::cerr << ">>> AztecOO-ish BiCGSTAB: iter = " << (iter - 1) << std::endl;
    if (brkdown_will_occur) {
      residual (v, b, A, x); // v = b - A*x
      actual_residual = norm (v);
      scaled_r_norm = actual_residual / b_norm;
      std::cerr << "Uh oh, breakdown" << std::endl;
      // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
      //return {scaled_r_norm, iter, false};
      return std::make_tuple (scaled_r_norm, iter, false);
    }

    beta = (rhon / rhonm1) * (alpha / omega);

    if (STS::magnitude (rhon) < brkdown_tol) {
      if (AZ_breakdown_f(r, r_tld, rhon)) {
        brkdown_will_occur = true;
      }
      else {
        brkdown_tol = 0.1 * STS::magnitude (rhon);
      }
    }

    rhonm1 = rhon;

    /* p    = r + beta*(p - omega*v)       */
    /* phat = M^-1 p                       */
    /* v    = A phat                       */

    dtemp = beta * omega;
    p.update (STS::one (), r, -dtemp, v, beta);
    Tpetra::deep_copy (phat, p);

    if (M != nullptr) {
      M->apply (p, phat);
    }
    A.apply (phat, v);
    sigma = dot (r_tld, v);

    if (STS::magnitude (sigma) < brkdown_tol) {
      if (AZ_breakdown_f(r_tld, v, sigma)) { // actual break down
        residual (v, b, A, x); // v = b - A*x;
        actual_residual = norm (v);
        scaled_r_norm = actual_residual / b_norm;
        // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
        //return {scaled_r_norm, iter, false};
        return std::make_tuple (scaled_r_norm, iter, false);
      }
      else {
        brkdown_tol = 0.1 * STS::magnitude (sigma);
      }
    }

    alpha = rhon / sigma;

    s.update (STS::one (), r, -alpha, v, STS::zero ());
    Tpetra::deep_copy (shat, s);

    if (M != nullptr) {
      M->apply (s, shat);
    }
    A.apply (shat, r);

    /* omega = (t,s)/(t,t) with r = t */

    const auto dot_vec_0 = dot (r, s);
    const auto dot_vec_1 = dot (r, r);
    if (STM::magnitude (dot_vec_1) < tol) {
      omega = STS::zero ();
      brkdown_will_occur = true;
    }
    else {
      omega = dot_vec_0 / dot_vec_1;
    }

    /* x = x + alpha*phat + omega*shat */
    /* r = s - omega*r */

    // DAXPY_F77(&N, &alpha, phat, &one, x, &one);
    // DAXPY_F77(&N, &omega, shat, &one, x, &one);
    x.update (alpha, phat, omega, shat, STS::one ());

    // for (i = 0; i < N; i++) r[i] = s[i] - omega * r[i];
    r.update (STS::one (), s, -omega);

    rec_residual = norm (r);
    scaled_r_norm = rec_residual / b_norm;
    rhon = dot (r, r_tld);
    if (scaled_r_norm <= tol) {
      // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
      //return {scaled_r_norm, 0, true};
      return std::make_tuple (scaled_r_norm, 0, true);
    }
  }

  // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
  //return {scaled_r_norm, iter, scaled_r_norm <= tol};
  return std::make_tuple (scaled_r_norm, iter, scaled_r_norm <= tol);
}

template<class MV, class OP>
std::tuple<typename MV::mag_type, int, bool>
bicgstab_no_prec_paper (MV& x,
                        const OP& A,
                        const MV& b,
                        const int max_it,
                        const typename MV::mag_type tol)
{
  using dot_type = typename MV::dot_type;
  using mag_type = typename MV::mag_type;
  using STS = Teuchos::ScalarTraits<dot_type>;
  using STM = Teuchos::ScalarTraits<mag_type>;

  x.putScalar (STS::zero ()); // just in case

  const mag_type b_norm = norm (b);
  if (b_norm == STS::zero ()) {
    x.putScalar (STS::zero ());
    // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
    //return {STM::zero (), 0, true};
    return std::make_tuple (STM::zero (), 0, true);
  }
  dot_type rho_cur = STS::one ();
  dot_type rho_prv = STS::one ();
  dot_type alpha = STS::one ();
  dot_type beta = STS::one ();
  dot_type omega = STS::one ();

  MV r (b.getMap (), b.getNumVectors (), false);
  residual (r, b, A, x);

  mag_type r_norm = norm (r);
  mag_type scaled_r_norm = r_norm / b_norm;
  if (r_norm / b_norm <= tol) {
    // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
    //return {scaled_r_norm, 0, true};
    return std::make_tuple (scaled_r_norm, 0, true);
  }
  MV r_hat (r, Teuchos::Copy);

  MV v (b.getMap (), b.getNumVectors ());
  MV p (b.getMap (), b.getNumVectors ());
  MV s (b.getMap (), b.getNumVectors ());
  MV t (b.getMap (), b.getNumVectors ());
  MV tmp (b.getMap (), b.getNumVectors ());

  for (int iter = 1; iter <= max_it; ++iter) {
    rho_cur = dot (r_hat, r);
    beta = (rho_cur / rho_prv) * (alpha / omega);

    // p = r + beta(p - omega*v); v is 0 on 1st iter
    if (iter > 1) {
      p.update (-omega, v, STS::one ()); // p = -omega*v + p
    }
    p.update (STS::one (), r, beta); // p = r + beta*p

    A.apply (p, v);

    alpha = rho_cur / dot (r_hat, v);

    // s = r - alpha*v
    Tpetra::deep_copy (s, r);
    s.update (-alpha, v, STS::one ()); // s = -alpha*v + s

    A.apply (s, t);

    omega = dot (t, s) / dot (t, t);

    // x = x + alpha*p + omega*s
    x.update (alpha, p, STS::one ()); // x = alpha*p + x
    x.update (omega, s, STS::one ()); // x = omega*s + s

    // r = s - omega*t
    Tpetra::deep_copy (r, s);
    r.update (-omega, t, STS::one ()); // r = -omega*t + r

    // check convergence
    r_norm = norm (r);
    scaled_r_norm = r_norm / b_norm;
    if (scaled_r_norm <= tol) {
      residual (tmp, b, A, x);
      const mag_type actual_r_norm = norm (tmp);
      const mag_type actual_scaled_r_norm = actual_r_norm / b_norm;
      std::cerr << "!!! Actual scaled r norm: " << actual_scaled_r_norm << std::endl;
      // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
      //return {scaled_r_norm, iter, true};
      return std::make_tuple (scaled_r_norm, iter, true);
    }

    rho_prv = rho_cur; // don't forget this!
  }

  {
    residual (tmp, b, A, x);
    const mag_type actual_r_norm = norm (tmp);
    const mag_type actual_scaled_r_norm = actual_r_norm / b_norm;
    std::cerr << "!!! Actual scaled r norm: " << actual_scaled_r_norm << std::endl;
  }
  // Clang accepts the commented-out line, but GCC 4.9.3 (-std=c++11) does not.
  //return {scaled_r_norm, max_it, scaled_r_norm <= tol};
  return std::make_tuple (scaled_r_norm, max_it, scaled_r_norm <= tol);
}

template<class SC, class LO, class GO, class NT>
std::pair<Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >,
          Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> > >
gatherCrsMatrixAndMultiVector (LO& errCode,
                               const Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                               const Tpetra::MultiVector<SC, LO, GO, NT>& B)
{
  using Tpetra::Details::computeGatherMap;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using export_type = Tpetra::Export<LO, GO, NT>;
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NT>;

  auto rowMap_gathered = computeGatherMap (A.getRowMap (), Teuchos::null);
  export_type exp (A.getRowMap (), rowMap_gathered);
  auto A_gathered =
    Teuchos::rcp (new crs_matrix_type (rowMap_gathered,
                                       A.getGlobalMaxNumRowEntries ()));
  A_gathered->doExport (A, exp, Tpetra::INSERT);
  auto domainMap_gathered = computeGatherMap (A.getDomainMap (), Teuchos::null);
  auto rangeMap_gathered = computeGatherMap (A.getRangeMap (), Teuchos::null);
  A_gathered->fillComplete (domainMap_gathered, rangeMap_gathered);

  auto B_gathered =
    Teuchos::rcp (new mv_type (rangeMap_gathered, B.getNumVectors ()));
  B_gathered->doExport (B, exp, Tpetra::ADD);

  return std::make_pair (A_gathered, B_gathered);
}

using host_device_type = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;

template<class SC = Tpetra::MultiVector<>::scalar_type,
         class LO = Tpetra::MultiVector<>::local_ordinal_type,
         class GO = Tpetra::MultiVector<>::global_ordinal_type,
         class NT = Tpetra::MultiVector<>::node_type>
using HostDenseMatrix =
  Kokkos::View<typename Tpetra::MultiVector<SC, LO, GO, NT>::impl_scalar_type**,
               Kokkos::LayoutLeft,
               host_device_type>;

template<class SC, class LO, class GO, class NT>
HostDenseMatrix<SC, LO, GO, NT>
densifyGatheredCrsMatrix (LO& errCode,
                          const Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                          const std::string& label)
{
  const LO numRows = LO (A.getRangeMap ()->getLocalNumElements ());
  const LO numCols = LO (A.getDomainMap ()->getLocalNumElements ());
  using lids_type = typename Tpetra::CrsMatrix<SC, LO, GO, NT>::local_inds_host_view_type;
  using vals_type = typename Tpetra::CrsMatrix<SC, LO, GO, NT>::values_host_view_type;

  using dense_matrix_type = HostDenseMatrix<SC, LO, GO, NT>;
  dense_matrix_type A_dense (label, numRows, numCols);

  for (LO lclRow = 0; lclRow < numRows; ++lclRow) {
    lids_type lclColInds;
    vals_type vals;
    A.getLocalRowView (lclRow, lclColInds, vals);
    LO numEnt = vals.size();
    for (LO k = 0; k < numEnt; ++k) {
      const LO lclCol = lclColInds[k];
      using impl_scalar_type =
        typename Tpetra::CrsMatrix<SC, LO, GO, NT>::impl_scalar_type;
      A_dense(lclRow, lclCol) += impl_scalar_type (vals[k]);
    }   
  }

  return A_dense;
}

template<class SC, class LO, class GO, class NT>
HostDenseMatrix<SC, LO, GO, NT>
copyGatheredMultiVector (Tpetra::MultiVector<SC, LO, GO, NT>& X,
                         const std::string& label)
{
  using dense_matrix_type = HostDenseMatrix<SC, LO, GO, NT>;

  auto X_lcl = X.getLocalViewHost (Tpetra::Access::ReadOnly);
  dense_matrix_type X_copy (label, X.getLocalLength (), X.getNumVectors ());
  Kokkos::deep_copy (X_copy, X_lcl);

  return X_copy;
}

template<class SC, class LO, class GO, class NT>
LO
gatherAndDensify (HostDenseMatrix<SC, LO, GO, NT>& A_dense,
                  HostDenseMatrix<SC, LO, GO, NT>& B_dense,
                  const Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                  Tpetra::MultiVector<SC, LO, GO, NT>& B)
{
  LO errCode = 0;
  auto A_and_B_gathered = gatherCrsMatrixAndMultiVector (errCode, A, B);

  if (errCode == 0) {
    A_dense = densifyGatheredCrsMatrix (errCode, * (A_and_B_gathered.first), "A_dense");
    B_dense = copyGatheredMultiVector (B, "B_dense");
  }
  return errCode;
}

template<class DenseMatrixType>
void
solveLeastSquaresProblemAndReport (DenseMatrixType A,
                                   DenseMatrixType B,
                                   const double RCOND = -1.0 /* negative means use machine precision */)
{
  const int numRows = int (A.extent (0));
  const int numCols = int (A.extent (1));
  const int NRHS = int (B.extent (1));

  DenseMatrixType A_copy ("A_copy", numRows, numCols);
  Kokkos::deep_copy (A_copy, A);

  DenseMatrixType B_copy ("B_copy", numRows, NRHS);
  Kokkos::deep_copy (B_copy, B);
  //DenseMatrixType X ("X", numCols, NRHS);

  const int LDA = (numRows == 0) ? 1 : int (A_copy.stride (1));
  const int LDB = (numRows == 0) ? 1 : int (B_copy.stride (1));

  std::vector<double> S (std::min (numRows, numCols));
  int RANK = 0;
  int LWORK = -1; // workspace query
  int INFO = 0;
  Teuchos::LAPACK<int, double> lapack;

  using std::cerr;
  using std::cout;
  using std::endl;
  cout << "Solver:" << endl
       << "  Solver type: LAPACK's DGELSS" << endl;

  std::vector<double> WORK (1);
  lapack.GELSS (numRows, numCols, NRHS, A_copy.data (), LDA,
                B_copy.data (), LDB,
                S.data (), RCOND, &RANK, WORK.data (), LWORK, &INFO);
  if (INFO != 0) {
    cerr << "DGELSS returned INFO = " << INFO << " != 0." << endl;
    return;
  }
  LWORK = int (WORK[0]);
  if (LWORK < 0) {
    cerr << "DGELSS reported LWORK = " << LWORK << " < 0." << endl;
    return;
  }
  WORK.resize (LWORK);
  lapack.GELSS (numRows, numCols, NRHS, A_copy.data (), LDA,
                B_copy.data (), LDB,
                S.data (), RCOND, &RANK, WORK.data (), LWORK, &INFO);

  cout << "Results:" << endl
       << "  INFO: " << INFO << endl
       << "  RCOND: " << RCOND << endl
       << "  Singular values: ";
  for (double sigma : S) {
    cout << sigma << " ";
  }
  cout << endl;

  if (numRows == numCols) {
    std::vector<double> B_norms (NRHS);
    for (int k = 0; k < NRHS; ++k) {
      B_norms[k] =
        KokkosBlas::nrm2 (Kokkos::subview (B, Kokkos::ALL (), k));
    }

    auto X = B_copy;
    DenseMatrixType R ("R", numRows, NRHS);
    Kokkos::deep_copy (R, B);
    KokkosBlas::gemm ("N", "N", -1.0, A, X, +1.0, R);

    std::vector<double> explicitResidualNorms (NRHS);
    for (int k = 0; k < NRHS; ++k) {
      explicitResidualNorms[k] =
        KokkosBlas::nrm2 (Kokkos::subview (R, Kokkos::ALL (), k));
    }

    for (int j = 0; j < NRHS; ++j) {
      cout << "  For right-hand side " << j
           << ": ||B-A*X||_2 = "
           << explicitResidualNorms[j]
           << ", ||B||_2 = " << B_norms[j] << endl;
    }
    cout << endl;
  }
}


template<class DenseMatrixType>
void
findEigenvaluesAndReport (DenseMatrixType A)
{
  using std::cerr;
  using std::cout;
  using std::endl;

  const int numRows = int (A.extent (0));
  const int numCols = int (A.extent (1));
  const int N = std::min (numRows, numCols);

  DenseMatrixType A_copy ("A_copy", numRows, numCols);
  Kokkos::deep_copy (A_copy, A);

  const int LDA = (numRows == 0) ? 1 : int (A_copy.stride (1));

  std::vector<double> realParts (N);
  std::vector<double> imagParts (N);
  std::vector<double> WORK (1);

  int INFO = 0;
  int LWORK = -1; // workspace query
  Teuchos::LAPACK<int, double> lapack;
  lapack.GEEV ('N', 'N', N, A_copy.data (), LDA, realParts.data (),
               imagParts.data (), nullptr, 1, nullptr, 1, WORK.data (),
               LWORK, nullptr, &INFO);
  if (INFO != 0) {
    cerr << "DGELSS returned INFO = " << INFO << " != 0." << endl;
    return;
  }
  LWORK = int (WORK[0]);
  if (LWORK < 0) {
    cerr << "DGEEV reported LWORK = " << LWORK << " < 0." << endl;
    return;
  }
  WORK.resize (LWORK);

  cout << "Solver:" << endl
       << "  Solver type: LAPACK's DGEEV" << endl;
  lapack.GEEV ('N', 'N', N, A_copy.data (), LDA, realParts.data (),
               imagParts.data (), nullptr, 1, nullptr, 1, WORK.data (),
               LWORK, nullptr, &INFO);

  cout << "Results:" << endl
       << "  INFO: " << INFO << endl
       << "  Eigenvalues: ";
  for (int k = 0; k < N; ++k) {
    cout << "(" << realParts[k] << "," << imagParts[k] << ")";
    if (k + 1 < N) {
      cout << ", ";
    }
  }
  cout << endl << endl;
}


template<class SC, class LO, class GO, class NT>
Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
deepCopyFillCompleteCrsMatrix (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
{
  using Teuchos::RCP;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "deepCopyFillCompleteCrsMatrix: Input matrix A must be fillComplete.");
  RCP<crs_matrix_type> A_copy (new crs_matrix_type (A.getCrsGraph ()));
  auto A_copy_lcl = A_copy->getLocalMatrixDevice ();
  auto A_lcl = A.getLocalMatrixDevice ();
  Kokkos::deep_copy (A_copy_lcl.values, A_lcl.values);
  A_copy->fillComplete (A.getDomainMap (), A.getRangeMap ());
  return A_copy;
}

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors =
           ! std::is_same<
               typename Kokkos::ArithTraits<
                 typename ViewType1::non_const_value_type
               >::mag_type,
               typename ViewType2::non_const_value_type
             >::value,
         const int rank = ViewType1::rank>
class ElementWiseMultiply {};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseMultiply<ViewType1,
                          ViewType2,
                          IndexType,
                          takeSquareRootsOfScalingFactors,
                          takeAbsoluteValueOfScalingFactors,
                          1> {
public:
  static_assert (ViewType1::rank == 1, "ViewType1 must be a rank-1 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseMultiply (const ViewType1& X,
                       const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    if (takeAbsoluteValueOfScalingFactors) {
      const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
      const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAM::sqrt (scalFactAbs) : scalFactAbs;
      X_(i) = X_(i) * scalFinalVal;
    }
    else {
      const val_type scalFact = scalingFactors_(i);
      const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAT::sqrt (scalFact) : scalFact;
      X_(i) = X_(i) * scalFinalVal;
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseMultiply<ViewType1,
                          ViewType2,
                          IndexType,
                          takeSquareRootsOfScalingFactors,
                          takeAbsoluteValueOfScalingFactors,
                          2> {
public:
  static_assert (ViewType1::rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseMultiply (const ViewType1& X,
                       const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    for (IndexType j = 0; j < static_cast<IndexType> (X_.extent (1)); ++j) {
      if (takeAbsoluteValueOfScalingFactors) {
        const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
        const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAM::sqrt (scalFactAbs) : scalFactAbs;
        X_(i,j) = X_(i,j) * scalFinalVal;
      }
      else {
        const val_type scalFact = scalingFactors_(i);
        const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAT::sqrt (scalFact) : scalFact;
        X_(i,j) = X_(i,j) * scalFinalVal;
      }
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class MultiVectorViewType,
         class ScalingFactorsViewType,
         class IndexType>
void
elementWiseMultiply (const MultiVectorViewType& X,
                     const ScalingFactorsViewType& scalingFactors,
                     const IndexType numRows,
                     const bool takeSquareRootsOfScalingFactors,
                     const bool takeAbsoluteValueOfScalingFactors =
                       ! std::is_same<
                           typename Kokkos::ArithTraits<
                             typename MultiVectorViewType::non_const_value_type
                           >::mag_type,
                           typename ScalingFactorsViewType::non_const_value_type
                         >::value)
{
  using execution_space = typename MultiVectorViewType::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

template<class MultiVectorType, class ScalingFactorsViewType>
void
elementWiseMultiplyMultiVector (MultiVectorType& X,
                                const ScalingFactorsViewType& scalingFactors,
                                const bool takeSquareRootsOfScalingFactors,
                                const bool takeAbsoluteValueOfScalingFactors =
                                  ! std::is_same<
                                      typename Kokkos::ArithTraits<
                                        typename MultiVectorType::scalar_type
                                      >::mag_type,
                                      typename ScalingFactorsViewType::non_const_value_type
                                    >::value)
{
  using index_type = typename MultiVectorType::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseMultiply (X_lcl_1d, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseMultiply (X_lcl, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
}

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors =
           ! std::is_same<
               typename Kokkos::ArithTraits<
                 typename ViewType1::non_const_value_type
               >::mag_type,
               typename ViewType2::non_const_value_type
             >::value,
         const int rank = ViewType1::rank>
class ElementWiseDivide {};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseDivide<ViewType1,
                        ViewType2,
                        IndexType,
                        takeSquareRootsOfScalingFactors,
                        takeAbsoluteValueOfScalingFactors,
                        1> {
public:
  static_assert (ViewType1::rank == 1, "ViewType1 must be a rank-1 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseDivide (const ViewType1& X,
                     const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    if (takeAbsoluteValueOfScalingFactors) {
      const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
      const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAM::sqrt (scalFactAbs) : scalFactAbs;
      X_(i) = X_(i) / scalFinalVal;
    }
    else {
      const val_type scalFact = scalingFactors_(i);
      const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAT::sqrt (scalFact) : scalFact;
      X_(i) = X_(i) / scalFinalVal;
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class ViewType1,
         class ViewType2,
         class IndexType,
         const bool takeSquareRootsOfScalingFactors,
         const bool takeAbsoluteValueOfScalingFactors>
class ElementWiseDivide<ViewType1,
                        ViewType2,
                        IndexType,
                        takeSquareRootsOfScalingFactors,
                        takeAbsoluteValueOfScalingFactors,
                        2> {
public:
  static_assert (ViewType1::rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseDivide (const ViewType1& X,
                     const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    for (IndexType j = 0; j < static_cast<IndexType> (X_.extent (1)); ++j) {
      if (takeAbsoluteValueOfScalingFactors) {
        const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
        const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAM::sqrt (scalFactAbs) : scalFactAbs;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
      else {
        const val_type scalFact = scalingFactors_(i);
        const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAT::sqrt (scalFact) : scalFact;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

template<class MultiVectorViewType,
         class ScalingFactorsViewType,
         class IndexType>
void
elementWiseDivide (const MultiVectorViewType& X,
                   const ScalingFactorsViewType& scalingFactors,
                   const IndexType numRows,
                   const bool takeSquareRootsOfScalingFactors,
                   const bool takeAbsoluteValueOfScalingFactors =
                     ! std::is_same<
                         typename Kokkos::ArithTraits<
                           typename MultiVectorViewType::non_const_value_type
                         >::mag_type,
                         typename ScalingFactorsViewType::non_const_value_type
                       >::value)
{
  using execution_space = typename MultiVectorViewType::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MultiVectorViewType,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

template<class MultiVectorType, class ScalingFactorsViewType>
void
elementWiseDivideMultiVector (MultiVectorType& X,
                              const ScalingFactorsViewType& scalingFactors,
                              const bool takeSquareRootsOfScalingFactors,
                              const bool takeAbsoluteValueOfScalingFactors =
                                ! std::is_same<
                                    typename Kokkos::ArithTraits<
                                      typename MultiVectorType::scalar_type
                                    >::mag_type,
                                    typename ScalingFactorsViewType::non_const_value_type
                                  >::value)
{
  using index_type = typename MultiVectorType::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseDivide (X_lcl_1d, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseDivide (X_lcl, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
}

std::vector<std::string>
splitIntoStrings (const std::string& s,
                  const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<std::string> strings;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      strings.push_back (stringToUpper (token));
    }
    last_pos = cur_pos + size_type (1);
  }
  return strings;
}

template<class T>
std::vector<T>
splitIntoValues (const std::string& s,
                 const char sep = ',')
{
  using size_type = std::string::size_type;

  size_type cur_pos;
  size_type last_pos = 0;
  size_type length = s.length ();

  std::vector<T> values;
  while (last_pos < length + size_type (1)) {
    cur_pos = s.find_first_of(sep, last_pos);
    if (cur_pos == std::string::npos) {
      cur_pos = length;
    }
    if (cur_pos != last_pos) {
      auto token = std::string (s.data () + last_pos,
                                static_cast<size_type> (cur_pos - last_pos));
      T val {};
      std::istringstream is (token);
      is >> val;
      if (is) {
        values.push_back (val);
      }
    }
    last_pos = cur_pos + size_type (1);
  }
  return values;
}

// Values of command-line arguments.
struct CmdLineArgs {
  std::string matrixFilename;
  std::string rhsFilename;
  std::string solverTypes = "GMRES";
  std::string orthogonalizationMethod = "ICGS";
  std::string convergenceToleranceValues = "1.0e-2";
  std::string maxIterValues = "100";
  std::string restartLengthValues = "20";
  std::string preconditionerTypes = "RELAXATION";
  bool solverVerbose = false;
  bool equilibrate = false;
  bool assumeSymmetric = false;
  bool assumeZeroInitialGuess = true;
  bool useDiagonalToEquilibrate = false;
  bool useLapack = false;
  bool useCustomBicgstab = false;
  bool useAztecOO = false;
  bool printSolution = false;
  bool useReorderingInIfpack2 = false;
};

// Read in values of command-line arguments.
bool
getCmdLineArgs (CmdLineArgs& args, int argc, char* argv[])
{
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("matrixFilename", &args.matrixFilename, "Name of Matrix "
                  "Market file with the sparse matrix A");
  cmdp.setOption ("rhsFilename", &args.rhsFilename, "Name of Matrix Market "
                  "file with the right-hand side vector(s) B");
  cmdp.setOption ("solverTypes", &args.solverTypes,
                  "One or more Belos solver types, "
                  "separated by commas");
  cmdp.setOption ("convergenceTolerances", &args.convergenceToleranceValues,
                  "One or more doubles, separated by commas; each value "
                  "is a convergence tolerance to try");
  cmdp.setOption ("orthogonalizationMethod", &args.orthogonalizationMethod,
                  "Orthogonalization method (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("maxIters", &args.maxIterValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum number of solver iterations to try");
  cmdp.setOption ("restartLengths", &args.restartLengthValues,
                  "One or more integers, separated by commas; each value "
                  "is a maximum restart length to try (for GMRES solver only; "
                  "ignored otherwise)");
  cmdp.setOption ("preconditionerTypes", &args.preconditionerTypes,
                  "One or more Ifpack2 preconditioner types, "
                  "separated by commas");
  cmdp.setOption ("solverVerbose", "solverQuiet", &args.solverVerbose,
                  "Whether the Belos solver should print verbose output");
  cmdp.setOption ("equilibrate", "no-equilibrate", &args.equilibrate,
                  "Whether to equilibrate the linear system before solving it");
  cmdp.setOption ("assumeSymmetric", "no-assumeSymmetric",
                  &args.assumeSymmetric, "Whether equilibration should assume "
                  "that the matrix is symmetric");
  cmdp.setOption ("assumeZeroInitialGuess", "assumeNonzeroInitialGuess",
                  &args.assumeZeroInitialGuess, "Whether equilibration should "
                  "assume that the initial guess (vector) is zero");
  cmdp.setOption ("useDiagonalToEquilibrate", "useOneNorms",
                  &args.useDiagonalToEquilibrate,
                  "Whether equilibration should use the matrix's diagonal; "
                  "default is to use row and column one norms");
  cmdp.setOption ("lapack", "no-lapack",
                  &args.useLapack,
                  "Whether to compare against LAPACK's LU factorization "
                  "(expert driver).  If --equilibration, then use "
                  "equilibration in LAPACK too.");
  cmdp.setOption ("custom-bicgstab", "no-custom-bicgstab",
                  &args.useCustomBicgstab,
                  "Whether to compare against a hand-rolled BiCGSTAB "
                  "implementation.");
  cmdp.setOption ("aztecoo", "no-aztecoo",
                  &args.useAztecOO,
                  "Whether to compare against AztecOO.");
  cmdp.setOption ("printSolution", "dontPrintSolution",
                  &args.printSolution,
                  "Whether to print the solution vector(s), "
                  "for each solution method.");
  cmdp.setOption ("useReorderingInIfpack2", "dontUseReorderingInIfpack2",
                  &args.useReorderingInIfpack2,
                  "Whether to use reordering in Ifpack2.");

  auto result = cmdp.parse (argc, argv);
  return result == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL;
}

template<class ScalarType>
struct BelosSolverResult {
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType achievedTolerance;
  int numIters;
  bool converged;
  bool lossOfAccuracy;
};

template<class CrsMatrixType>
class BelosIfpack2Solver {
private:
  using scalar_type = typename CrsMatrixType::scalar_type;
  using local_ordinal_type = typename CrsMatrixType::local_ordinal_type;
  using global_ordinal_type = typename CrsMatrixType::global_ordinal_type;
  using node_type = typename CrsMatrixType::node_type;
  using row_matrix_type =
    Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using OP = Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
  using problem_type = Belos::LinearProblem<scalar_type, MV, OP>;
  using solver_type = Belos::SolverManager<scalar_type, MV, OP>;
  using preconditioner_type =
    Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;

  void createPreconditioner ()
  {
    rightPrec_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    if (precType_ != "NONE") {
      rightPrec_ =
        Ifpack2::Factory::create<row_matrix_type> (precType_, A_);
      rightPrec_->setParameters (* (precParams_));
    }
  }

  void createSolver ()
  {
    Belos::SolverFactory<scalar_type, MV, OP> belosFactory;
    solver_ = Teuchos::null; // destroy old one first, to reduce peak memory use
    solver_ = belosFactory.create (solverType_, solverParams_);
  }

  void setPreconditionerMatrix (const Teuchos::RCP<const row_matrix_type>& A)
  {
    if (precType_ != "NONE" && rightPrec_.get () != nullptr) {
      // Not all Ifpack2 preconditioners have a setMatrix method.  Use
      // it if it exists; otherwise, start over with a new instance.
      using can_change_matrix = Ifpack2::Details::CanChangeMatrix<row_matrix_type>;
      can_change_matrix* rightPrec = dynamic_cast<can_change_matrix*> (rightPrec_.get ());
      if (rightPrec != nullptr) {
        rightPrec_->setMatrix (A);
      }
      else {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
    }
  }

  void initializePreconditioner ()
  {
    if (precType_ != "NONE") {
      if (rightPrec_.get () == nullptr) {
        createPreconditioner ();
      }
      rightPrec_->initialize ();
    }
  }

  void computePreconditioner ()
  {
    if (precType_ != "NONE") {
      if (rightPrec_.get () == nullptr) {
        createPreconditioner ();
        rightPrec_->initialize ();
      }
      rightPrec_->compute ();
    }
  }

  void equilibrateMatrix ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call compute.");
    if (equilibrate_) {
      using Tpetra::computeRowAndColumnOneNorms;
      using Tpetra::leftAndOrRightScaleCrsMatrix;

      equibResult_ = computeRowAndColumnOneNorms (*A_, assumeSymmetric_);
      if (useDiagonalToEquilibrate_) {
        using device_type = typename node_type::device_type;
        using mag_type = typename Kokkos::ArithTraits<scalar_type>::mag_type;
        using view_type = Kokkos::View<mag_type*, device_type>;

        view_type rowDiagAbsVals ("rowDiagAbsVals",
                                  equibResult_.rowDiagonalEntries.extent (0));
        KokkosBlas::abs (rowDiagAbsVals, equibResult_.rowDiagonalEntries);
        view_type colDiagAbsVals ("colDiagAbsVals",
                                  equibResult_.colDiagonalEntries.extent (0));
        KokkosBlas::abs (colDiagAbsVals, equibResult_.colDiagonalEntries);

        leftAndOrRightScaleCrsMatrix (*A_, rowDiagAbsVals, colDiagAbsVals,
                                      true, true, equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        leftAndOrRightScaleCrsMatrix (*A_, equibResult_.rowNorms,
                                      colScalingFactors, true, true,
                                      equibResult_.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
    } // if equilibrate_
  }

  void
  preScaleRightHandSides (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type>& B) const
  {
    if (equilibrate_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (B, equibResult_.rowDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        const bool takeSquareRootsOfScalingFactors = equibResult_.assumeSymmetric;
        elementWiseDivideMultiVector (B, equibResult_.rowNorms,
                                      takeSquareRootsOfScalingFactors);
      }
    }
  }

  void
  preScaleInitialGuesses (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type>& X) const
  {
    if (equilibrate_ && ! assumeZeroInitialGuess_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseMultiplyMultiVector (X, equibResult_.colDiagonalEntries,
                                        takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult_.assumeSymmetric;
        elementWiseMultiplyMultiVector (X, colScalingFactors,
                                        takeSquareRootsOfScalingFactors);
      }
    }
  }

  void
  postScaleSolutionVectors (Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type>& X) const
  {
    if (equilibrate_) {
      if (useDiagonalToEquilibrate_) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (X, equibResult_.colDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult_.assumeSymmetric ?
          equibResult_.colNorms :
          equibResult_.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult_.assumeSymmetric;
        elementWiseDivideMultiVector (X, colScalingFactors,
                                      takeSquareRootsOfScalingFactors);
      }
    }
  }

public:
  BelosIfpack2Solver () = default;

  BelosIfpack2Solver (const Teuchos::RCP<CrsMatrixType>& A,
                      const std::string& solverType = "GMRES",
                      const std::string& precType = "NONE") :
    A_ (A),
    solverType_ (solverType),
    precType_ (precType),
    equilibrate_ (false),
    useDiagonalToEquilibrate_ (false)
  {}

  void setMatrix (const Teuchos::RCP<const CrsMatrixType>& A)
  {
    if (A_.get () != A.get ()) {
      setPreconditionerMatrix (A);
      // Belos solvers don't deal well with a complete change of the matrix.
      solver_ = Teuchos::null;
    }
    A_ = A;
  }

  void
  setPreconditionerTypeAndParameters (const std::string& precType,
                                      const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    precType_ = precType;
    precParams_ = params;
    if (rightPrec_.get () != nullptr) {
      if (precType_ != precType) {
        rightPrec_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        rightPrec_->setParameters (*params);
      }
    }
  }

  void
  setSolverTypeAndParameters (const std::string& solverType,
                              const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    solverType_ = solverType;
    solverParams_ = params;
    if (solver_.get () != nullptr) {
      if (solverType_ != solverType) {
        solver_ = Teuchos::null; // blow it away; make a new one only on demand
      }
      else {
        solver_->setParameters (params);
      }
    }
  }

  void
  setEquilibrationParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    if (params.get () != nullptr) {
      equilibrate_ = params->get ("Equilibrate", equilibrate_);
      assumeSymmetric_ = params->get ("Assume symmetric", assumeSymmetric_);
      assumeZeroInitialGuess_ = params->get ("Assume zero initial guess",
                                             assumeZeroInitialGuess_);
      useDiagonalToEquilibrate_ = params->get ("Use diagonal to equilibrate",
                                               useDiagonalToEquilibrate_);
    }
  }

  void initialize ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call initialize.");
    // Calling this implies that the matrix's graph has changed.
    // Belos' solvers don't handle that very well, so best practice is
    // to recreate them in this case.
    solver_ = Teuchos::null;
    initializePreconditioner ();
  }

  void compute ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.get () == nullptr, std::runtime_error, "Solver: You must call "
       "setMatrix with a nonnull matrix before you may call compute.");
    equilibrateMatrix ();
    // equilibration changes the matrix, so don't compute the
    // preconditioner until after doing that.
    computePreconditioner ();
  }

  BelosSolverResult<scalar_type>
  solve (Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
         Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& B)
  {
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;

    if (solver_.get () == nullptr) {
      createSolver ();
    }
    if (rightPrec_.get () == nullptr) {
      createPreconditioner ();
      initializePreconditioner ();
      computePreconditioner ();
    }

    preScaleRightHandSides (B);
    preScaleInitialGuesses (X);

    RCP<problem_type> problem (new problem_type (A_, rcpFromRef (X), rcpFromRef (B)));
    if (rightPrec_.get () != nullptr) {
      problem->setRightPrec (rightPrec_);
    }
    problem->setProblem ();
    solver_->setProblem (problem);
    const Belos::ReturnType solveResult = solver_->solve ();

    postScaleSolutionVectors (X);

    typename Teuchos::ScalarTraits<scalar_type>::magnitudeType tol {-1.0};
    try { // not every Belos solver implements this
      tol = solver_->achievedTol ();
    }
    catch (...) {}
    const int numIters = solver_->getNumIters ();
    const bool converged = (solveResult == Belos::Converged);
    const bool lossOfAccuracy = solver_->isLOADetected ();
    return {tol, numIters, converged, lossOfAccuracy};
  }

private:
  Teuchos::RCP<CrsMatrixType> A_;
  Teuchos::RCP<solver_type> solver_;
  Teuchos::RCP<preconditioner_type> rightPrec_;
  using equilibration_result_type = decltype (Tpetra::computeRowAndColumnOneNorms (*A_, false));
  equilibration_result_type equibResult_;

  std::string solverType_;
  Teuchos::RCP<Teuchos::ParameterList> solverParams_;
  std::string precType_;
  Teuchos::RCP<Teuchos::ParameterList> precParams_;

  bool equilibrate_;
  bool assumeSymmetric_;
  bool assumeZeroInitialGuess_;
  bool useDiagonalToEquilibrate_;
};

class TpetraInstance {
public:
  TpetraInstance (int* argc, char*** argv) {
    Tpetra::initialize (argc, argv);
  }

  ~TpetraInstance () {
    Tpetra::finalize ();
  }
};

template<class CrsMatrixType, class MultiVectorType>
void
solveAndReport (BelosIfpack2Solver<CrsMatrixType>& solver,
                const CrsMatrixType& A_original, // before scaling
                MultiVectorType& X,
                MultiVectorType& B,
                const std::string& solverType,
                /*const std::string& */ std::string precType,
                const typename MultiVectorType::mag_type convergenceTolerance,
                const int maxIters,
                const int restartLength,
                const CmdLineArgs& args)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  X.putScalar (0.0);

  RCP<ParameterList> solverParams (new ParameterList ("Belos"));
  if (args.solverVerbose) {
    solverParams->set ("Verbosity",
                       Belos::IterationDetails |
                       Belos::FinalSummary |
                       Belos::StatusTestDetails);
  }
  solverParams->set ("Convergence Tolerance", convergenceTolerance);
  solverParams->set ("Maximum Iterations", maxIters);
  if (solverType == "GMRES") {
    solverParams->set ("Num Blocks", restartLength);
    solverParams->set ("Maximum Restarts", restartLength * maxIters);
    solverParams->set ("Orthogonalization", args.orthogonalizationMethod);
  }

  RCP<ParameterList> precParams (new ParameterList ("Ifpack2"));

  if (args.useReorderingInIfpack2) {
    precType = "SCHWARZ";
    precParams->set ("schwarz: subdomain solver name", "RELAXATION");
    {
      ParameterList& relaxParams =
        precParams->sublist ("schwarz: subdomain solver parameters", false);
      relaxParams.set ("relaxation: type", "Symmetric Gauss-Seidel");
    }
    precParams->set ("schwarz: overlap level", int (0));
    precParams->set ("schwarz: use reordering", true);
    //precParams->set ("schwarz: filter singletons", true);
  }
  else {
    if (precType == "RELAXATION") {
      precParams->set ("relaxation: type", "Symmetric Gauss-Seidel");
    }
  }

  RCP<ParameterList> equibParams (new ParameterList ("Equilibration"));
  equibParams->set ("Equilibrate", args.equilibrate);
  equibParams->set ("Assume symmetric", args.assumeSymmetric);
  equibParams->set ("Assume zero initial guess",
                    args.assumeZeroInitialGuess);
  equibParams->set ("Use diagonal to equilibrate",
                    args.useDiagonalToEquilibrate);

  solver.setSolverTypeAndParameters (solverType, solverParams);
  solver.setPreconditionerTypeAndParameters (precType, precParams);
  solver.setEquilibrationParameters (equibParams);

  solver.initialize ();
  solver.compute ();

  // Keep this around for later computation of the explicit residual
  // norm.  If the solver equilibrates, it will modify the original B.
  MultiVectorType R (B, Teuchos::Copy);

  // Compute ||B||_2.
  using mag_type = typename MultiVectorType::mag_type;
  Teuchos::Array<mag_type> norms (R.getNumVectors ());
  R.norm2 (norms ());
  mag_type B_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < B.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    B_norm2_max = norms[j] < B_norm2_max ? B_norm2_max : norms[j];
  }

  // Solve the linear system AX=B.
  auto result = solver.solve (X, B);

  // Compute the actual residual norm ||B - A*X||_2.
  using scalar_type = typename MultiVectorType::scalar_type;
  const scalar_type ONE = Teuchos::ScalarTraits<scalar_type>::one ();
  A_original.apply (X, R, Teuchos::NO_TRANS, -ONE, ONE); // R := -A*X + B
  R.norm2 (norms ());

  mag_type R_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < R.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    R_norm2_max = norms[j] < R_norm2_max ? R_norm2_max : norms[j];
  }

  X.norm2 (norms ());
  mag_type X_norm2_max = Kokkos::ArithTraits<mag_type>::zero ();
  for (std::size_t j = 0; j < R.getNumVectors (); ++j) {
    // Any NaN will persist (since the first test will fail);
    // this is what we want
    X_norm2_max = norms[j] < X_norm2_max ? X_norm2_max : norms[j];
  }

  const int myRank = X.getMap ()->getComm ()->getRank ();
  if (myRank == 0) {
    using std::cout;
    using std::endl;
    cout << "Solver:" << endl
         << "  Solver type: " << solverType << endl
         << "  Preconditioner type: " << precType << endl
         << "  Convergence tolerance: " << convergenceTolerance << endl
         << "  Maximum number of iterations: " << maxIters << endl;
    if (solverType == "GMRES") {
      cout << "  Restart length: " << restartLength << endl
           << "  Orthogonalization method: " << args.orthogonalizationMethod << endl;
    }
    cout << "Results:" << endl
         << "  Converged: " << (result.converged ? "true" : "false") << endl
         << "  Number of iterations: " << result.numIters << endl
         << "  Achieved tolerance: " << result.achievedTolerance << endl
         << "  Loss of accuracy: " << result.lossOfAccuracy << endl
         << "  ||B-A*X||_2: " << R_norm2_max << endl
         << "  ||B||_2: " << B_norm2_max << endl
         << "  ||X||_2: " << X_norm2_max << endl;
    if (B_norm2_max != Kokkos::ArithTraits<mag_type>::zero ()) {
      cout << "  ||B-A*X||_2 / ||B||_2: " << (R_norm2_max / B_norm2_max)
           << endl;
    }
    cout << endl;
  }

  if (args.printSolution) {
    using writer_type = Tpetra::MatrixMarket::Writer<CrsMatrixType>;
    writer_type::writeDense (std::cout, X, "X", "Belos solution");
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::cerr;
  using std::endl;
  using crs_matrix_type = Tpetra::CrsMatrix<>;
  using MV = Tpetra::MultiVector<>;
  // using mag_type = MV::mag_type;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;
  using writer_type = Tpetra::MatrixMarket::Writer<crs_matrix_type>;

  TpetraInstance tpetraInstance (&argc, &argv);
  auto comm = Tpetra::getDefaultComm ();

  // Get command-line arguments.
  CmdLineArgs args;
  const bool gotCmdLineArgs = getCmdLineArgs (args, argc, argv);
  if (! gotCmdLineArgs) {
    if (comm->getRank () == 0) {
      cerr << "Failed to get command-line arguments!" << endl;
    }
    return EXIT_FAILURE;
  }

#if ! defined(HAVE_IFPACK2_AZTECOO) || ! defined(HAVE_IFPACK2_EPETRA)
  if (args.useAztecOO) {
    if (comm->getRank () == 0) {
      cerr << "In order to solve AztecOO, you must have built with the AztecOO "
        "and Epetra packages enabled.  Please reconfigure Trilinos (i.e., "
        "rerun CMake) with the CMake options Trilinos_ENABLE_AztecOO:BOOL=ON "
        "and Trilinos_ENABLE_Epetra:BOOL=ON set." << endl;
    }
    return EXIT_FAILURE;
  }
#endif

  if (args.matrixFilename == "") {
    if (comm->getRank () == 0) {
      cerr << "Must specify sparse matrix filename!" << endl;
    }
    return EXIT_FAILURE;
  }

  std::vector<std::string> solverTypes;
  if (args.solverTypes == "") {
    solverTypes = {"GMRES", "TFQMR", "BICGSTAB"};
  }
  else {
    solverTypes = splitIntoStrings (args.solverTypes);
  }

  std::vector<std::string> preconditionerTypes;
  if (args.preconditionerTypes == "") {
    preconditionerTypes = {"RELAXATION"};
  }
  else {
    preconditionerTypes = splitIntoStrings (args.preconditionerTypes);
  }

  std::vector<int> maxIterValues;
  if (args.maxIterValues == "") {
    maxIterValues = {100};
  }
  else {
    maxIterValues = splitIntoValues<int> (args.maxIterValues);
  }

  std::vector<int> restartLengthValues;
  if (args.restartLengthValues == "") {
    restartLengthValues = {20};
  }
  else {
    restartLengthValues = splitIntoValues<int> (args.restartLengthValues);
  }

  std::vector<double> convergenceToleranceValues;
  if (args.convergenceToleranceValues == "") {
    convergenceToleranceValues = {20};
  }
  else {
    convergenceToleranceValues =
      splitIntoValues<double> (args.convergenceToleranceValues);
  }

  // Read sparse matrix A from Matrix Market file.
  RCP<crs_matrix_type> A =
    reader_type::readSparseFile (args.matrixFilename, comm);
  if (A.get () == nullptr) {
    if (comm->getRank () == 0) {
      cerr << "Failed to load sparse matrix A from file "
        "\"" << args.matrixFilename << "\"!" << endl;
    }
    return EXIT_FAILURE;
  }

  // Read right-hand side vector(s) B from Matrix Market file, or
  // generate B if file not specified.
  RCP<MV> B;
  RCP<MV> X;
  if (args.rhsFilename == "") {
    B = Teuchos::rcp (new MV (A->getRangeMap (), 1));
    X = Teuchos::rcp (new MV (A->getDomainMap (), 1));
    X->putScalar (1.0);
    A->apply (*X, *B);
    X->putScalar (0.0);

    double norm {0.0};
    using host_device_type = Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>;
    Kokkos::View<double*, host_device_type> normView (&norm, B->getNumVectors ());
    B->norm2 (normView);
    if (norm != 0.0) {
      B->scale (1.0 / norm);
    }
  }
  else {
    auto map = A->getRangeMap ();
    B = reader_type::readDenseFile (args.rhsFilename, comm, map);
    if (B.get () == nullptr) {
      if (comm->getRank () == 0) {
        cerr << "Failed to load right-hand side vector(s) from file \""
             << args.rhsFilename << "\"!" << endl;
      }
      return EXIT_FAILURE;
    }
    X = Teuchos::rcp (new MV (A->getDomainMap (), B->getNumVectors ()));
  }

  auto A_original = deepCopyFillCompleteCrsMatrix (*A);

  if (args.useLapack) {
    using std::cout;

    using lapack_matrix_type = HostDenseMatrix<>;
    lapack_matrix_type A_lapack;
    lapack_matrix_type B_lapack;
    int errCode = gatherAndDensify (A_lapack, B_lapack, *A, *B);

    int lapackSolveOk = (errCode == 0) ? 1 : 0;
    Teuchos::broadcast (*comm, 0, Teuchos::inOutArg (lapackSolveOk));
    if (lapackSolveOk != 1) {
      if (comm->getRank () == 0) {
        cerr << "Gathering and densification for LAPACK solve FAILED!" << endl;
      }
    }

    if (lapackSolveOk && comm->getRank () == 0) {
      const int numRows = int (A_lapack.extent (0));
      const int numCols = int (A_lapack.extent (1));
      const int LDA = (numRows == 0) ? 1 : A_lapack.stride (1);
      const int NRHS = int (B_lapack.extent (1));
      const int LDB = (numRows == 0) ? 1 : B_lapack.stride (1);

      // cout << "The matrix A, in dense format:" << endl;
      // for (int i = 0; i < numRows; ++i) {
      //        for (int j = 0; j < numCols; ++j) {
      //          cout << A_lapack(i,j);
      //          if (j + 1 < numCols) {
      //            cout << " ";
      //          }
      //        }
      //        cout << endl;
      // }
      // cout << endl;

      // cout << "The right-hand side(s) B, in dense format:" << endl;
      // for (int i = 0; i < numRows; ++i) {
      //        for (int j = 0; j < NRHS; ++j) {
      //          cout << B_lapack(i,j);
      //          if (j + 1 < NRHS) {
      //            cout << " ";
      //          }
      //        }
      //        cout << endl;
      // }
      // cout << endl;

      if (numRows != numCols) {
        cerr << "Matrix is not square, so we can't use LAPACK's LU factorization on it!" << endl;
        lapackSolveOk = 0;
      }
      else {
        findEigenvaluesAndReport (A_lapack);
        solveLeastSquaresProblemAndReport (A_lapack, B_lapack);

        Teuchos::LAPACK<int, double> lapack;

        int INFO = 0;

        lapack_matrix_type AF_lapack ("AF", numRows, numCols);
        lapack_matrix_type X_lapack ("X_lapack", numCols, NRHS);
        const int LDAF = (numRows == 0) ? 1 : AF_lapack.stride (1);
        const int LDX = (numRows == 0) ? 1 : X_lapack.stride (1);

        // Save norms of columns of B, since _GESVX may modify B.
        std::vector<double> B_norms (NRHS);
        for (int k = 0; k < NRHS; ++k) {
          B_norms[k] =
            KokkosBlas::nrm2 (Kokkos::subview (B_lapack, Kokkos::ALL (), k));
        }

        // Save a copy of A for later.
        lapack_matrix_type A_copy ("A_copy", numRows, numCols);
        Kokkos::deep_copy (A_copy, A_lapack);

        // Save a copy of B for later.
        lapack_matrix_type R_lapack ("R_lapack", numRows, NRHS);
        Kokkos::deep_copy (R_lapack, B_lapack);

        std::vector<int> IPIV (numCols);
        std::vector<double> R (numRows);
        std::vector<double> C (numCols);
        std::vector<double> FERR (NRHS);
        std::vector<double> BERR (NRHS);
        std::vector<double> WORK (4*numRows);
        std::vector<int> IWORK (numRows);
        INFO = 0;

        const char FACT ('E');
        const char TRANS ('N');
        char EQUED[1];
        EQUED[0] = args.equilibrate ? 'B' : 'N';
        double RCOND = 1.0;

        cout << "Solver:" << endl
             << "  Solver type: LAPACK's _GESVX" << endl
             << "  Equilibrate: " << (args.equilibrate ? "YES" : "NO") << endl;
        lapack.GESVX (FACT, TRANS, numRows, NRHS, A_lapack.data (), LDA,
                      AF_lapack.data (), LDAF, IPIV.data (), EQUED, R.data (),
                      C.data (), B_lapack.data (), LDB, X_lapack.data (), LDX,
                      &RCOND, FERR.data (), BERR.data (), WORK.data (),
                      IWORK.data (), &INFO);

        cout << "Results:" << endl
             << "  INFO: " << INFO << endl
             << "  RCOND: " << RCOND << endl;
        if (NRHS > 0) {
          cout << "  Pivot growth factor: " << WORK[0] << endl;
        }
        for (int j = 0; j < NRHS; ++j) {
          cout << "  For right-hand side " << j
               << ": forward error is " << FERR[j]
               << ", backward error is " << BERR[j] << endl;
        }

        // Compute the explicit residual norm(s).
        KokkosBlas::gemm ("N", "N", -1.0, A_copy, X_lapack, +1.0, R_lapack);
        std::vector<double> explicitResidualNorms (NRHS);
        for (int k = 0; k < NRHS; ++k) {
          explicitResidualNorms[k] =
            KokkosBlas::nrm2 (Kokkos::subview (R_lapack, Kokkos::ALL (), k));
        }

        for (int j = 0; j < NRHS; ++j) {
          cout << "  For right-hand side " << j
               << ": ||B-A*X||_2 = "
               << explicitResidualNorms[j]
               << ", ||B||_2 = " << B_norms[j] << endl;
        }
        cout << endl;
      }
    }

    Teuchos::broadcast (*comm, 0, Teuchos::inOutArg (lapackSolveOk));
    if (lapackSolveOk != 1) {
      if (comm->getRank () == 0) {
        cerr << "LAPACK solve FAILED!" << endl;
      }
    }
  }

  if (args.useCustomBicgstab) {
    for (double convTol : convergenceToleranceValues) {
      for (int maxIters : maxIterValues) {
        MV X_copy (*X, Teuchos::Copy);
        MV B_copy (*B, Teuchos::Copy);
        const Tpetra::Operator<>* M = nullptr;
        const Tpetra::Operator<>& A_ref = static_cast<const Tpetra::Operator<>& > (*A);
        auto result = bicgstab_aztecoo (X_copy, A_ref, M, B_copy, maxIters, convTol);

        MV R (*B, Teuchos::Copy);
        A->apply (X_copy, R, Teuchos::NO_TRANS, -1.0, 1.0);
        const double B_norm = norm (B_copy);
        const double R_norm = norm (R);
        const double explicitRelResNorm = R_norm / B_norm;

        if (comm->getRank () == 0) {
          using std::cout;
          using std::endl;
          cout << "Solver:" << endl
               << "  Solver type: AztecOO-imitating custom BiCGSTAB" << endl
               << "  Preconditioner type: NONE" << endl
               << "  Convergence tolerance: " << convTol << endl
               << "  Maximum number of iterations: " << maxIters << endl
               << "Results:" << endl
               << "  Converged: " << std::get<2> (result) << endl
               << "  Number of iterations: " << std::get<1> (result) << endl
               << "  Achieved tolerance: " << std::get<0> (result) << endl
               << "  ||B-A*X||_2 / ||B||_2: " << explicitRelResNorm << endl
               << endl;
        }
        if (args.printSolution) {
          writer_type::writeDense (std::cout, X_copy, "X", "AztecOO-imitating custom BiCGSTAB solution");
        }
      }
    }

    for (double convTol : convergenceToleranceValues) {
      for (int maxIters : maxIterValues) {
        MV X_copy (*X, Teuchos::Copy);
        MV B_copy (*B, Teuchos::Copy);
        const Tpetra::Operator<>& A_ref = static_cast<const Tpetra::Operator<>& > (*A);
        auto result = bicgstab_no_prec_paper (X_copy, A_ref, B_copy, maxIters, convTol);

        MV R (*B, Teuchos::Copy);
        A->apply (X_copy, R, Teuchos::NO_TRANS, -1.0, 1.0);
        const double B_norm = norm (B_copy);
        const double R_norm = norm (R);
        const double explicitRelResNorm = R_norm / B_norm;

        if (comm->getRank () == 0) {
          using std::cout;
          using std::endl;
          cout << "Solver:" << endl
               << "  Solver type: Original paper BiCGSTAB" << endl
               << "  Preconditioner type: NONE" << endl
               << "  Convergence tolerance: " << convTol << endl
               << "  Maximum number of iterations: " << maxIters << endl
               << "Results:" << endl
               << "  Converged: " << std::get<2> (result) << endl
               << "  Number of iterations: " << std::get<1> (result) << endl
               << "  Achieved tolerance: " << std::get<0> (result) << endl
               << "  ||B-A*X||_2 / ||B||_2: " << explicitRelResNorm << endl
               << endl;
        }
        if (args.printSolution) {
          writer_type::writeDense (std::cout, X_copy, "X", "Original paper BiCGSTAB solution");
        }
      }
    }
  }

#if defined(HAVE_IFPACK2_AZTECOO) && defined(HAVE_IFPACK2_EPETRA)
  if (args.useAztecOO) {
    MV X_copy (*X, Teuchos::Copy);
    MV B_copy (*B, Teuchos::Copy);
    auto epetraLinSys = tpetraToEpetraLinearSystem (*A, X_copy, B_copy);

    for (std::string solverType : solverTypes) {
      for (double convTol : convergenceToleranceValues) {
        for (int maxIters : maxIterValues) {
          auto result =
            solveEpetraLinearSystemWithAztecOO (epetraLinSys, solverType,
                                                convTol, maxIters);
          const auto& X_epetra = std::get<1> (epetraLinSys);
          deep_copy (X_copy, X_epetra);

          MV R (B_copy, Teuchos::Copy);
          A->apply (X_copy, R, Teuchos::NO_TRANS, -1.0, 1.0);
          const double B_norm = norm (B_copy);
          const double R_norm = norm (R);
          const double explicitRelResNorm = R_norm / B_norm;

          if (comm->getRank () == 0) {
            using std::cout;
            using std::endl;
            cout << "Solver:" << endl
                 << "  Solver type: AztecOO " << solverType << endl
                 << "  Preconditioner type: NONE" << endl
                 << "  Convergence tolerance: " << convTol << endl
                 << "  Maximum number of iterations: " << maxIters << endl
                 << "Results:" << endl
                 << "  Converged: " << std::get<2> (result) << endl
                 << "  Number of iterations: " << std::get<1> (result) << endl
                 << "  Achieved tolerance: " << std::get<0> (result) << endl
                 << "  ||B-A*X||_2 / ||B||_2: " << explicitRelResNorm << endl
                 << endl;
          }
          if (args.printSolution) {
            writer_type::writeDense (std::cout, X_copy, "X", "AztecOO solution");
          }
        }
      }
    }
  }
#endif // defined(HAVE_IFPACK2_AZTECOO) && defined(HAVE_IFPACK2_EPETRA)

  // Create the solver.
  BelosIfpack2Solver<crs_matrix_type> solver (A);

  // Solve the linear system using various solvers and preconditioners.
  for (std::string solverType : solverTypes) {
    for (std::string precType : preconditionerTypes) {
      for (int maxIters : maxIterValues) {
        for (int restartLength : restartLengthValues) {
          for (double convTol : convergenceToleranceValues) {
            solveAndReport (solver, *A_original, *X, *B,
                            solverType,
                            precType,
                            convTol,
                            maxIters,
                            restartLength,
                            args);
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
