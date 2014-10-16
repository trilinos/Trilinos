// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef BELOS_TPETRA_ADAPTER_MP_VECTOR_HPP
#define BELOS_TPETRA_ADAPTER_MP_VECTOR_HPP

#include <Kokkos_NodeTrace.hpp>

/*! \file BelosTpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Tpetra concrete classes.
*/

// TODO: the assumption is made that the solver, multivector and operator are templated on the same scalar. this will need to be modified.

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

#ifdef HAVE_BELOS_TSQR
#  include <Tpetra_TsqrAdaptor_MP_Vector.hpp>
#endif // HAVE_BELOS_TSQR


namespace Belos {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of Belos::MultiVecTraits for Tpetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of MultiVecTraits for MV = Tpetra::MultiVector.
  ///
  /// This interface lets Belos' solvers work directly with
  /// Tpetra::MultiVector objects as the multivector type
  /// (corresponding to the MV template parameter).
  ///
  /// The four template parameters of this partial specialization
  /// correspond exactly to the four template parameters of
  /// Tpetra::MultiVector.  See the Tpetra::MultiVector documentation
  /// for more information.
  template<class BaseScalar, class Storage, class LO, class GO, class Node>
  class MultiVecTraits<BaseScalar,
                       Tpetra::MultiVector< Sacado::MP::Vector<Storage>,
                                            LO,GO,Node> > {
    typedef typename Storage::ordinal_type s_ordinal;
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:

    /// \brief Create a new multivector with \c numvecs columns.
    ///
    /// The returned Tpetra::MultiVector has the same Tpetra::Map
    /// (distribution over one or more parallel processes) as \c mv.
    /// Its entries are not initialized and have undefined values.
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    Clone (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs)
    {
      return Teuchos::rcp (new MV (mv.getMap (), numvecs, false));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& X)
    {
      // Make a deep copy of X.  The one-argument copy constructor
      // does a shallow copy by default; the second argument tells it
      // to do a deep copy.
      Teuchos::RCP<MV> X_copy (new MV (X, Teuchos::Copy));
      // Make Tpetra::MultiVector use the new view semantics.  This is
      // a no-op for the Kokkos refactor version of Tpetra; it only
      // does something for the "classic" version of Tpetra.  This
      // shouldn't matter because Belos only handles MV through RCP
      // and through this interface anyway, but it doesn't hurt to set
      // it and make sure that it works.
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits::CloneCopy(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::runtime_error, fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors () << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subCopy, so we
      // don't have to check here.
      Teuchos::RCP<MV> X_copy = mv.subCopy (columns ());
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < GetNumberVecs(mv);
      if (! validRange) { // invalid range; generate error message
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneCopy(mv,index=["
           << index.lbound() << "," << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= GetNumberVecs(mv), std::invalid_argument,
          os.str() << "Index range exceeds number of vectors "
          << mv.getNumVectors() << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_copy = mv.subCopy (index);
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }


    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits::CloneViewNonConst(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in
      // MultiVector::subViewNonConst, so we don't have to check here.
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (columns ());
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }


    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneViewNonConst(mv,index=["
           << index.lbound() << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative inde{x,ices}, which is "
          "not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }


    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar, "
        "Tpetra::MultiVector<...> >::CloneView(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subView, so we
      // don't have to check here.
      Teuchos::RCP<const MV> X_view = mv.subView (columns);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneView(mv, index=["
           << index.lbound () << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<const MV> X_view = mv.subView (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static int
    GetVecLength (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.getGlobalLength ();
    }

    static int
    GetNumberVecs (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.getNumVectors ();
    }

    static bool
    HasConstantStride (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.isConstantStride ();
    }

    static void
    MvTimesMatAddMv (const Scalar& alpha,
                     const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
                     const Teuchos::SerialDenseMatrix<int,BaseScalar>& B,
                     const Scalar& beta,
                     Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      Teuchos::SerialDenseMatrix<int,Scalar> B_mp(B.numRows(), B.numCols());
      for (int i=0; i<B.numRows(); i++)
        for (int j=0; j<B.numCols(); j++)
          B_mp(i,j) = B(i,j);
      MvTimesMatAddMv(alpha, A, B_mp, beta, mv);
    }

    static void
    MvTimesMatAddMv (const Scalar& alpha,
                     const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
                     const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                     const Scalar& beta,
                     Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      using Teuchos::ArrayView;
      using Teuchos::Comm;
      using Teuchos::rcpFromRef;
      typedef Tpetra::Map<LO, GO, Node> map_type;

#ifdef HAVE_BELOS_TPETRA_TIMERS
      const std::string timerName ("Belos::MVT::MvTimesMatAddMv");
      Teuchos::RCP<Teuchos::Time> timer =
        Teuchos::TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = Teuchos::TimeMonitor::getNewCounter (timerName);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        timer.is_null (), std::logic_error,
        "Belos::MultiVecTraits::MvTimesMatAddMv: "
        "Failed to look up timer \"" << timerName << "\".  "
        "Please report this bug to the Belos developers.");

      // This starts the timer.  It will be stopped on scope exit.
      Teuchos::TimeMonitor timeMon (*timer);
#endif // HAVE_BELOS_TPETRA_TIMERS

      // Check if B is 1-by-1, in which case we can just call update()
      if (B.numRows () == 1 && B.numCols () == 1) {
        mv.update (alpha*B(0,0), A, beta);
        return;
      }

      // Create local map
      Teuchos::SerialComm<int> serialComm;
      map_type LocalMap (B.numRows (), A.getMap ()->getIndexBase (),
                         rcpFromRef<const Comm<int> > (serialComm),
                         Tpetra::LocallyReplicated, A.getMap ()->getNode ());
      // encapsulate Teuchos::SerialDenseMatrix data in ArrayView
      ArrayView<const Scalar> Bvalues (B.values (), B.stride () * B.numCols ());
      // create locally replicated MultiVector with a copy of this data
      MV B_mv (rcpFromRef (LocalMap), Bvalues, B.stride (), B.numCols ());
      mv.multiply (Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A, B_mv, beta);
    }

    /// \brief <tt>mv := alpha*A + beta*B</tt>
    ///
    /// The Tpetra specialization of this method ignores and
    /// completely overwrites any NaN or Inf entries in A.  Thus, it
    /// does <i>not</i> mean the same thing as <tt>mv := 0*mv +
    /// alpha*A + beta*B</tt> in IEEE 754 floating-point arithmetic.
    /// (Remember that NaN*0 = NaN.)
    static void
    MvAddMv (Scalar alpha,
             const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
             Scalar beta,
             const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
             Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      mv.update (alpha, A, beta, B, Teuchos::ScalarTraits<Scalar>::zero ());
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const Scalar& alpha)
    {
      mv.scale (alpha);
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const std::vector<BaseScalar>& alphas)
    {
      std::vector<Scalar> alphas_mp(alphas.size());
      const size_t sz = alphas.size();
      for (size_t i=0; i<sz; ++i)
        alphas_mp[i] = alphas[i];
      mv.scale (alphas_mp);
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const std::vector<Scalar>& alphas)
    {
      mv.scale (alphas);
    }

    static void
    MvTransMv (Scalar alpha,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
               Teuchos::SerialDenseMatrix<int,BaseScalar>& C)
    {
      Teuchos::SerialDenseMatrix<int,Scalar> C_mp(C.numRows(), C.numCols());
      MvTransMv(alpha, A, B, C_mp);
      for (int i=0; i<C.numRows(); i++) {
        for (int j=0; j<C.numCols(); j++) {
          BaseScalar v = 0.0;
          const s_ordinal sz = C_mp(i,j).size();
          for (s_ordinal k=0; k<sz; ++k)
            v += C_mp(i,j).fastAccessCoeff(k);
          C(i,j) = v;
        }
      }
    }

    static void
    MvTransMv (Scalar alpha,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
      using Tpetra::LocallyReplicated;
      using Teuchos::Comm;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcpFromRef;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      using Teuchos::SerialComm;
      typedef Tpetra::Map<LO,GO,Node> map_type;
      typedef Tpetra::MultiVector<Scalar,LO,GO,Node> mv_type;

#ifdef HAVE_BELOS_TPETRA_TIMERS
      const std::string timerName ("Belos::MVT::MvTransMv");
      RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = Teuchos::TimeMonitor::getNewCounter (timerName);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        timer.is_null (), std::logic_error, "Belos::MvTransMv: "
        "Failed to look up timer \"" << timerName << "\".  "
        "Please report this bug to the Belos developers.");

      // This starts the timer.  It will be stopped on scope exit.
      Teuchos::TimeMonitor timeMon (*timer);
#endif // HAVE_BELOS_TPETRA_TIMERS

      // Form alpha * A^H * B, then copy into the SerialDenseMatrix.
      // We will create a multivector C_mv from a a local map.  This
      // map has a serial comm, the purpose being to short-circuit the
      // MultiVector::reduce() call at the end of
      // MultiVector::multiply().  Otherwise, the reduced multivector
      // data would be copied back to the GPU, only to turn around and
      // have to get it back here.  This saves us a round trip for
      // this data.
      const int numRowsC = C.numRows ();
      const int numColsC = C.numCols ();
      const int strideC  = C.stride ();

      // Check if numRowsC == numColsC == 1, in which case we can call dot()
      if (numRowsC == 1 && numColsC == 1) {
        A.dot (B, Teuchos::ArrayView<Scalar> (C.values (), 1));
        return;
      }

      RCP<const Comm<int> > serialComm (new SerialComm<int> ());
      // create local map with serial comm
      RCP<const map_type> LocalMap =
        rcp (new map_type (numRowsC, 0, serialComm, LocallyReplicated,
                           A.getMap ()->getNode ()));
      // create local multivector to hold the result
      const bool INIT_TO_ZERO = true;
      mv_type C_mv (LocalMap, numColsC, INIT_TO_ZERO);

      // multiply result into local multivector
      C_mv.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, alpha, A, B,
                     Teuchos::ScalarTraits<Scalar>::zero ());
      // get comm
      RCP<const Comm<int> > pcomm = A.getMap ()->getComm ();
      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> C_view (C.values (), strideC*numColsC);
      if (pcomm->getSize () == 1) {
        // No accumulation to do; simply extract the multivector data
        // into C.  Extract a copy of the result into the array view
        // (and therefore, the SerialDenseMatrix).
        C_mv.get1dCopy (C_view, strideC);
      }
      else {
        // get a const host view of the data in C_mv
        Teuchos::ArrayRCP<const Scalar> C_mv_view = C_mv.get1dView ();
        if (strideC == numRowsC) {
          // sum all into C
          reduceAll<int,Scalar> (*pcomm, REDUCE_SUM, numColsC*numRowsC,
                                 C_mv_view.getRawPtr (), C_view.getRawPtr ());
        }
        else {
          // sum-all into temp, copy into C
          Teuchos::Array<Scalar> destBuff (numColsC * numRowsC);
          reduceAll<int,Scalar> (*pcomm, REDUCE_SUM, numColsC*numRowsC,
                                 C_mv_view.getRawPtr (), destBuff.getRawPtr ());
          for (int j = 0; j < numColsC; ++j) {
            for (int i = 0; i < numRowsC; ++i) {
              C_view[strideC*j+i] = destBuff[numRowsC*j+i];
            }
          }
        }
      }
    }

    //! For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    static void
    MvDot (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
           const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
           std::vector<BaseScalar>& dots)
    {
      const size_t numVecs = A.getNumVectors ();

      TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs != B.getNumVectors (), std::invalid_argument,
        "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): "
        "A and B must have the same number of columns.  "
        "A has " << numVecs << " column(s), "
        "but B has " << B.getNumVectors () << " column(s).");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        dots.size() < numVecs, std::invalid_argument,
        "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): "
        "The output array 'dots' must have room for all dot products.  "
        "A and B each have " << numVecs << " column(s), "
        "but 'dots' only has " << dots.size() << " entry(/ies).");
#endif // HAVE_TPETRA_DEBUG

      // Teuchos::ArrayView<Scalar> av (dots);
      // A.dot (B, av (0, numVecs));
      Teuchos::Array<Scalar> mp_dots(numVecs);
      A.dot(B, mp_dots);
      for (size_t i=0; i<numVecs; i++) {
        const s_ordinal sz = mp_dots[i].size();
        dots[i] = 0.0;
        for (s_ordinal j=0; j<sz; ++j)
          dots[i] += mp_dots[i].fastAccessCoeff(j);
      }
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
            std::vector<typename Teuchos::ScalarTraits<BaseScalar>::magnitudeType> &normvec,
            NormType type=TwoNorm)
    {

#ifdef HAVE_TPETRA_DEBUG
      typedef std::vector<int>::size_type size_type;
      TEUCHOS_TEST_FOR_EXCEPTION(
        normvec.size () < static_cast<size_type> (mv.getNumVectors ()),
        std::invalid_argument,
        "Belos::MultiVecTraits::MvNorm(mv,normvec): The normvec output "
        "argument must have at least as many entries as the number of vectors "
        "(columns) in the MultiVector mv.  normvec.size() = " << normvec.size ()
        << " < mv.getNumVectors() = " << mv.getNumVectors () << ".");
#endif
      const size_t num_mv = mv.getNumVectors();
      Teuchos::Array<Scalar> mp_norms(num_mv);
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(mp_norms);
      switch (type) {
      case OneNorm:
        mv.norm1(av(0,mv.getNumVectors()));

        for (size_t col=0; col<num_mv; ++col) {
          BaseScalar v = 0.0;
          const s_ordinal sz = mp_norms[col].size();
          for (s_ordinal i=0; i<sz; ++i) {
            const BaseScalar a = mp_norms[col].fastAccessCoeff(i);
            v += a;
          }
          normvec[col] = v;
        }
        break;
      case TwoNorm:
        mv.norm2(av(0,mv.getNumVectors()));

        for (size_t col=0; col<num_mv; ++col) {
          BaseScalar v = 0.0;
          const s_ordinal sz = mp_norms[col].size();
          for (s_ordinal i=0; i<sz; ++i) {
            const BaseScalar a = mp_norms[col].fastAccessCoeff(i);
            v += a*a;
          }
          normvec[col] = std::sqrt(v);
        }
        break;
      case InfNorm:
        mv.normInf(av(0,mv.getNumVectors()));

        for (size_t col=0; col<num_mv; ++col) {
          BaseScalar v = 0.0;
          const s_ordinal sz = mp_norms[col].size();
          for (s_ordinal i=0; i<sz; ++i) {
            const BaseScalar a = mp_norms[col].fastAccessCoeff(i);
            if (a > v) v = a;
          }
          normvec[col] = v;
        }
        break;
      default:
        // Throw logic_error rather than invalid_argument, because if
        // we get here, it's probably the fault of a Belos solver,
        // rather than a user giving Belos an invalid input.
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Belos::MultiVecTraits::MvNorm: Invalid NormType value " << type
          << ".  Valid values are OneNorm=" << OneNorm << ", TwoNorm="
          << TwoNorm <<", and InfNorm=" << InfNorm << ".  If you are a Belos "
          "user and have not modified Belos in any way, and you get this "
          "message, then this is probably a bug in the Belos solver you were "
          "using.  Please report this to the Belos developers.");
      }
    }

    static void
    SetBlock (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const std::vector<int>& index,
              Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      typedef std::vector<int>::size_type size_type;

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_type> (A.getNumVectors ()) < index.size (),
        std::invalid_argument,
        "Belos::MultiVecTraits::SetBlock(A,index,mv): The index argument must "
        "have the same number of entries as the number of columns in A.  "
        "index.size() = " << index.size () << " != A.getNumVectors() = "
        << A.getNumVectors () << ".");
#endif
      RCP<MV> mvsub = CloneViewNonConst (mv, index);
      if (static_cast<size_type> (A.getNumVectors ()) > index.size ()) {
        RCP<const MV> Asub = A.subView (Range1D (0, index.size () - 1));
        Tpetra::deep_copy (*mvsub, *Asub);
      }
      else {
        Tpetra::deep_copy (*mvsub, A);
      }
      mvsub = Teuchos::null;
    }

    static void
    SetBlock (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(maxInt < A.getNumVectors(), std::range_error,
                                   os.str() << "Number of columns in the input multi"
                                   "vector 'A' (a size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(maxInt < mv.getNumVectors(), std::range_error,
                                   os.str() << "Number of columns in the output multi"
                                   "vector 'mv' (a size_t) overflows int.");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors());
      const int numColsMv = static_cast<int> (mv.getNumVectors());
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Range lower bound must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numColsMv, std::invalid_argument,
          os.str() << "Range upper bound must be less than the number of "
          "columns " << numColsA << " in the 'mv' output argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() > numColsA, std::invalid_argument,
          os.str() << "Range must have no more elements than the number of "
          "columns " << numColsA << " in the 'A' input argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      Teuchos::RCP<MV> mv_view;
      if (index.lbound () == 0 && index.ubound () + 1 == numColsMv) {
        mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      } else {
        mv_view = CloneViewNonConst (mv, index);
      }

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      Teuchos::RCP<const MV> A_view;
      if (index.size () == numColsA) {
        A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      } else {
        A_view = CloneView (A, Teuchos::Range1D (0, index.size () - 1));
      }

      Tpetra::deep_copy (*mv_view, *A_view);
    }

    static void
    Assign (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      const char errPrefix[] = "Belos::MultiVecTraits::Assign(A, mv): ";

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the input multivector 'A' (a "
          "size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the output multivector 'mv' (a "
          "size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors ());
      const int numColsMv = static_cast<int> (mv.getNumVectors ());
      if (numColsA > numColsMv) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          numColsA > numColsMv, std::invalid_argument,
          errPrefix << "Input multivector 'A' has " << numColsA << " columns, "
          "but output multivector 'mv' has only " << numColsMv << " columns.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
      }
      // Assignment of Tpetra::MultiVector objects via operator=
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator= may
      // throw an std::runtime_error if the Maps are incompatible.
      if (numColsA == numColsMv) {
        Tpetra::deep_copy (mv, A);
      } else {
        Teuchos::RCP<MV> mv_view =
          CloneViewNonConst (mv, Teuchos::Range1D (0, numColsA-1));
        Tpetra::deep_copy (*mv_view, A);
      }
    }


    static void
    MvRandom (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      mv.randomize ();
    }

    static void
    MvInit (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
            Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero ())
    {
      mv.putScalar (alpha);
    }

    static void
    MvPrint (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os)
    {
      Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
      mv.describe (fos, Teuchos::VERB_EXTREME);
    }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Tpetra::MultiVector
    ///
    typedef Tpetra::TsqrAdaptor< Tpetra::MultiVector< Scalar, LO, GO, Node > > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)

  /// \brief Partial specialization of MultiVecTraits for MV = Tpetra::MultiVector.
  ///
  /// This interface lets Belos' solvers work directly with
  /// Tpetra::MultiVector objects as the multivector type
  /// (corresponding to the MV template parameter).
  ///
  /// The four template parameters of this partial specialization
  /// correspond exactly to the four template parameters of
  /// Tpetra::MultiVector.  See the Tpetra::MultiVector documentation
  /// for more information.
  template<class BaseScalar, class Storage, class LO, class GO, class Device>
  class MultiVecTraits<BaseScalar,
                       Tpetra::MultiVector< Sacado::MP::Vector<Storage>,
                                            LO,GO,
                                            Kokkos::Compat::KokkosDeviceWrapperNode<Device> > > {
  public:
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
  private:
    typedef typename Storage::ordinal_type s_ordinal;
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:
    typedef typename Tpetra::MultiVector<Scalar,LO,GO,Node>::dot_type dot_type;
    typedef typename Tpetra::MultiVector<Scalar,LO,GO,Node>::mag_type mag_type;

    /// \brief Create a new multivector with \c numvecs columns.
    ///
    /// The returned Tpetra::MultiVector has the same Tpetra::Map
    /// (distribution over one or more parallel processes) as \c mv.
    /// Its entries are not initialized and have undefined values.
    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    Clone (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const int numvecs)
    {
      return Teuchos::rcp (new MV (mv.getMap (), numvecs, false));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& X)
    {
      // Make a deep copy of X.  The one-argument copy constructor
      // does a shallow copy by default; the second argument tells it
      // to do a deep copy.
      Teuchos::RCP<MV> X_copy (new MV (X, Teuchos::Copy));
      // Make Tpetra::MultiVector use the new view semantics.  This is
      // a no-op for the Kokkos refactor version of Tpetra; it only
      // does something for the "classic" version of Tpetra.  This
      // shouldn't matter because Belos only handles MV through RCP
      // and through this interface anyway, but it doesn't hurt to set
      // it and make sure that it works.
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits::CloneCopy(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::runtime_error, fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors () << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subCopy, so we
      // don't have to check here.
      Teuchos::RCP<MV> X_copy = mv.subCopy (columns ());
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < GetNumberVecs(mv);
      if (! validRange) { // invalid range; generate error message
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneCopy(mv,index=["
           << index.lbound() << "," << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= GetNumberVecs(mv), std::invalid_argument,
          os.str() << "Index range exceeds number of vectors "
          << mv.getNumVectors() << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_copy = mv.subCopy (index);
      X_copy->setCopyOrView (Teuchos::View);
      return X_copy;
    }


    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits::CloneViewNonConst(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in
      // MultiVector::subViewNonConst, so we don't have to check here.
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (columns ());
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }


    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneViewNonConst(mv,index=["
           << index.lbound() << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative inde{x,ices}, which is "
          "not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<MV> X_view = mv.subViewNonConst (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }


    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar, "
        "Tpetra::MultiVector<...> >::CloneView(mv,index)";
      const size_t numVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= numVecs,
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << numVecs << " in the input MultiVector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Tpetra wants an array of size_t, not of int.
      Teuchos::Array<size_t> columns (index.size ());
      for (std::vector<int>::size_type j = 0; j < index.size (); ++j) {
        columns[j] = index[j];
      }
      // mfh 14 Aug 2014: Tpetra already detects and optimizes for a
      // continuous column index range in MultiVector::subView, so we
      // don't have to check here.
      Teuchos::RCP<const MV> X_view = mv.subView (columns);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::CloneView(mv, index=["
           << index.lbound () << ", " << index.ubound() << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numCols, std::invalid_argument,
          os.str() << "Index range exceeds number of vectors " << numCols
          << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      Teuchos::RCP<const MV> X_view = mv.subView (index);
      // FIXME (mfh 14 Aug 2014) For some reason I currently don't
      // understand, the Belos MVOpTester and/or OrthoManager tests
      // fail if I uncomment the line below.  This is true for both
      // the "classic" and Kokkos refactor versions of Tpetra.  I
      // don't know why this is the case.  Belos shouldn't care
      // whether Tpetra uses copy or view semantics, and the Kokkos
      // refactor version of Tpetra _always_ uses view semantics.
      // Nevertheless, the tests fail, so I will leave the following
      // line commented out for now.  This is true for both CloneView
      // overloads as well as both CloneViewNonConst overloads.

      //Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static int
    GetVecLength (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.getGlobalLength ();
    }

    static int
    GetNumberVecs (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.getNumVectors ();
    }

    static bool
    HasConstantStride (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv) {
      return mv.isConstantStride ();
    }

    static void
    MvTimesMatAddMv (const dot_type& alpha,
                     const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
                     const Teuchos::SerialDenseMatrix<int,dot_type>& B,
                     const dot_type& beta,
                     Tpetra::MultiVector<Scalar,LO,GO,Node>& C)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Kokkos::Compat::persistingView;

      // Check if numRowsB == numColsB == 1, in which case we can call update()
      const int numRowsB = B.numRows ();
      const int numColsB = B.numCols ();
      const int strideB  = B.stride ();
      if (numRowsB == 1 && numColsB == 1) {
        C.update (alpha*B(0,0), A, beta);
        return;
      }

      // Ensure A and C have constant stride
      RCP<const MV> Atmp, Ctmp;
      if (A.isConstantStride() == false) Atmp = rcp (new MV (A));
      else Atmp = rcp(&A,false);

      if (C.isConstantStride() == false) Ctmp = rcp (new MV (C));
      else Ctmp = rcp(&C,false);

      // Create flattened Kokkos::MultiVector's
      typedef Tpetra::MultiVector<dot_type,LO,GO,Node> FMV;
      typedef typename FMV::dual_view_type::t_dev view_type;
      typedef typename view_type::device_type device_type;
      typedef KokkosClassic::MultiVector<dot_type,Node> DKMV;
      view_type A_view = Atmp->template getLocalView<device_type>();
      view_type C_view = Ctmp->template getLocalView<device_type>();
      DKMV A_mv(A.getMap()->getNode());
      DKMV C_mv(A.getMap()->getNode());
      size_t A_stride[8], C_stride[8];
      A_view.stride(A_stride);
      C_view.stride(C_stride);
      const size_t LDA =
        A_view.dimension_1() > 1 ? A_stride[1] : A_view.dimension_0();
      const size_t LDC =
        C_view.dimension_1() > 1 ? C_stride[1] : C_view.dimension_0();
      A_mv.initializeValues(A_view.dimension_0(),
                            A_view.dimension_1(),
                            persistingView(A_view),
                            LDA);
      C_mv.initializeValues(C_view.dimension_0(),
                            C_view.dimension_1(),
                            persistingView(C_view),
                            LDC);

      // Create a view for B
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, device_type> b_view_type;
      typedef typename b_view_type::HostMirror b_host_view_type;
      b_host_view_type B_view(Kokkos::view_without_managing,
                              B.values(), strideB, numColsB);

      // Create view for B on the device -- need to be careful to get the
      // right stride to match B
      typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft, device_type> b_1d_view_type;
      b_1d_view_type B_1d_view_dev(Kokkos::allocate_without_initializing,
                                   "B", strideB*numColsB);
      b_view_type B_view_dev(Kokkos::view_without_managing,
                             B_1d_view_dev.ptr_on_device(), strideB, numColsB);
      Kokkos::deep_copy(B_view_dev, B_view);

      // Make a Kokkos::MultiVector for B
      DKMV B_mv(A.getMap()->getNode());
      B_mv.initializeValues(numRowsB, numColsB, persistingView(B_view_dev),
                            strideB);

      // Do local multiply
      typedef KokkosClassic::DefaultArithmetic<DKMV> MVT;
      MVT::GEMM(C_mv, Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha,
                A_mv, B_mv, beta);

      // Copy back into C if necessary
      if (C.isConstantStride() == false)
        Tpetra::deep_copy(C, *Ctmp);

      // Release A and C if necessary
      Atmp = Teuchos::null;
      Ctmp = Teuchos::null;
    }

    /// \brief <tt>mv := alpha*A + beta*B</tt>
    ///
    /// The Tpetra specialization of this method ignores and
    /// completely overwrites any NaN or Inf entries in A.  Thus, it
    /// does <i>not</i> mean the same thing as <tt>mv := 0*mv +
    /// alpha*A + beta*B</tt> in IEEE 754 floating-point arithmetic.
    /// (Remember that NaN*0 = NaN.)
    static void
    MvAddMv (Scalar alpha,
             const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
             Scalar beta,
             const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
             Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      mv.update (alpha, A, beta, B, Teuchos::ScalarTraits<Scalar>::zero ());
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const Scalar& alpha)
    {
      mv.scale (alpha);
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const std::vector<BaseScalar>& alphas)
    {
      std::vector<Scalar> alphas_mp(alphas.size());
      const size_t sz = alphas.size();
      for (size_t i=0; i<sz; ++i)
        alphas_mp[i] = alphas[i];
      mv.scale (alphas_mp);
    }

    static void
    MvScale (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
             const std::vector<Scalar>& alphas)
    {
      mv.scale (alphas);
    }

    static void
    MvTransMv (dot_type alpha,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
               Teuchos::SerialDenseMatrix<int,dot_type>& C)
    {
      using Teuchos::Comm;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      using Kokkos::Compat::persistingView;

      // Check if numRowsC == numColsC == 1, in which case we can call dot()
      const int numRowsC = C.numRows ();
      const int numColsC = C.numCols ();
      const int strideC  = C.stride ();
      if (numRowsC == 1 && numColsC == 1) {
        A.dot (B, Teuchos::ArrayView<dot_type> (C.values (), 1));
        return;
      }

      // Ensure A and B have constant stride
      RCP<const MV> Atmp, Btmp;
      if (A.isConstantStride() == false) Atmp = rcp (new MV (A));
      else Atmp = rcp(&A,false);

      if (B.isConstantStride() == false) Btmp = rcp (new MV (B));
      else Btmp = rcp(&B,false);

      // Create flattened Kokkos::MultiVector's
      typedef Tpetra::MultiVector<dot_type,LO,GO,Node> FMV;
      typedef typename FMV::dual_view_type::t_dev view_type;
      typedef typename view_type::device_type device_type;
      typedef KokkosClassic::MultiVector<dot_type,Node> DKMV;
      view_type A_view = Atmp->template getLocalView<device_type>();
      view_type B_view = Btmp->template getLocalView<device_type>();
      DKMV A_mv(A.getMap()->getNode());
      DKMV B_mv(A.getMap()->getNode());
      size_t A_stride[8], B_stride[8];
      A_view.stride(A_stride);
      B_view.stride(B_stride);
      const size_t LDA =
        A_view.dimension_1() > 1 ? A_stride[1] : A_view.dimension_0();
      const size_t LDB =
        B_view.dimension_1() > 1 ? B_stride[1] : B_view.dimension_0();
      A_mv.initializeValues(A_view.dimension_0(),
                            A_view.dimension_1(),
                            persistingView(A_view),
                            LDA);
      B_mv.initializeValues(B_view.dimension_0(),
                            B_view.dimension_1(),
                            persistingView(B_view),
                            LDB);

      // Create a view for C
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, device_type> c_view_type;
      typedef typename c_view_type::HostMirror c_host_view_type;
      c_host_view_type C_view(Kokkos::view_without_managing,
                              C.values(), strideC, numColsC);

      // Create view for C on the device -- need to be careful to get the
      // right stride to match C (allow setting to 0 for first-touch)
      typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft, device_type> c_1d_view_type;
      c_1d_view_type C_1d_view_dev("C", strideC*numColsC);
      c_view_type C_view_dev(Kokkos::view_without_managing,
                             C_1d_view_dev.ptr_on_device(), strideC, numColsC);

      // Make a Kokkos::MultiVector for C
      DKMV C_mv(A.getMap()->getNode());
      C_mv.initializeValues(numRowsC, numColsC, persistingView(C_view_dev),
                            strideC);

      // Do local multiply
      typedef KokkosClassic::DefaultArithmetic<DKMV> MVT;
      MVT::GEMM(C_mv, Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, alpha,
                A_mv, B_mv, Kokkos::Details::ArithTraits<dot_type>::zero());

      // Release A and B if necessary
      Atmp = Teuchos::null;
      Btmp = Teuchos::null;

      // reduce across processors -- could check for RDMA
      RCP<const Comm<int> > pcomm = A.getMap()->getComm ();
      if (pcomm->getSize () == 1)
        Kokkos::deep_copy(C_view, C_view_dev);
      else {
        typedef typename c_1d_view_type::HostMirror c_1d_host_view_type;
        c_1d_host_view_type C_1d_view_tmp(Kokkos::allocate_without_initializing,
                                          "C_tmp", strideC*numColsC);
        c_host_view_type C_view_tmp(Kokkos::view_without_managing,
                                    C_1d_view_tmp.ptr_on_device(),
                                    strideC, numColsC);
        Kokkos::deep_copy(C_view_tmp, C_view_dev);
        reduceAll<int> (*pcomm, REDUCE_SUM, strideC*numColsC,
                        C_view_tmp.ptr_on_device(), C_view.ptr_on_device());
      }
    }

    //! For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    static void
    MvDot (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
           const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
           std::vector<dot_type>& dots)
    {
      const size_t numVecs = A.getNumVectors ();

      TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs != B.getNumVectors (), std::invalid_argument,
        "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): "
        "A and B must have the same number of columns.  "
        "A has " << numVecs << " column(s), "
        "but B has " << B.getNumVectors () << " column(s).");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        dots.size() < numVecs, std::invalid_argument,
        "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvDot(A,B,dots): "
        "The output array 'dots' must have room for all dot products.  "
        "A and B each have " << numVecs << " column(s), "
        "but 'dots' only has " << dots.size() << " entry(/ies).");
#endif // HAVE_TPETRA_DEBUG

      Teuchos::ArrayView<dot_type> av (dots);
      A.dot (B, av (0, numVecs));
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
            std::vector<mag_type> &normvec,
            NormType type=TwoNorm)
    {

#ifdef HAVE_TPETRA_DEBUG
      typedef std::vector<int>::size_type size_type;
      TEUCHOS_TEST_FOR_EXCEPTION(
        normvec.size () < static_cast<size_type> (mv.getNumVectors ()),
        std::invalid_argument,
        "Belos::MultiVecTraits::MvNorm(mv,normvec): The normvec output "
        "argument must have at least as many entries as the number of vectors "
        "(columns) in the MultiVector mv.  normvec.size() = " << normvec.size ()
        << " < mv.getNumVectors() = " << mv.getNumVectors () << ".");
#endif

      Teuchos::ArrayView<mag_type> av(normvec);
      switch (type) {
      case OneNorm:
        mv.norm1(av(0,mv.getNumVectors()));
        break;
      case TwoNorm:
        mv.norm2(av(0,mv.getNumVectors()));
        break;
      case InfNorm:
        mv.normInf(av(0,mv.getNumVectors()));
        break;
      default:
        // Throw logic_error rather than invalid_argument, because if
        // we get here, it's probably the fault of a Belos solver,
        // rather than a user giving Belos an invalid input.
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Belos::MultiVecTraits::MvNorm: Invalid NormType value " << type
          << ".  Valid values are OneNorm=" << OneNorm << ", TwoNorm="
          << TwoNorm <<", and InfNorm=" << InfNorm << ".  If you are a Belos "
          "user and have not modified Belos in any way, and you get this "
          "message, then this is probably a bug in the Belos solver you were "
          "using.  Please report this to the Belos developers.");
      }
    }

    static void
    SetBlock (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const std::vector<int>& index,
              Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      typedef std::vector<int>::size_type size_type;

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_type> (A.getNumVectors ()) < index.size (),
        std::invalid_argument,
        "Belos::MultiVecTraits::SetBlock(A,index,mv): The index argument must "
        "have the same number of entries as the number of columns in A.  "
        "index.size() = " << index.size () << " != A.getNumVectors() = "
        << A.getNumVectors () << ".");
#endif
      RCP<MV> mvsub = CloneViewNonConst (mv, index);
      if (static_cast<size_type> (A.getNumVectors ()) > index.size ()) {
        RCP<const MV> Asub = A.subView (Range1D (0, index.size () - 1));
        Tpetra::deep_copy (*mvsub, *Asub);
      }
      else {
        Tpetra::deep_copy (*mvsub, A);
      }
      mvsub = Teuchos::null;
    }

    static void
    SetBlock (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(maxInt < A.getNumVectors(), std::range_error,
                                   os.str() << "Number of columns in the input multi"
                                   "vector 'A' (a size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(maxInt < mv.getNumVectors(), std::range_error,
                                   os.str() << "Number of columns in the output multi"
                                   "vector 'mv' (a size_t) overflows int.");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors());
      const int numColsMv = static_cast<int> (mv.getNumVectors());
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.lbound() < 0, std::invalid_argument,
          os.str() << "Range lower bound must be nonnegative.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= numColsMv, std::invalid_argument,
          os.str() << "Range upper bound must be less than the number of "
          "columns " << numColsA << " in the 'mv' output argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.size() > numColsA, std::invalid_argument,
          os.str() << "Range must have no more elements than the number of "
          "columns " << numColsA << " in the 'A' input argument.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      Teuchos::RCP<MV> mv_view;
      if (index.lbound () == 0 && index.ubound () + 1 == numColsMv) {
        mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      } else {
        mv_view = CloneViewNonConst (mv, index);
      }

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      Teuchos::RCP<const MV> A_view;
      if (index.size () == numColsA) {
        A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      } else {
        A_view = CloneView (A, Teuchos::Range1D (0, index.size () - 1));
      }

      Tpetra::deep_copy (*mv_view, *A_view);
    }

    static void
    Assign (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      const char errPrefix[] = "Belos::MultiVecTraits::Assign(A, mv): ";

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the input multivector 'A' (a "
          "size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the output multivector 'mv' (a "
          "size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Should never get here!");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors ());
      const int numColsMv = static_cast<int> (mv.getNumVectors ());
      if (numColsA > numColsMv) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          numColsA > numColsMv, std::invalid_argument,
          errPrefix << "Input multivector 'A' has " << numColsA << " columns, "
          "but output multivector 'mv' has only " << numColsMv << " columns.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
      }
      // Assignment of Tpetra::MultiVector objects via operator=
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator= may
      // throw an std::runtime_error if the Maps are incompatible.
      if (numColsA == numColsMv) {
        Tpetra::deep_copy (mv, A);
      } else {
        Teuchos::RCP<MV> mv_view =
          CloneViewNonConst (mv, Teuchos::Range1D (0, numColsA-1));
        Tpetra::deep_copy (*mv_view, A);
      }
    }


    static void
    MvRandom (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      mv.randomize ();
    }

    static void
    MvInit (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
            Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero ())
    {
      mv.putScalar (alpha);
    }

    static void
    MvPrint (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os)
    {
      Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
      mv.describe (fos, Teuchos::VERB_EXTREME);
    }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Tpetra::MultiVector
    ///
    typedef Tpetra::TsqrAdaptor< Tpetra::MultiVector< Scalar, LO, GO, Node > > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

#endif

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of OperatorTraits for Tpetra::Operator.
  template <class BaseScalar, class Storage, class LO, class GO, class Node>
  class OperatorTraits <BaseScalar,
                        Tpetra::MultiVector<Sacado::MP::Vector<Storage>,
                                             LO,GO,Node>,
                        Tpetra::Operator<Sacado::MP::Vector<Storage>,
                                         LO,GO,Node> >
  {
  public:
    typedef Sacado::MP::Vector<Storage> Scalar;
    static void
    Apply (const Tpetra::Operator<Scalar,LO,GO,Node>& Op,
           const Tpetra::MultiVector<Scalar,LO,GO,Node>& X,
           Tpetra::MultiVector<Scalar,LO,GO,Node>& Y,
           ETrans trans=NOTRANS)
    {
      switch (trans) {
        case NOTRANS:
          Op.apply(X,Y,Teuchos::NO_TRANS);
          break;
        case TRANS:
          Op.apply(X,Y,Teuchos::TRANS);
          break;
        case CONJTRANS:
          Op.apply(X,Y,Teuchos::CONJ_TRANS);
          break;
      default:
        const std::string scalarName = Teuchos::TypeNameTraits<Scalar>::name();
        const std::string loName = Teuchos::TypeNameTraits<LO>::name();
        const std::string goName = Teuchos::TypeNameTraits<GO>::name();
        const std::string nodeName = Teuchos::TypeNameTraits<Node>::name();
        const std::string otName = "Belos::OperatorTraits<" + scalarName
          + "," + loName + "," + goName + "," + nodeName + ">";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, otName << ": Should never "
                           "get here; fell through a switch statement.  "
                           "Please report this bug to the Belos developers.");
      }
    }

    static bool
    HasApplyTranspose (const Tpetra::Operator<Scalar,LO,GO,Node>& Op)
    {
      return Op.hasTransposeApply ();
    }
  };

  // Partial specialization for MV=Tpetra::MultiVector.
  template<class BaseScalar, class Storage, class LO, class GO, class Node>
  class MultiVecTraitsExt<BaseScalar,
                          Tpetra::MultiVector<Sacado::MP::Vector<Storage>,
                                              LO, GO, Node> > {
  public:
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    static ptrdiff_t GetGlobalLength( const MV& mv ) {
      return Teuchos::as<ptrdiff_t> (mv.getGlobalLength ());
    }
  };

} // end of Belos namespace

#endif
// end of file BELOS_TPETRA_ADAPTER_HPP
