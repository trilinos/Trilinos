//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOSMULTIVECTRAITS_TPETRA_HPP
#define BELOSMULTIVECTRAITS_TPETRA_HPP

/// \file BelosMultiVecTraits_Tpetra.hpp
/// \brief Partial specialization of Belos::MultiVecTraits for Tpetra objects.

#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_BELOS_TSQR
#  include "Tpetra_TsqrAdaptor.hpp"
#endif // HAVE_BELOS_TSQR

namespace Belos {

  /// \brief Specialization of MultiVecTraits for MV = Tpetra::MultiVector.
  ///
  /// This interface lets Belos' solvers work directly with
  /// Tpetra::MultiVector objects as the MultiVector type.  That type
  /// corresponds to the MV template parameter, which is the second
  /// template parameter (after Scalar) of most Belos classes.
  ///
  /// The four template parameters of this partial specialization
  /// correspond exactly to the four template parameters of
  /// Tpetra::MultiVector.  See the Tpetra::MultiVector documentation
  /// for more information.
  ///
  /// \tparam Scalar Same as template parameter 1 of Tpetra::MultiVector.
  /// \tparam LO Same as template parameter 2 of Tpetra::MultiVector.
  /// \tparam GO Same as template parameter 3 of Tpetra::MultiVector.
  /// \tparam Node Same as template parameter 4 of Tpetra::MultiVector.
  template<class Scalar, class LO, class GO, class Node>
  class MultiVecTraits<Scalar, ::Tpetra::MultiVector<Scalar,LO,GO,Node> > {
    typedef ::Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:
    /// \brief Create a new MultiVector with \c numVecs columns.
    ///
    /// The returned Tpetra::MultiVector has the same Tpetra::Map
    /// (distribution over one or more parallel processes) as \c X.
    /// Its entries are not initialized and have undefined values.
    static Teuchos::RCP<MV> Clone (const MV& X, const int numVecs) {
      Teuchos::RCP<MV> Y (new MV (X.getMap (), numVecs, false));
      Y->setCopyOrView (Teuchos::View);
      return Y;
    }

    //! Create and return a deep copy of X.
    static Teuchos::RCP<MV> CloneCopy (const MV& X)
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

    /// \brief Create and return a deep copy of the given columns of mv.
    ///
    /// \pre \code mv.getNumVectors() != 0 || index.size() == 0 \endcode
    /// \pre For all k such that <tt>0 <= k < index.size()</tt>,
    ///   \code
    ///   0 <= index[k] < mv.getNumVectors();
    ///   \endcode
    /// \post If this method returns Y:
    ///   \code
    ///   Y->isConstantStride() && Y->getNumVectors() == index.size();
    ///   \endcode
    static Teuchos::RCP<MV>
    CloneCopy (const MV& mv, const std::vector<int>& index)
    {
#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits::CloneCopy(mv,index)";
      const size_t inNumVecs = mv.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 && *std::min_element (index.begin (), index.end ()) < 0,
        std::runtime_error, fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        index.size () > 0 &&
        static_cast<size_t> (*std::max_element (index.begin (), index.end ())) >= inNumVecs,
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << inNumVecs << " of the input multivector mv.");
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

    /// \brief Create and return a deep copy of the given columns of mv.
    ///
    /// \post If this method returns Y:
    ///   \code
    ///   Y->isConstantStride() && Y->getNumVectors() == index.size();
    ///   \endcode
    static Teuchos::RCP<MV>
    CloneCopy (const MV& mv, const Teuchos::Range1D& index)
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

    static Teuchos::RCP<MV>
    CloneViewNonConst (MV& mv, const std::vector<int>& index)
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
      X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<MV>
    CloneViewNonConst (MV& mv, const Teuchos::Range1D& index)
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
      X_view->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const MV>
    CloneView (const MV& mv, const std::vector<int>& index)
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
      Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static Teuchos::RCP<const MV>
    CloneView (const MV& mv, const Teuchos::Range1D& index)
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
      Teuchos::rcp_const_cast<MV> (X_view)->setCopyOrView (Teuchos::View);
      return X_view;
    }

    static ptrdiff_t GetGlobalLength (const MV& mv) {
      return static_cast<ptrdiff_t> (mv.getGlobalLength ());
    }

    static int GetNumberVecs (const MV& mv) {
      return static_cast<int> (mv.getNumVectors ());
    }

    static bool HasConstantStride (const MV& mv) {
      return mv.isConstantStride ();
    }

    static void
    MvTimesMatAddMv (Scalar alpha,
                     const MV& A,
                     const Teuchos::SerialDenseMatrix<int, Scalar>& B,
                     Scalar beta,
                     MV& mv)
    {
      using Teuchos::ArrayView;
      using Teuchos::Comm;
      using Teuchos::rcpFromRef;
      typedef ::Tpetra::Map<LO, GO, Node> map_type;

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
                         ::Tpetra::LocallyReplicated);
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
             const MV& A,
             Scalar beta,
             const MV& B,
             MV& mv)
    {
      mv.update (alpha, A, beta, B, Teuchos::ScalarTraits<Scalar>::zero ());
    }

    static void MvScale (MV& mv, Scalar alpha) {
      mv.scale (alpha);
    }

    static void MvScale (MV& mv, const std::vector<Scalar>& alphas) {
      mv.scale (alphas);
    }

    static void
    MvTransMv (const Scalar alpha,
               const MV& A,
               const MV& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
      using ::Tpetra::LocallyReplicated;
      using Teuchos::Comm;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::TimeMonitor;
      typedef ::Tpetra::Map<LO,GO,Node> map_type;

#ifdef HAVE_BELOS_TPETRA_TIMERS
      const std::string timerName ("Belos::MVT::MvTransMv");
      RCP<Teuchos::Time> timer = TimeMonitor::lookupCounter (timerName);
      if (timer.is_null ()) {
        timer = TimeMonitor::getNewCounter (timerName);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        timer.is_null (), std::logic_error, "Belos::MvTransMv: "
        "Failed to look up timer \"" << timerName << "\".  "
        "Please report this bug to the Belos developers.");

      // This starts the timer.  It will be stopped on scope exit.
      TimeMonitor timeMon (*timer);
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
        if (alpha == Teuchos::ScalarTraits<Scalar>::zero ()) {
          // Short-circuit, as required by BLAS semantics.
          C(0,0) = alpha;
          return;
        }
        A.dot (B, Teuchos::ArrayView<Scalar> (C.values (), 1));
        if (alpha != Teuchos::ScalarTraits<Scalar>::one ()) {
          C(0,0) *= alpha;
        }
        return;
      }

      // get comm
      RCP<const Comm<int> > pcomm = A.getMap ()->getComm ();

      // create local map with comm
      RCP<const map_type> LocalMap =
        rcp (new map_type (numRowsC, 0, pcomm, LocallyReplicated));

      // create local multivector to hold the result
      const bool INIT_TO_ZERO = true;
      MV C_mv (LocalMap, numColsC, INIT_TO_ZERO);

      // multiply result into local multivector
      C_mv.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, alpha, A, B,
                     Teuchos::ScalarTraits<Scalar>::zero ());

      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> C_view (C.values (), strideC * numColsC);

      // No accumulation to do (since Tpetra has already done it);
      // simply extract the multivector data into C.
      // Extract a copy of the result into the array view
      // (and therefore, the SerialDenseMatrix).
      C_mv.get1dCopy (C_view, strideC);
    }

    //! For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    static void
    MvDot (const MV& A, const MV& B, std::vector<Scalar> &dots)
    {
      const size_t numVecs = A.getNumVectors ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs != B.getNumVectors (), std::invalid_argument,
        "Belos::MultiVecTraits::MvDot(A,B,dots): "
        "A and B must have the same number of columns.  "
        "A has " << numVecs << " column(s), "
        "but B has " << B.getNumVectors () << " column(s).");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        dots.size() < numVecs, std::invalid_argument,
        "Belos::MultiVecTraits::MvDot(A,B,dots): "
        "The output array 'dots' must have room for all dot products.  "
        "A and B each have " << numVecs << " column(s), "
        "but 'dots' only has " << dots.size() << " entry(/ies).");
#endif // HAVE_TPETRA_DEBUG

      Teuchos::ArrayView<Scalar> av (dots);
      A.dot (B, av (0, numVecs));
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const MV& mv,
            std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& normvec,
            NormType type=TwoNorm)
    {
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        normvec.size () < static_cast<std::vector<int>::size_type> (mv.getNumVectors ()),
        std::invalid_argument,
        "Belos::MultiVecTraits::MvNorm(mv,normvec): The normvec output "
        "argument must have at least as many entries as the number of vectors "
        "(columns) in the MultiVector mv.  normvec.size() = " << normvec.size ()
        << " < mv.getNumVectors() = " << mv.getNumVectors () << ".");
#endif // HAVE_TPETRA_DEBUG
      Teuchos::ArrayView<magnitude_type> av (normvec);
      switch (type) {
      case OneNorm:
        mv.norm1 (av (0, mv.getNumVectors ()));
        break;
      case TwoNorm:
        mv.norm2 (av (0, mv.getNumVectors ()));
        break;
      case InfNorm:
        mv.normInf (av (0,mv.getNumVectors ()));
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
    SetBlock (const MV& A, const std::vector<int>& index, MV& mv)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      const size_t inNumVecs = A.getNumVectors ();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        inNumVecs < static_cast<size_t> (index.size ()), std::invalid_argument,
        "Belos::MultiVecTraits::SetBlock(A,index,mv): 'index' argument must "
        "have no more entries as the number of columns in the input MultiVector"
        " A.  A.getNumVectors() = " << inNumVecs << " < index.size () = "
        << index.size () << ".");
#endif // HAVE_TPETRA_DEBUG
      RCP<MV> mvsub = CloneViewNonConst (mv, index);
      if (inNumVecs > static_cast<size_t> (index.size ())) {
        RCP<const MV> Asub = A.subView (Range1D (0, index.size () - 1));
        ::Tpetra::deep_copy (*mvsub, *Asub);
      } else {
        ::Tpetra::deep_copy (*mvsub, A);
      }
    }

    static void
    SetBlock (const MV& A, const Teuchos::Range1D& index, MV& mv)
    {
      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt =
        static_cast<size_t> (Teuchos::OrdinalTraits<int>::max ());
      const bool overflow =
        maxInt < A.getNumVectors () && maxInt < mv.getNumVectors ();
      if (overflow) {
        std::ostringstream os;
        os << "Belos::MultiVecTraits::SetBlock(A, index=[" << index.lbound ()
           << ", " << index.ubound () << "], mv): ";
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors (), std::range_error, os.str () << "Number "
          "of columns (size_t) in the input MultiVector 'A' overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors (), std::range_error, os.str () << "Number "
          "of columns (size_t) in the output MultiVector 'mv' overflows int.");
      }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors ());
      const int numColsMv = static_cast<int> (mv.getNumVectors ());
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex =
        index.lbound () >= 0 && index.ubound () < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size () <= numColsA;

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

      ::Tpetra::deep_copy (*mv_view, *A_view);
    }

    static void Assign (const MV& A, MV& mv)
    {
      const char errPrefix[] = "Belos::MultiVecTraits::Assign(A, mv): ";

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt =
        static_cast<size_t> (Teuchos::OrdinalTraits<int>::max ());
      const bool overflow =
        maxInt < A.getNumVectors () && maxInt < mv.getNumVectors ();
      if (overflow) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < A.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the input multivector 'A' "
          "(a size_t) overflows int.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          maxInt < mv.getNumVectors(), std::range_error,
          errPrefix << "Number of columns in the output multivector 'mv' "
          "(a size_t) overflows int.");
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
      if (numColsA == numColsMv) {
        ::Tpetra::deep_copy (mv, A);
      } else {
        Teuchos::RCP<MV> mv_view =
          CloneViewNonConst (mv, Teuchos::Range1D (0, numColsA - 1));
        ::Tpetra::deep_copy (*mv_view, A);
      }
    }

    static void MvRandom (MV& mv) {
      mv.randomize ();
    }

    static void
    MvInit (MV& mv, const Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero ())
    {
      mv.putScalar (alpha);
    }

    static void MvPrint (const MV& mv, std::ostream& os) {
      Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
      mv.describe (fos, Teuchos::VERB_EXTREME);
    }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Tpetra::MultiVector
    typedef ::Tpetra::TsqrAdaptor< ::Tpetra::MultiVector<Scalar, LO, GO, Node> > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

} // namespace Belos

#endif // BELOSMULTIVECTRAITS_TPETRA_HPP
