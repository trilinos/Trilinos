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

#ifndef BELOS_TPETRA_ADAPTER_HPP
#define BELOS_TPETRA_ADAPTER_HPP

#include <Kokkos_NodeTrace.hpp>

/*! \file BelosTpetraAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Tpetra concrete classes.
*/

// TODO: the assumption is made that the solver, multivector and operator are templated on the same scalar. this will need to be modified.

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>
#include <Kokkos_NodeAPIConfigDefs.hpp>

#ifdef HAVE_BELOS_TSQR
#  include <Tpetra_TsqrAdaptor.hpp>
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
  template<class Scalar, class LO, class GO, class Node>
  class MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node> > {
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:
#ifdef HAVE_BELOS_TPETRA_TIMERS
    static Teuchos::RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
#endif

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

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      KOKKOS_NODE_TRACE("Belos::MVT::CloneCopy(MV)")
      return Teuchos::rcp (new MV (mv));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > CloneCopy( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, const std::vector<int>& index )
    {
      using Teuchos::as;
      using Teuchos::Array;
      using Teuchos::Range1D;
      KOKKOS_NODE_TRACE("Belos::MVT::CloneCopy(MV,ind)")

#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneCopy(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0, std::runtime_error,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        as<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::runtime_error,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors() << " of the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Detect whether the index range is contiguous.
      // If it is, use the more efficient Range1D version of CloneCopy.
      for (std::vector<int>::size_type j = 1; j < index.size (); ++j) {
        if (index[j] != index[j-1] + 1) {
          // not contiguous; short circuit
          Array<size_t> stinds (index.begin (), index.end ());
          return mv.subCopy (stinds);
        }
      }
      // contiguous
      return mv.subCopy (Range1D (index.front (),index.back ()));
    }

    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneCopy (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const Teuchos::Range1D& index)
    {
      KOKKOS_NODE_TRACE("Belos::MVT::CloneCopy(MV,ind)")
      const bool validRange = index.size() > 0 &&
        index.lbound() >= 0 &&
        index.ubound() < GetNumberVecs(mv);
      if (! validRange) { // invalid range; generate error message
        std::ostringstream os;
        os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<...> >::"
          "CloneCopy(mv,index=[" << index.lbound() << ", " << index.ubound()
           << "]): ";
        TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
          os.str() << "Empty index range is not allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
          os.str() << "Index range includes negative index/ices, which is not "
          "allowed.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          index.ubound() >= GetNumberVecs(mv), std::invalid_argument,
          os.str() << "Index range exceeds number of vectors "
          << mv.getNumVectors() << " in the input multivector.");
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          os.str() << "Should never get here!");
      }
      return mv.subCopy (index);
    }


    static Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneViewNonConst (Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
                       const std::vector<int>& index)
    {
      using Teuchos::as;
      using Teuchos::Array;
      using Teuchos::Range1D;

#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneViewNonConst(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0, std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        as<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors () << " in the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Detect whether the index range is contiguous.
      // If it is, use the more efficient Range1D version of CloneViewNonConst.
      for (typename std::vector<int>::size_type j = 1; j<index.size (); ++j) {
        if (index[j] != index[j-1] + 1) {
          // not contiguous; short circuit
          Array<size_t> stinds (index.begin (), index.end ());
          return mv.subViewNonConst (stinds);
        }
      }
      // contiguous
      return mv.subViewNonConst (Range1D (index.front (), index.back ()));
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
      if (! validRange)
        {
          std::ostringstream os;
          os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<...> >::"
            "CloneViewNonConst(mv,index=[" << index.lbound() << ", "
             << index.ubound() << "]): ";
          TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
                             os.str() << "Empty index range is not allowed.");
          TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                             os.str() << "Index range includes negative "
                             "index/ices, which is not allowed.");
          TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numCols, std::invalid_argument,
                             os.str() << "Index range exceeds number of "
                             "vectors " << numCols << " in the input "
                             "multivector.");
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                             os.str() << "Should never get here!");
        }
      return mv.subViewNonConst (index);
    }


    static Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> >
    CloneView (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
               const std::vector<int>& index)
    {
      using Teuchos::as;
      using Teuchos::Array;
      using Teuchos::Range1D;

#ifdef HAVE_TPETRA_DEBUG
      const char fnName[] = "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::CloneView(mv,index)";
      TEUCHOS_TEST_FOR_EXCEPTION(
        *std::min_element (index.begin (), index.end ()) < 0,
        std::invalid_argument,
        fnName << ": All indices must be nonnegative.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        as<size_t> (*std::max_element (index.begin (), index.end ())) >= mv.getNumVectors (),
        std::invalid_argument,
        fnName << ": All indices must be strictly less than the number of "
        "columns " << mv.getNumVectors () << " in the input multivector mv.");
#endif // HAVE_TPETRA_DEBUG

      // Detect whether the index range is contiguous.
      // If it is, use the more efficient Range1D version of CloneView.
      for (typename std::vector<int>::size_type j = 1; j < index.size (); ++j) {
        if (index[j] != index[j-1] + 1) {
          // not contiguous; short circuit
          Array<size_t> stinds (index.begin (), index.end ());
          return mv.subView (stinds);
        }
      }
      // contiguous
      return mv.subView (Range1D (index.front (), index.back ()));
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
        os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<...> >::"
          "CloneView(mv, index=[" << index.lbound() << ", "
           << index.ubound() << "]): ";
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
      return mv.subView (index);
    }

    static int GetVecLength( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.getGlobalLength(); }

    static int GetNumberVecs( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.getNumVectors(); }

    static bool HasConstantStride( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    { return mv.isConstantStride(); }

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
      KOKKOS_NODE_TRACE("Belos::MVT::MvTimesMatAddMv()")
#ifdef HAVE_BELOS_TPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTimesMatAddMvTimer_);
#endif
      // create local map
      Teuchos::SerialComm<int> serialComm;
      // FIXME (mfh 07 Mar 2013) Shouldn't we use the same index base
      // for this Map as for the Maps of A and B?
      map_type LocalMap (B.numRows (), 0, rcpFromRef<const Comm<int> > (serialComm),
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
             const std::vector<Scalar>& alphas)
    {
      mv.scale (alphas);
    }

    static void
    MvTransMv (Scalar alpha,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
               const Tpetra::MultiVector<Scalar,LO,GO,Node>& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
      KOKKOS_NODE_TRACE("Belos::MVT::MvTransMv()")
#ifdef HAVE_BELOS_TPETRA_TIMERS
      Teuchos::TimeMonitor lcltimer(*mvTransMvTimer_);
#endif
      // form alpha * A^H * B, then copy into SDM
      // we will create a multivector C_mv from a a local map
      // this map has a serial comm, the purpose being to short-circuit the MultiVector::reduce() call at the end of MultiVector::multiply()
      // otherwise, the reduced multivector data would be copied back to the GPU, only to turn around and have to get it back here.
      // this saves us a round trip for this data.
      const int numRowsC = C.numRows(),
                numColsC = C.numCols(),
                strideC  = C.stride();
      Teuchos::SerialComm<int> scomm;
      // create local map with serial comm
      Tpetra::Map<LO,GO,Node> LocalMap(numRowsC, 0, Teuchos::rcpFromRef< const Teuchos::Comm<int> >(scomm), Tpetra::LocallyReplicated, A.getMap()->getNode());
      // create local multivector to hold the result
      const bool INIT_TO_ZERO = true;
      Tpetra::MultiVector<Scalar,LO,GO,Node> C_mv(Teuchos::rcpFromRef(LocalMap),numColsC, INIT_TO_ZERO);
      // multiply result into local multivector
      C_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,B,Teuchos::ScalarTraits<Scalar>::zero());
      // get comm
      Teuchos::RCP< const Teuchos::Comm<int> > pcomm = A.getMap()->getComm();
      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> C_view(C.values(),strideC*numColsC);
      if (pcomm->getSize() == 1) {
        // no accumulation to do; simply extract the multivector data into C
        // extract a copy of the result into the array view (and therefore, the SerialDenseMatrix)
        C_mv.get1dCopy(C_view,strideC);
      }
      else {
        // get a const host view of the data in C_mv
        Teuchos::ArrayRCP<const Scalar> C_mv_view = C_mv.get1dView();
        if (strideC == numRowsC) {
          // sumall into C
          Teuchos::reduceAll<int,Scalar>(*pcomm,Teuchos::REDUCE_SUM,numColsC*numRowsC,C_mv_view.getRawPtr(),C_view.getRawPtr());
        }
        else {
          // sumall into temp, copy into C
          Teuchos::Array<Scalar> destBuff(numColsC*numRowsC);
          Teuchos::reduceAll<int,Scalar>(*pcomm,Teuchos::REDUCE_SUM,numColsC*numRowsC,C_mv_view.getRawPtr(),destBuff.getRawPtr());
          for (int j=0; j < numColsC; ++j) {
            for (int i=0; i < numRowsC; ++i) {
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
           std::vector<Scalar>& dots)
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

      Teuchos::ArrayView<Scalar> av (dots);
      A.dot (B, av (0, numVecs));
    }

    //! For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv,
            std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
            NormType type=TwoNorm)
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.getNumVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
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
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::MvNorm: "
          "Invalid NormType value " << type << ".  Valid values are OneNorm="
          << OneNorm << ", TwoNorm=" << TwoNorm <<", and InfNorm=" << InfNorm
          << ".  If you are a Belos user and have not modified Belos in any "
          "way, and you get this message, then this is probably a bug in the "
          "Belos solver you were using.  Please report this to the Belos "
          "developers.");
      }
    }

    static void SetBlock( const Tpetra::MultiVector<Scalar,LO,GO,Node>& A, const std::vector<int>& index, Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      KOKKOS_NODE_TRACE("Belos::MVT::SetBlock()")
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.getNumVectors() < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,Tpetra::MultiVector>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > mvsub = CloneViewNonConst(mv,index);
      if ((typename std::vector<int>::size_type)A.getNumVectors() > index.size()) {
        Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > Asub = A.subView(Teuchos::Range1D(0,index.size()-1));
        (*mvsub) = (*Asub);
      }
      else {
        (*mvsub) = A;
      }
      mvsub = Teuchos::null;
    }

    static void
    SetBlock (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
              const Teuchos::Range1D& index,
              Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      KOKKOS_NODE_TRACE("Belos::MVT::SetBlock()")

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow)
        {
          std::ostringstream os;
          os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar, ..."
            "> >::SetBlock(A, index=[" << index.lbound() << ", "
             << index.ubound() << "], mv): ";
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

      if (! validIndex || ! validSource)
        {
          std::ostringstream os;
          os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar, ..."
            "> >::SetBlock(A, index=[" << index.lbound() << ", "
             << index.ubound() << "], mv): ";
          TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
                             os.str() << "Range lower bound must be nonnegative.");
          TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numColsMv, std::invalid_argument,
                             os.str() << "Range upper bound must be less than "
                             "the number of columns " << numColsA << " in the "
                             "'mv' output argument.");
          TEUCHOS_TEST_FOR_EXCEPTION(index.size() > numColsA, std::invalid_argument,
                             os.str() << "Range must have no more elements than"
                             " the number of columns " << numColsA << " in the "
                             "'A' input argument.");
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
        }
      typedef Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > MV_ptr;
      typedef Teuchos::RCP<const Tpetra::MultiVector<Scalar,LO,GO,Node> > const_MV_ptr;

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      MV_ptr mv_view;
      if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
        mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else
        mv_view = CloneViewNonConst (mv, index);

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      const_MV_ptr A_view;
      if (index.size() == numColsA)
        A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
        A_view = CloneView (A, Teuchos::Range1D(0, index.size()-1));

      // Assignment of Tpetra::MultiVector objects via operator=()
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator=()
      // will throw an std::runtime_error if the Maps are
      // incompatible.
      *mv_view = *A_view;
    }

    static void
    Assign (const Tpetra::MultiVector<Scalar,LO,GO,Node>& A,
            Tpetra::MultiVector<Scalar,LO,GO,Node>& mv)
    {
      KOKKOS_NODE_TRACE("Belos::MVT::Assign()")

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of Tpetra::MultiVector is a deep copy.

      // Tpetra::MultiVector::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow)
        {
          std::ostringstream os;
          os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar, ..."
            "> >::Assign(A, mv): ";
          TEUCHOS_TEST_FOR_EXCEPTION(maxInt < A.getNumVectors(), std::range_error,
                             os.str() << "Number of columns in the input multi"
                             "vector 'A' (a size_t) overflows int.");
          TEUCHOS_TEST_FOR_EXCEPTION(maxInt < mv.getNumVectors(), std::range_error,
                             os.str() << "Number of columns in the output multi"
                             "vector 'mv' (a size_t) overflows int.");
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
        }
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors());
      const int numColsMv = static_cast<int> (mv.getNumVectors());
      if (numColsA > numColsMv)
        {
          std::ostringstream os;
          os << "Belos::MultiVecTraits<Scalar, Tpetra::MultiVector<Scalar, ..."
            "> >::Assign(A, mv): ";
          TEUCHOS_TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
                             os.str() << "Input multivector 'A' has "
                             << numColsA << " columns, but output multivector "
                             "'mv' has only " << numColsMv << " columns.");
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
        }
      // Assignment of Tpetra::MultiVector objects via operator=()
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator=()
      // will throw an std::runtime_error if the Maps are
      // incompatible.
      if (numColsA == numColsMv)
        mv = A;
      else
        {
          Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node> > mv_view =
            CloneViewNonConst (mv, Teuchos::Range1D(0, numColsA-1));
          *mv_view = A;
        }
    }


    static void MvRandom( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv )
    {
      KOKKOS_NODE_TRACE("Belos::MVT::randomize()")
      mv.randomize();
    }

    static void MvInit( Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    static void MvPrint( const Tpetra::MultiVector<Scalar,LO,GO,Node>& mv, std::ostream& os )
    {
      Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
      mv.describe(fos,Teuchos::VERB_EXTREME);
    }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Tpetra::MultiVector
    ///
    typedef Tpetra::TsqrAdaptor< Tpetra::MultiVector< Scalar, LO, GO, Node > > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of OperatorTraits for Tpetra::Operator.
  template <class Scalar, class LO, class GO, class Node>
  class OperatorTraits <Scalar, Tpetra::MultiVector<Scalar,LO,GO,Node>, Tpetra::Operator<Scalar,LO,GO,Node> >
  {
  public:
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

} // end of Belos namespace

#endif
// end of file BELOS_TPETRA_ADAPTER_HPP
