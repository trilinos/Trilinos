// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_ADAPTER_MP_VECTOR_HPP
#define BELOS_TPETRA_ADAPTER_MP_VECTOR_HPP

#include "BelosTpetraAdapter.hpp"
#include "Belos_Tpetra_MP_Vector.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "KokkosBlas.hpp"

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
  template<class Storage, class LO, class GO, class Node>
  class MultiVecTraits<typename Storage::value_type,
                       Tpetra::MultiVector< Sacado::MP::Vector<Storage>,
                                            LO, GO, Node > > {
  public:
    typedef typename Storage::ordinal_type s_ordinal;
    typedef typename Storage::value_type BaseScalar;
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
  public:
    typedef typename Tpetra::MultiVector<Scalar,LO,GO,Node>::dot_type dot_type;
    typedef typename Tpetra::MultiVector<Scalar,LO,GO,Node>::mag_type mag_type;

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

    template <class DstType, class SrcType>
    static void deep_copy_2d_view_with_intercessory_space(
      const DstType& dst,
      const SrcType& src
    )
    {
      if (std::is_same<typename SrcType::memory_space, Kokkos::HostSpace>::value)
      {
        auto dst_i = Kokkos::create_mirror(dst);
        for (size_t i = 0; i < src.extent(0); i++)
        {
          for (size_t j = 0; j < src.extent(1); j++)
          {
            dst_i(i, j) = src(i, j);
          }
        }
        Kokkos::deep_copy(dst, dst_i);
      }
      else
      {
        auto src_i = Kokkos::create_mirror(src);
        Kokkos::deep_copy(src_i, src);
        for (size_t i = 0; i < src.extent(0); i++)
        {
          for (size_t j = 0; j < src.extent(1); j++)
          {
            dst(i, j) = src_i(i, j);
          }
        }
      }
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

      // Check if numRowsB == numColsB == 1, in which case we can call update()
      const int numRowsB = B.numRows ();
      const int numColsB = B.numCols ();
      const int strideB  = B.stride ();
      if (numRowsB == 1 && numColsB == 1) {
        C.update (alpha*B(0,0), A, beta);
        return;
      }

      // Ensure A and C have constant stride
      RCP<const MV> Atmp;
      RCP<      MV> Ctmp;
      if (A.isConstantStride() == false) Atmp = rcp (new MV (A, Teuchos::Copy));
      else Atmp = rcp(&A,false);

      if (C.isConstantStride() == false) Ctmp = rcp (new MV (C, Teuchos::Copy));
      else Ctmp = rcp(&C,false);

      // Create flattened view's
      typedef Tpetra::MultiVector<dot_type,LO,GO,Node> FMV;
      typedef typename FMV::dual_view_type::t_dev flat_view_type;
      typedef typename flat_view_type::execution_space execution_space;
      typename flat_view_type::const_type flat_A_view = Atmp->getLocalViewDevice(Tpetra::Access::ReadOnly);
      flat_view_type flat_C_view = Ctmp->getLocalViewDevice(Tpetra::Access::OverwriteAll);

      // Create a view for B on the host
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, Kokkos::HostSpace> b_host_view_type;
      b_host_view_type B_view_host_input( B.values(), strideB, numColsB);
      auto B_view_host = Kokkos::subview(B_view_host_input,
                                         Kokkos::pair<int,int>(0,numRowsB),
                                         Kokkos::pair<int,int>(0,numColsB));

      // Create view for B on the device -- need to be careful to get the
      // right stride to match B
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, execution_space> b_view_type;
      typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft, execution_space> b_1d_view_type;
      b_1d_view_type B_1d_view_dev(Kokkos::ViewAllocateWithoutInitializing("B"), numRowsB*numColsB);
      b_view_type B_view_dev( B_1d_view_dev.data(), numRowsB, numColsB);

      if (Kokkos::SpaceAccessibility<Kokkos::HostSpace, typename execution_space::memory_space>::accessible)
      {
        Kokkos::deep_copy(B_view_dev, B_view_host);
      }
      else
      {
        deep_copy_2d_view_with_intercessory_space(B_view_dev, B_view_host);
      }

      // Do local multiply
      {
        const char ctransA = 'N', ctransB = 'N';
        KokkosBlas::gemm (
          &ctransA, &ctransB,
          alpha, flat_A_view, B_view_dev, beta, flat_C_view);
      }
      // Copy back to C if we made a copy
      if (C.isConstantStride() == false)
        C.assign(*Ctmp);
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

      // Check if numRowsC == numColsC == 1, in which case we can call dot()
      const int numRowsC = C.numRows ();
      const int numColsC = C.numCols ();
      const int strideC  = C.stride ();
      if (numRowsC == 1 && numColsC == 1) {
        if (alpha == Teuchos::ScalarTraits<Scalar>::zero ()) {
          // Short-circuit, as required by BLAS semantics.
          C(0,0) = alpha;
          return;
        }
        A.dot (B, Teuchos::ArrayView<dot_type> (C.values (), 1));
        if (alpha != Teuchos::ScalarTraits<Scalar>::one ()) {
          C(0,0) *= alpha;
        }
        return;
      }

      // Ensure A and B have constant stride
      RCP<const MV> Atmp, Btmp;
      if (A.isConstantStride() == false) Atmp = rcp (new MV (A, Teuchos::Copy));
      else Atmp = rcp(&A,false);

      if (B.isConstantStride() == false) Btmp = rcp (new MV (B, Teuchos::Copy));
      else Btmp = rcp(&B,false);

      // Create flattened Kokkos::MultiVector's
      typedef Tpetra::MultiVector<dot_type,LO,GO,Node> FMV;
      typedef typename FMV::dual_view_type::t_dev flat_view_type;
      typedef typename flat_view_type::execution_space execution_space;
      typename flat_view_type::const_type flat_A_view = Atmp->getLocalViewDevice(Tpetra::Access::ReadOnly);
      typename flat_view_type::const_type flat_B_view = Btmp->getLocalViewDevice(Tpetra::Access::ReadOnly);

      // Create a view for C on the host
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, Kokkos::HostSpace> c_host_view_type;
      c_host_view_type C_view_host_input( C.values(), strideC, numColsC);
      auto C_view_host = Kokkos::subview(C_view_host_input, 
                                         Kokkos::pair<int,int>(0,numRowsC), 
                                         Kokkos::pair<int,int>(0,numColsC));

      // Create view for C on the device -- need to be careful to get the
      // right stride to match C (allow setting to 0 for first-touch)
      typedef Kokkos::View<dot_type**, Kokkos::LayoutLeft, execution_space> c_view_type;
      typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft, execution_space> c_1d_view_type;
      c_1d_view_type C_1d_view_dev("C", numRowsC*numColsC);
      c_view_type C_view_dev( C_1d_view_dev.data(), numRowsC, numColsC);

      // Do local multiply
      {
        const char ctransA = 'C', ctransB = 'N';
        KokkosBlas::gemm (
          &ctransA, &ctransB,
          alpha, flat_A_view, flat_B_view,
          Kokkos::ArithTraits<dot_type>::zero(),
          C_view_dev);
      }

      // reduce across processors -- could check for RDMA
      RCP<const Comm<int> > pcomm = A.getMap()->getComm ();
      if (pcomm->getSize () == 1)
      {
        if (Kokkos::SpaceAccessibility<execution_space, Kokkos::HostSpace>::accessible)
        {
          Kokkos::deep_copy(C_view_host, C_view_dev);
        }
        else
        {
          deep_copy_2d_view_with_intercessory_space(C_view_host, C_view_dev);
        }
      }
      else
      {
        typedef Kokkos::View<dot_type*, Kokkos::LayoutLeft, Kokkos::HostSpace> c_1d_host_view_type;
        c_1d_host_view_type C_1d_view_tmp(Kokkos::ViewAllocateWithoutInitializing("C_tmp"), strideC*numColsC);
        c_host_view_type C_view_tmp_input( C_1d_view_tmp.data(), strideC, numColsC);
        auto C_view_tmp = Kokkos::subview(C_view_tmp_input,
                                          Kokkos::pair<int,int>(0,numRowsC),
                                          Kokkos::pair<int,int>(0,numColsC));
        
        if (Kokkos::SpaceAccessibility<execution_space, Kokkos::HostSpace>::accessible)
        {
          Kokkos::deep_copy(C_view_tmp, C_view_dev);
        }
        else
        {
          deep_copy_2d_view_with_intercessory_space(C_view_tmp, C_view_dev);
        }

        reduceAll<int> (*pcomm, REDUCE_SUM, strideC*numColsC,
                        C_view_tmp.data(),
                        C_view_host.data());
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
  template <class Storage, class LO, class GO, class Node>
  class OperatorTraits <typename Storage::value_type,
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

} // end of Belos namespace

#endif
// end of file BELOS_TPETRA_ADAPTER_HPP
