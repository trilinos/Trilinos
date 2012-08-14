//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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

#ifndef KOKKOS_DEFAULTSPARSEOPS_HPP
#define KOKKOS_DEFAULTSPARSEOPS_HPP

#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Describable.hpp>
#include <iterator>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseSolveKernelOps.hpp"
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

namespace Kokkos {

  namespace details {

    template <class O>
    struct rowPtrsInitKernel {
      O *rowPtrs;
      inline KERNEL_PREFIX void execute(int i) const {
        rowPtrs[i] = 0;
      }
    };

    template <class O, class T>
    struct valsInitKernel {
      const O *rowPtrs;
      T * entries;
      inline KERNEL_PREFIX void execute(int i) const {
        T *beg = entries+rowPtrs[i],
          *end = entries+rowPtrs[i+1];
        while (beg != end) {
          *beg++ = Teuchos::ScalarTraits<T>::zero();
        }
      }
    };

    class FirstTouchCRSAllocator {
      public:
      //! \brief Allocate and initialize the storage for the matrix values.
      template <class Ordinal, class Node>
      static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &numEntriesPerRow)
      {
        const Ordinal numrows = numEntriesPerRow.size();
        details::rowPtrsInitKernel<Ordinal> kern;
        // allocate
        kern.rowPtrs = new Ordinal[numrows+1];
        // parallel first touch
        node->parallel_for(0,numrows+1,kern);
        // encapsulate
        ArrayRCP<Ordinal> ptrs = arcp<Ordinal>(kern.rowPtrs,0,numrows+1,true);
        // compute in serial. parallelize later, perhaps; it's only O(N)
        ptrs[0] = 0;
        std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
        return ptrs;
      }

      //! \brief Allocate and initialize the storage for a sparse graph.
      template <class T, class Ordinal, class Node>
      static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const Ordinal> &rowPtrs)
      {
        const Ordinal totalNumEntries = *(rowPtrs.end()-1);
        const Ordinal numRows = rowPtrs.size() - 1;
        details::valsInitKernel<Ordinal,T> kern;
        ArrayRCP<T> vals;
        if (totalNumEntries > 0) {
          // allocate
          kern.entries = new T[totalNumEntries];
          kern.rowPtrs = rowPtrs.getRawPtr();
          // first touch
          node->parallel_for(0,numRows,kern);
          // encapsulate
          vals = arcp<T>(kern.entries,0,totalNumEntries,true);
        }
        return vals;
      }
    };

    class DefaultCRSAllocator {
      public:
      //! \brief Allocate and initialize the storage for the matrix values.
      template <class Ordinal, class Node>
      static ArrayRCP<Ordinal> allocRowPtrs(const RCP<Node> &node, const ArrayView<const Ordinal> &numEntriesPerRow)
      {
        ArrayRCP<Ordinal> ptrs = arcp<Ordinal>( numEntriesPerRow.size() + 1 );
        ptrs[0] = 0;
        std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
        return ptrs;
      }

      //! \brief Allocate and initialize the storage for a sparse graph.
      template <class T, class Ordinal, class Node>
      static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
      {
        const Ordinal totalNumEntries = *(rowPtrs.end()-1);
        // alloc data
        ArrayRCP<T> vals;
        if (totalNumEntries > 0) vals = arcp<T>(totalNumEntries);
        std::fill( vals.begin(), vals.end(), Teuchos::ScalarTraits<T>::zero() );
        return vals;
      }
    };

  }

  //! \class DefaultCrsGraph
  /** \brief Default implementation of CRS sparse graph, using generic kernels and suitable for host-based nodes.
  */
  template <class Ordinal,
            class Node>
  class DefaultCrsGraph : public CrsGraphBase<Ordinal,Node>
  {
    public:
      DefaultCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params);
      bool isEmpty() const;
      void setStructure(const ArrayRCP<const size_t> &ptrs,
                        const ArrayRCP<const Ordinal> &inds);
      inline ArrayRCP<const size_t>  getPointers() const;
      inline void setSmallPointers(const ArrayRCP<const Ordinal> &ptrs);
      inline ArrayRCP<const Ordinal> getSmallPointers() const;
      inline ArrayRCP<const Ordinal> getIndices() const;
      inline bool isInitialized() const;
      inline void setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag);
      inline void getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;
    private:
      ArrayRCP<const  size_t> big_ptrs_;
      ArrayRCP<const Ordinal> sml_ptrs_;
      ArrayRCP<const Ordinal> inds_;
      bool isInitialized_;
      bool isEmpty_;
      Teuchos::EUplo  tri_uplo_;
      Teuchos::EDiag unit_diag_;
  };

  //! \class DefaultCrsMatrix
  /** \brief Default implementation of CRS sparse matrix, using generic kernels and suitable for host-based nodes.
  */
  template <class Scalar,
            class Ordinal,
            class Node>
  class DefaultCrsMatrix : public CrsMatrixBase<Scalar,Ordinal,Node>
  {
    public:
      DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params);
      void setValues(const ArrayRCP<const Scalar> &vals);
      inline ArrayRCP<const Scalar> getValues() const;
      inline bool isInitialized() const;
    private:
      ArrayRCP<const Scalar> vals_;
      bool isInitialized_;
  };

  template <class Ordinal, class Node>
  DefaultCrsGraph<Ordinal,Node>::DefaultCrsGraph(Ordinal numRows, Ordinal numCols, const RCP<Node> &node, const RCP<ParameterList> &params)
  : CrsGraphBase<Ordinal,Node>(numRows,numCols,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isEmpty() const
  { return isEmpty_; }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setStructure(
                      const ArrayRCP<const size_t>  &ptrs,
                      const ArrayRCP<const Ordinal> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const Ordinal numrows = this->getNumRows();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptrs.is_null (),
      std::runtime_error,
      ": The input array 'ptrs' must be nonnull, even for a matrix with zero "
      "rows.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptrs.size() != (size_t) numrows+1,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptrs.size() = " << ptrs.size() << " != numrows+1 = "
      << (numrows+1) << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptrs[0] != 0,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptrs[0] = " << ptrs[0] << " != 0.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) inds.size() != (size_t) ptrs[numrows],
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- inds.size() = " << inds.size() << " != ptrs[numrows="
      << numrows << "] = " << ptrs[numrows] << ".");

    const Ordinal numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_,
      std::runtime_error,
      " matrix has already been initialized."
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    big_ptrs_ = ptrs;
    inds_ = inds;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  ArrayRCP<const size_t> DefaultCrsGraph<Ordinal,Node>::getPointers() const
  { return big_ptrs_; }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsGraph<Ordinal,Node>::getSmallPointers() const
  { return sml_ptrs_; }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setSmallPointers(const ArrayRCP<const Ordinal> &ptrs)
  {
    sml_ptrs_ = ptrs;
    big_ptrs_ = null;
  }

  template <class Ordinal, class Node>
  ArrayRCP<const Ordinal> DefaultCrsGraph<Ordinal,Node>::getIndices() const
  { return inds_; }

  template <class Ordinal, class Node>
  bool DefaultCrsGraph<Ordinal,Node>::isInitialized() const
  { return isInitialized_; }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::setMatDesc(Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    tri_uplo_ = uplo;
    unit_diag_ = diag;
  }

  template <class Ordinal, class Node>
  void DefaultCrsGraph<Ordinal,Node>::getMatDesc(Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = tri_uplo_;
    diag = unit_diag_;
  }

  template <class Scalar, class Ordinal, class Node>
  DefaultCrsMatrix<Scalar,Ordinal,Node>::DefaultCrsMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph, const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,Ordinal,Node>(graph,params)
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultCrsMatrix<Scalar,Ordinal,Node>::setValues(const ArrayRCP<const Scalar> &vals)
  {
    std::string tfecfFuncName("setValues(vals)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix is already initialized."
    )
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  ArrayRCP<const Scalar> DefaultCrsMatrix<Scalar,Ordinal,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Ordinal, class Node>
  bool DefaultCrsMatrix<Scalar,Ordinal,Node>::isInitialized() const
  {
    return isInitialized_;
  }

  /// \class DefaultHostSparseOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Ordinal, class Node, class Allocator = details::DefaultCRSAllocator>
  class DefaultHostSparseOps : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Local sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for local
    /// sparse operations for a scalar type S2, which may be different
    /// from \c Scalar.
    ///
    /// This class' typedef is used by Tpetra::CrsMatrix to bind a
    /// potentially "void" scalar type to the appropriate scalar.
    /// Other local sparse ops implementations (especially those that
    /// wrap third-party libraries implementing sparse kernels) might
    /// use this to provide a "fall-back" sparse ops implementation of
    /// a possibly different type, if the third-party library does not
    /// support scalar type S2.
    ///
    /// In the case of DefaultHostSparseOps, this class' typedef
    /// always specifies a specialization of \c DefaultHostSparseOps,
    /// regardless of the scalar type S2.  This is not necessarily
    /// true of other implementations of local sparse ops.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef DefaultHostSparseOps<S2,Ordinal,Node,Allocator> other_type;
    };

    /// \brief Local sparse operations type for a different ordinal type.
    ///
    /// The bind_scalar struct defines the type responsible for local
    /// sparse operations for an ordinal type O2, which may be
    /// different from \c Ordinal.
    ///
    /// In the case of DefaultHostSparseOps, this class' typedef
    /// always specifies a specialization of \c DefaultHostSparseOps,
    /// regardless of the ordinal type S2.  This is not necessarily
    /// true of other implementations of local sparse ops.  Other
    /// local sparse ops implementations (especially those that wrap
    /// third-party libraries implementing sparse kernels) might use
    /// this to provide a "fall-back" sparse ops implementation of a
    /// possibly different type, if the third-party library does not
    /// support ordinal type O2.
    ///
    /// \tparam O2 An ordinal type possibly different from \c Ordinal.
    template <class O2>
    struct bind_ordinal {
      typedef DefaultHostSparseOps<Scalar,O2,Node,Allocator> other_type;
    };

    //@}
    //! \name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    DefaultHostSparseOps(const RCP<Node> &node);

    //! Constructor accepting and retaining a node object, and taking parameters.
    DefaultHostSparseOps(const RCP<Node> &node, Teuchos::ParameterList& params);

    //! Destructor
    ~DefaultHostSparseOps();

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os <<  "Kokkos::DefaultHostSparseOps<"
         << "Scalar=" << TypeNameTraits<Scalar>::name()
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::EVerbosityLevel;
      using Teuchos::includesVerbLevel;
      using Teuchos::OSTab;
      using Teuchos::rcpFromRef;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      using std::endl;

      // Interpret the default verbosity level as VERB_MEDIUM.
      const EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
        out << this->description();

        if (includesVerbLevel (vl, VERB_MEDIUM)) { // vl >= VERB_MEDIUM
          out << ":" << endl;
          OSTab tab1 (rcpFromRef (out));

          out << "isInitialized_ = " << isInitialized_ << endl;
          if (isInitialized_) {
            std::string triUplo ("INVALID");
            if (tri_uplo_ == Teuchos::UNDEF_TRI) {
              triUplo = "UNDEF_TRI";
            }
            else if (tri_uplo_ == Teuchos::LOWER_TRI) {
              triUplo = "LOWER_TRI";
            }
            else if (tri_uplo_ == Teuchos::UPPER_TRI) {
              triUplo = "UPPER_TRI";
            }
            std::string unitDiag ("INVALID");
            if (unit_diag_ == Teuchos::NON_UNIT_DIAG) {
              unitDiag = "NON_UNIT_DIAG";
            }
            else if (unit_diag_ == Teuchos::UNIT_DIAG) {
              unitDiag = "UNIT_DIAG";
            }

            out << "numRows_ = " << numRows_ << endl
                << "numCols_ = " << numCols_ << endl
                << "isEmpty_ = " << isEmpty_ << endl
                << "tri_uplo_ = " << triUplo << endl
                << "unit_diag_ = " << unitDiag << endl;
            if (big_ptrs_.size() > 0) {
              out << "numEntries = " << big_ptrs_[big_ptrs_.size()-1] << endl;
            }
            else if (sml_ptrs_.size() > 0) {
              out << "numEntries = " << sml_ptrs_[sml_ptrs_.size()-1] << endl;
            }
            else {
              out << "numEntries = 0" << endl;
            }

            if (includesVerbLevel (vl, VERB_EXTREME)) { // vl >= VERB_EXTREME
              // Only print out all the sparse matrix's data in
              // extreme verbosity mode.
              out << "ptrs = [";
              if (big_ptrs_.size() > 0) {
                std::copy (big_ptrs_.begin(), big_ptrs_.end(),
                           std::ostream_iterator<Ordinal> (out, " "));
              }
              else {
                std::copy (sml_ptrs_.begin(), sml_ptrs_.end(),
                           std::ostream_iterator<Ordinal> (out, " "));
              }
              out << "]" << endl << "inds_ = [";
              std::copy (inds_.begin(), inds_.end(),
                         std::ostream_iterator<Ordinal> (out, " "));
              out << "]" << endl << "vals_ = [";
              std::copy (vals_.begin(), vals_.end(),
                         std::ostream_iterator<Scalar> (out, " "));
              out << "]" << endl;
            } // vl >= VERB_EXTREME
          } // if is initialized
        } // vl >= VERB_MEDIUM
      } // vl >= VERB_LOW
    }

    /// \brief Convert to dense matrix and return.
    ///
    /// \warning This method is for debugging only.  It uses a lot of
    ///   memory.  Users should never call this method.  Do not rely
    ///   on this method continuing to exist in future releases.
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, scalar_type> >
    asDenseMatrix () const
    {
      using Teuchos::ArrayRCP;
      using Teuchos::RCP;
      using Teuchos::rcp;
      typedef Teuchos::OrdinalTraits<ordinal_type> OTO;
      typedef Teuchos::ScalarTraits<scalar_type> STS;
      typedef Teuchos::SerialDenseMatrix<int, scalar_type> dense_matrix_type;

      RCP<dense_matrix_type> A_ptr =
        rcp (new dense_matrix_type (numRows_, numCols_));
      dense_matrix_type& A = *A_ptr; // for notational convenience

      if (big_ptrs_.size() > 0) {
        ArrayRCP<const size_t> ptr = big_ptrs_;
        for (ordinal_type i = OTO::zero(); i < numRows_; ++i) {
          for (size_t k = ptr[i]; k < ptr[i+1]; ++k) {
            const ordinal_type j = inds_[k];
            const scalar_type A_ij = vals_[k];
            A(i,j) += A_ij;
          }
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            // Respect whatever is in the sparse matrix, even if it is wrong.
            // This is helpful for debugging.
            A(i,i) += STS::one ();
          }
        }
      }
      else if (sml_ptrs_.size() > 0) {
        ArrayRCP<const ordinal_type> ptr = sml_ptrs_;
        for (ordinal_type i = OTO::zero(); i < numRows_; ++i) {
          for (ordinal_type k = ptr[i]; k < ptr[i+1]; ++k) {
            const ordinal_type j = inds_[k];
            const scalar_type A_ij = vals_[k];
            A(i,j) += A_ij;
          }
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            // Respect whatever is in the sparse matrix, even if it is wrong.
            // This is helpful for debugging.
            A(i,i) += STS::one ();
          }
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Kokkos::DefaultHost"
          "SparseOps::asDenseMatrix: both big_ptrs_ and sml_ptrs_ are empty.  "
          "Please report this bug to the Kokkos developers.");
      }
      return A_ptr;
    }

    //@}
    //! \name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
    {
      return Allocator::allocRowPtrs(node,numEntriesPerRow);
    }

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
    {
      return Allocator::template allocStorage<T,size_t>(node,rowPtrs);
    }

    //! Finalize a graph
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    ///
    /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
    ///   else LOWER_TRI if the matrix is lower triangular.
    ///
    /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
    ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
    ///
    void setGraphAndMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph,
                           const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &matrix);

    //@}
    //! @name Computational methods
    //@{

    /// \brief Y := alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, overwriting Y with the result.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \note This method does not respect the implicit unit diagonal
    ///   indication.  If you want to simulate having an implicitly
    ///   stored unit diagonal for the operation Y := A*X, you must
    ///   compute Y := X + A*X instead.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector. Contents will be overwritten.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := beta * Y + alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, accumulating the result into Y.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \note This method does not respect the implicit unit diagonal
    ///   indication.  If you want to simulate having an implicitly
    ///   stored unit diagonal for the operation Y := A*X, you must
    ///   compute Y := X + A*X instead.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param beta [in] Scalar constant \f$\beta\f$ by which to
    ///   multiply Y when summing with the result of the sparse
    ///   matrix-(multi)vector multiply.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
    ///
    /// Solve the (upper or lower) triangular system Y = Op(A) X.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to solve with the matrix, its
    ///   transpose, or its conjugate transpose (if applicable).
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultHostSparseOps(const DefaultHostSparseOps& source);

    void getOffsets(ArrayRCP<const size_t> &ptrs) const {
      ptrs = big_ptrs_;
    }

    void getOffsets(ArrayRCP<const Ordinal> &ptrs) const {
      ptrs = sml_ptrs_;
    }

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void solvePrivate(Teuchos::ETransp trans,
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector< RangeScalar,Node> &X) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                               MultiVector<RangeScalar,Node> &Y) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                         RangeScalar beta,
                         MultiVector<RangeScalar,Node> &Y) const;

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    // packed CRS: array of row pointers, array of indices, array of values
    // pointers are EITHER size_t or Ordinal; one will be null
    ArrayRCP<const Ordinal> inds_;
    ArrayRCP<const size_t>  big_ptrs_;
    ArrayRCP<const Ordinal> sml_ptrs_;
    ArrayRCP<const Scalar>  vals_;

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;

    Ordinal numRows_;
    Ordinal numCols_;
    bool isInitialized_;
    bool isEmpty_;
  };


  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params)
  {
    graph.setMatDesc(uplo,diag);
    std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeGraph(graph,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    // determine how many non-zeros, so that we can decide whether to reduce the offset pointer type
    ArrayRCP<const size_t> bigptrs = graph.getPointers();
    const size_t numrows = bigptrs.size() - 1,
                   numnz = bigptrs[numrows];
    if (numnz < (size_t)Teuchos::OrdinalTraits<Ordinal>::max()) {
      ArrayRCP<Ordinal> smallptrs = arcp<Ordinal>(numrows+1);
      std::copy( bigptrs.begin(), bigptrs.end(), smallptrs.begin() );
      graph.setSmallPointers(smallptrs);
    }
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // nothing much to do here
    std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // finalize them individually; no benefit to doing them together
    finalizeGraph(uplo,diag,graph,params);
    finalizeMatrix(graph,matrix,params);
  }


  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::DefaultHostSparseOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , numCols_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::
  DefaultHostSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params)
  : node_(node)
  , numRows_(0)
  , numCols_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    (void) params; // Not using this yet.

    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::~DefaultHostSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  RCP<Node> DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::setGraphAndMatrix(
              const RCP<const DefaultCrsGraph<Ordinal,Node> > &opgraph,
              const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " operators already initialized.");
    numRows_ = opgraph->getNumRows();
    numCols_ = opgraph->getNumCols();
    if (opgraph->isEmpty() || numRows_ == 0 || numCols_ == 0) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      big_ptrs_ = opgraph->getPointers();
      sml_ptrs_ = opgraph->getSmallPointers();
      inds_ = opgraph->getIndices();
      vals_ = opmatrix->getValues();
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        big_ptrs_ != null && sml_ptrs_ != null,
        std::logic_error, " Internal logic error: graph has small and big pointers. Please notify Kokkos team."
      )
      const size_t lenptrs = (big_ptrs_ != null ? big_ptrs_.size() : sml_ptrs_.size());
      // these checks just about the most that we can perform
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          lenptrs != (size_t)numRows_+1 || inds_.size() != vals_.size(),
          std::runtime_error, " matrix and graph seem incongruent.");
    }
    opgraph->getMatDesc( tri_uplo_, unit_diag_ );
    isInitialized_ = true;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::solvePrivate(Teuchos::ETransp trans,
                                                        const MultiVector<DomainScalar,Node> &Y,
                                                              MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    typedef DefaultSparseSolveOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar>   Op;
    typedef DefaultSparseTransposeSolveOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar>  TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(unit_diag_ != Teuchos::UNIT_DIAG, std::runtime_error,
          " solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.offs    = rbh.template addConstBuffer< OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.offs    = rbh.template addConstBuffer< OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::solve(Teuchos::ETransp trans,
                                                        const MultiVector<DomainScalar,Node> &Y,
                                                              MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " this solve was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumCols() != (size_t)Y.getNumCols(),
        std::runtime_error, " Left hand side and right hand side multivectors have differing numbers of vectors."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumRows() < (size_t)numRows_,
        std::runtime_error, " Left-hand-side multivector does not have enough rows. "
                            "Likely cause is that the column map was not provided to "
                            "the Tpetra::CrsMatrix in the case of an implicit unit diagonal."
    );
    if (big_ptrs_ != null) {
      solvePrivate<DomainScalar,RangeScalar,size_t>(trans,Y,X);
    }
    else {
      solvePrivate<DomainScalar,RangeScalar,Ordinal>(trans,Y,X);
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiplyPrivate(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X,
                                      MultiVector<RangeScalar ,Node> &Y) const
  {
    // the 1 template parameter below means that beta is not used in computations
    // and the output multivector enjoys overwrite semantics (i.e., will overwrite data/NaNs in Y)
    typedef DefaultSparseMultiplyOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 1>  Op;
    typedef DefaultSparseTransposeMultiplyOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 1> TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X,
                                      MultiVector<RangeScalar ,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    if (big_ptrs_ != null) {
      multiplyPrivate<DomainScalar,RangeScalar,size_t>(trans,alpha,X,Y);
    }
    else {
      multiplyPrivate<DomainScalar,RangeScalar,Ordinal>(trans,alpha,X,Y);
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiplyPrivate(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    // the 0 template parameter below means that beta is used in computations
    // and the output multivector enjoys accumulation semantics (i.e., will not overwrite data/NaNs in Y)
    typedef DefaultSparseMultiplyOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 0>  Op;
    typedef DefaultSparseTransposeMultiplyOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 0> TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    if (big_ptrs_ != null) {
      multiplyPrivate<DomainScalar,RangeScalar,size_t>(trans,alpha,X,beta,Y);
    }
    else {
      multiplyPrivate<DomainScalar,RangeScalar,Ordinal>(trans,alpha,X,beta,Y);
    }
  }

  /// \brief Partial specialization of DefaultHostSparseOps for Scalar=void.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  /// \tparam Allocator The allocator type.
  ///
  /// \warning This partial specialization is _not_ for users.  Kokkos
  ///   developers should see the discussion below explaining why we
  ///   need this partial specialization.
  ///
  /// Developer documentation
  /// =======================
  ///
  /// We include a partial specialization as a work-around for a
  /// Windows MSVC compilation problem reported on 08 Aug 2012 by
  /// Brent Perschbacher.  The issue is that MSVC is attempting to
  /// compile the generic methods for Scalar=void, since we do refer
  /// to the type for Scalar=void in e.g., Tpetra::CrsGraph.  However,
  /// whenever we refer to the Scalar=void case, we only reference the
  /// typedefs and inner classes inside, not the methods.  Other
  /// compilers do not attempt to compile methods of a template class
  /// that aren't called; MSVC apparently does.
  ///
  /// Kokkos developers must imitate DefaultHostSparseOps by providing
  /// their own partial specializations of their local sparse kernels
  /// classes for the Scalar=void case.
  ///
  /// gcc 4.5.1 says that "default template arguments may not be used
  /// in partial specializations," so we aren't allowed to specify a
  /// default Allocator.
  template <class Ordinal, class Node, class Allocator>
  class DefaultHostSparseOps<void, Ordinal, Node, Allocator> : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef void scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;
    //! The type of this object, the sparse operator object
    typedef DefaultHostSparseOps<void, Ordinal, Node, Allocator> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// scalar_type.  (In fact, it _should_ be different than
    /// scalar_type=void in this case.  The intended use case is to
    /// rebind scalar_type=void to the different scalar type S2 of a
    /// matrix.)
    ///
    /// This always specifies a specialization of \c
    /// DefaultHostSparseOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c scalar_type.
    template <class S2>
    struct bind_scalar {
      typedef DefaultHostSparseOps<S2,Ordinal,Node,Allocator> other_type;
    };

    /// \brief Sparse operations type for a different ordinal type.
    ///
    /// The bind_ordinal struct defines the type responsible for
    /// sparse operations for an ordinal type O2, which may be
    /// different from Ordinal.
    ///
    /// This always specifies a specialization of \c
    /// DefaultHostSparseOps, regardless of the ordinal type O2.
    ///
    /// \tparam O2 An ordinal type possibly different from \c Ordinal.
    template <class O2>
    struct bind_ordinal {
      typedef DefaultHostSparseOps<void, O2, Node, Allocator> other_type;
    };

    //@}
    //! \name Constructors/Destructor
    //@{

    //! Constructor that takes a Kokkos Node: DO NOT CALL (Scalar=void specialization).
    DefaultHostSparseOps (const RCP<Node> &node) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Someone attempted to "
        "instantiate Kokkos::DefaultHostSparseOps with Scalar=void.  "
        "This is not allowed.  "
        "The Scalar=void specialization exists only for its typedefs.  "
        "Please report this bug to the Kokkos developers.");
    }

    //! Constructor that takes a Kokkos Node and parameters: DO NOT CALL (Scalar=void specialization).
    DefaultHostSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Someone attempted to "
        "instantiate Kokkos::DefaultHostSparseOps with Scalar=void.  "
        "This is not allowed.  "
        "The Scalar=void specialization exists only for its typedefs.  "
        "Please report this bug to the Kokkos developers.");
    }

    //! Destructor.
    ~DefaultHostSparseOps() {
      // We don't throw an exception here, because throwing exceptions
      // in a destructor may cause the application to terminate [1].
      // However, it's impossible that execution will reach this
      // point, since all the constructors throw exceptions.
      //
      // [1] http://www.parashift.com/c++-faq/dtors-shouldnt-throw.html
    }

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os <<  "Kokkos::DefaultHostSparseOps<"
         << "Scalar=void"
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::includesVerbLevel;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;

      // Interpret the default verbosity level as VERB_MEDIUM.
      const Teuchos::EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
        out << this->description() << std::endl;
      }
    }

    //@}
    //! \name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const {
      // You're not supposed to instantiate this object, so we always
      // return null here.
      return Teuchos::null;
    }

    //@}
    //! @name Initialization of graph and matrix
    //@{

    /// \brief Allocate and initialize the storage for the row offsets.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   allocating its row offsets.  Since it's a class method, we
    ///   may call it without needing to instantiate a
    ///   DefaultHostSparseOps instance.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
    {
      return Allocator::allocRowPtrs(node,numEntriesPerRow);
    }

    /// \brief Allocate and initialize the storage for graph or matrix storage.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   allocating the column indices (T=Ordinal).  Since it's a
    ///   class method, we may call it without needing to instantiate
    ///   a DefaultHostSparseOps instance.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
    {
      return Allocator::template allocStorage<T,size_t>(node,rowPtrs);
    }

    /// \brief Finalize a graph.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   finalizing the graph structure.  Since it's a class method,
    ///   we may call it without needing to instantiate a
    ///   DefaultHostSparseOps instance.
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params) {
      using Teuchos::ArrayRCP;
      using Teuchos::as;

      std::string FuncName("Kokkos::DefaultHostSparseOps::finalizeGraph(graph,params)");

      graph.setMatDesc(uplo,diag);
      TEUCHOS_TEST_FOR_EXCEPTION(graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized.");

      // determine how many non-zeros, so that we can decide whether to reduce the offset pointer type
      ArrayRCP<const size_t> bigptrs = graph.getPointers();
      const size_t numrows = bigptrs.size() - 1,
        numnz = bigptrs[numrows];
      if (numnz < as<size_t> (Teuchos::OrdinalTraits<Ordinal>::max())) {
        ArrayRCP<Ordinal> smallptrs = arcp<Ordinal>(numrows+1);
        std::copy( bigptrs.begin(), bigptrs.end(), smallptrs.begin() );
        graph.setSmallPointers(smallptrs);
      }
    }

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    DefaultHostSparseOps(const DefaultHostSparseOps& source);
  };

  /** \example CrsMatrix_DefaultMultiplyTests.hpp
    * This is an example that unit tests and demonstrates the implementation requirements for the DefaultSparseOps class.
    */

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */

