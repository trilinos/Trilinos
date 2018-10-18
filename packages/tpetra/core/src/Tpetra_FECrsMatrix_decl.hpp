// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_FECRSMATRIX_DECL_HPP
#define TPETRA_FECRSMATRIX_DECL_HPP

/// \file Tpetra_FECrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::FECrsMatrix class
///
/// If you want to use Tpetra::CrsMatrix, include
/// "Tpetra_FECrsMatrix.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::CrsMatrix,
/// include this file (Tpetra_CrsMatrix_decl.hpp).

#include "Tpetra_CrsMatrix_decl.hpp"


namespace Tpetra {


// \class FECrsMatrix
// \brief Sparse matrix that presents a row-oriented interface that lets
//        users read or modify entries.
template<class Scalar        = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class FECrsMatrix :
   public CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{

    private:
        friend class CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    public:

    //! @name Typedefs
    //@{

    /// \brief This class' first template parameter; the type of each
    ///   entry in the matrix.
    typedef Scalar scalar_type;

    /// \brief The type used internally in place of \c Scalar.
    ///
    /// Some \c Scalar types might not work with Kokkos on all
    /// execution spaces, due to missing CUDA device macros or
    /// volatile overloads.  The C++ standard type std::complex<T> has
    /// this problem.  To fix this, we replace std::complex<T> values
    /// internally with the (usually) bitwise identical type
    /// Kokkos::complex<T>.  The latter is the \c impl_scalar_type
    /// corresponding to \c Scalar = std::complex.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type impl_scalar_type;

    //! This class' second template parameter; the type of local indices.
    typedef LocalOrdinal local_ordinal_type;

    //! This class' third template parameter; the type of global indices.
    typedef GlobalOrdinal global_ordinal_type;

    //! This class' fourth template parameter; the Kokkos device type.
    typedef Node node_type;

    //! The Kokkos device type.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::device_type device_type;

    //! The Kokkos execution space.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::execution_space execution_space;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>Scalar</tt>, but may differ for
    /// certain <tt>Scalar</tt> types.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type mag_type;

    //! The Map specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type map_type;

    //! The Import specialization suitable for this CrsMatrix specialization
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::import_type import_type;

    //! The Export specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::export_type export_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_graph_type crs_graph_type;

    //! The part of the sparse matrix's graph on each MPI process.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_graph_type local_graph_type;

    /// \brief The specialization of Kokkos::CrsMatrix that represents
    ///   the part of the sparse matrix on each MPI process.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;

    /// \brief Parent CrsMatrix type using the same scalars
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

    //@}
    //! @name Constructors and destructor
    //@{


    /// \brief Constructor specifying one or two previously constructed graphs.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graphs must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param graph [in] The graph structure of the on-rank sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param offRankGraph [in] The graph structure of the off-rank sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit FECrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
              const Teuchos::RCP<const crs_graph_type>& offRankGraph,
                          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor.
    virtual ~CrsMatrix ();

    //@}
    //! @name Methods for inserting, modifying, or removing entries
    //@{


    /// \brief Sum into one or more sparse matrix entries, using
    ///   global indices.
    ///
    /// This is a local operation; it does not involve communication.
    /// However, if you sum into rows not owned by the calling
    /// process, it may result in future communication in
    /// globalAssemble() (which is called by fillComplete()).
    ///
    /// If \c globalRow is owned by the calling process, then this
    /// method performs the sum-into operation right away.  Otherwise,
    /// if the row is <i>not</i> owned by the calling process, this
    /// method defers the sum-into operation until globalAssemble().
    /// That method communicates data for nonowned rows to the
    /// processes that own those rows.  Then, globalAssemble() does
    /// one of the following:
    /// <ul>
    /// <li> It calls insertGlobalValues() for that data if the matrix
    ///      has a dynamic graph. </li>
    /// <li> It calls sumIntoGlobalValues() for that data if the matrix
    ///      has a static graph.  The matrix silently ignores
    ///      (row,column) pairs that do not exist in the graph.
    /// </ul>
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param cols [in] One or more column indices.
    /// \param vals [in] One or more values corresponding to those
    ///   column indices.  <tt>vals[k]</tt> corresponds to
    ///   <tt>cols[k]</tt>.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceGlobalValues() (which see).
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                         const Teuchos::ArrayView<const Scalar>& vals,
                         const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Epetra compatibility version of sumIntoGlobalValues
    ///   (see above), that takes input as raw pointers instead of
    ///   Kokkos::View.
    ///
    /// Arguments are the same and in the same order as those of
    /// Epetra_CrsMatrix::SumIntoGlobalValues, except for \c atomic,
    /// which is as above.
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param numEnt [in] Number of valid entries in \c vals and
    ///   \c cols.  This has type \c LocalOrdinal because we assume
    ///   that users will never want to insert more column indices
    ///   in one call than the matrix has columns.
    /// \param vals [in] \c numEnt values corresponding to the column
    ///   indices in \c cols.  That is, \c vals[k] is the value
    ///   corresponding to \c cols[k].
    /// \param cols [in] \c numEnt global column indices.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const LocalOrdinal numEnt,
                         const Scalar vals[],
                         const GlobalOrdinal cols[],
                         const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Sum into one or more sparse matrix entries, using global
    ///   row and column indices.
    ///
    /// For global row index \c globalRow and globaal column indices
    /// <tt>cols</tt>, perform the update <tt>A(globalRow, cols[k]) +=
    /// vals[k]</tt>.  The row index and column indices must be valid
    /// on the calling process, and all matrix entries <tt>A(globalRow,
    /// cols[k])</tt> must already exist.  (This method does
    /// <i>not</i> change the matrix's structure.)  If the row index
    /// is valid, any invalid column indices are ignored, but counted
    /// in the return value.
    ///
    /// This overload of the method takes the column indices and
    /// values as Kokkos::View.
    ///
    /// \tparam GlobalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of GlobalOrdinal.
    /// \tparam ImplScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of impl_scalar_type (usually the same as Scalar,
    ///   unless Scalar is std::complex<T> for some T, in which case
    ///   it is Kokkos::complex<T>).
    ///
    /// \param globalRow [in] Global index of a row.  This row
    ///   <i>must</i> be owned by the calling process.
    /// \param cols [in] Global indices of the columns whose entries we
    ///   want to modify.
    /// \param vals [in] Values corresponding to the above column
    ///   indices.  <tt>vals(k)</tt> corresponds to <tt>cols(k)</tt>.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceGlobalValues() (which see).
    template<class GlobalIndicesViewType,
         class ImplScalarViewType>
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
            const typename UnmanagedView<GlobalIndicesViewType>::type& inputInds,
            const typename UnmanagedView<ImplScalarViewType>::type& inputVals,
            const bool atomic = useAtomicUpdatesByDefault) const;

    //! Set all matrix entries equal to \c alpha.
    void setAllToScalar (const Scalar& alpha);

    //@}
    //! @name Transformational methods
    //@{

    /// \brief Communicate nonlocal contributions to other processes.
    ///
    /// Users do not normally need to call this method.  fillComplete
    /// always calls this method, unless you specifically tell
    /// fillComplete to do otherwise by setting its "No Nonlocal
    /// Changes" parameter to \c true.  Thus, it suffices to call
    /// fillComplete.
    ///
    /// Methods like insertGlobalValues and sumIntoGlobalValues let
    /// you add or modify entries in rows that are not owned by the
    /// calling process.  These entries are called "nonlocal
    /// contributions."  The methods that allow nonlocal contributions
    /// store the entries on the calling process, until globalAssemble
    /// is called.  globalAssemble sends these nonlocal contributions
    /// to the process(es) that own them, where they then become part
    /// of the matrix.
    ///
    /// This method only does global assembly if there are nonlocal
    /// entries on at least one process.  It does an all-reduce to
    /// find that out.  If not, it returns early, without doing any
    /// more communication or work.
    ///
    /// If you previously inserted into a row which is not owned by
    /// <i>any</i> process in the row Map, the behavior of this method
    /// is undefined.  It may detect the invalid row indices and throw
    /// an exception, or it may silently drop the entries inserted
    /// into invalid rows.  Behavior may vary, depending on whether
    /// Tpetra was built with debug checking enabled.
    void globalAssemble();

    /// \brief Resume operations that may change the values or
    ///   structure of the matrix.
    ///
    /// This method must be called as a collective operation.
    ///
    /// Calling fillComplete "freezes" both the values and the
    /// structure of the matrix.  If you want to modify the matrix
    /// again, you must first call resumeFill.  You then may not call
    /// resumeFill again on that matrix until you first call
    /// fillComplete.  You may make sequences of fillComplete,
    /// resumeFill calls as many times as you wish.
    ///
    /// \post <tt>isFillActive() && ! isFillComplete()</tt>
    void resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


  private:
    // We forbid copy construction by declaring this method private
    // and not implementing it.
    FECrsMatrix (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs);

    // We forbid assignment (operator=) by declaring this method
    // private and not implementing it.
    FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs);


    ///! Off-rank graph
    Teuchos::RCP<Graph> offRankGraph_;
    //@}

    //! The offRank local sparse matrix.
    Teuchos::RCP<crs_matrix_type> offRankMatrix_;


}; // end class FECrsMatrix


};  // end namespace Tpetra

#endif // TPETRA_FECRSMATRIX_DECL_HPP
