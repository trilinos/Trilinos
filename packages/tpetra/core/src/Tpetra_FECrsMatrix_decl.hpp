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
/// If you want to use Tpetra::FECrsMatrix, include
/// "Tpetra_FECrsMatrix.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::FECrsMatrix,
/// include this file (Tpetra_FECrsMatrix_decl.hpp).

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_FECrsGraph.hpp"

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
  friend class CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
public:
    //! @name Typedefs
    //@{

    /// \brief This class' first template parameter; the type of each
    ///        entry in the matrix.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::scalar_type scalar_type;

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
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_ordinal_type local_ordinal_type;

    //! This class' third template parameter; the type of global indices.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::global_ordinal_type global_ordinal_type;

    //! This class' fourth template parameter; the Kokkos device type.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::node_type node_type;

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
    ///        the part of the sparse matrix on each MPI process.
    typedef typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;

    /// \brief Parent CrsMatrix type using the same scalars
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    typedef FECrsGraph<LocalOrdinal, GlobalOrdinal, Node> fe_crs_graph_type;

    //@}
    //! @name Constructors and destructor
    //@{


    /// \brief Constructor specifying one or two previously constructed graphs.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a FECrsMatrix with this
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
    /// create a FECrsMatrix by accident when passing a FECrsGraph into a
    /// function that takes a FECrsMatrix.
    ///
    /// \param graph [in] The FECrsGraph structure of the sparse matrix.
    ///   The graph <i>must</i> have endFill() called
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit FECrsMatrix (const Teuchos::RCP<const fe_crs_graph_type>& graph,
                          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor.
    virtual ~FECrsMatrix () {}

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
    void globalAssemble () {endFill();}

    //! Migrates data to the owned mode
    void endFill();

    //! Activates the owned+shared mode for assembly
    void beginFill();

  protected:
    /// \brief Migrate data from the owned+shared to the owned matrix
    /// Since this is non-unique -> unique, we need a combine mode.
    /// Precondition: Must be FE_ACTIVE_OWNED_PLUS_SHARED mode
    void doOwnedPlusSharedToOwned(const CombineMode CM=Tpetra::ADD);

    /// \brief Migrate data from the owned to the owned+shared matrix
    /// Precondition: Must be FE_ACTIVE_OWNED mode
    void doOwnedToOwnedPlusShared(const CombineMode CM=Tpetra::ADD);

    //! Switches which CrsGraph is active (without migrating data)
    void switchActiveCrsMatrix();
    //@}
  
  private:
    // We forbid assignment (operator=) and copy construction
    FECrsMatrix<Scalar,LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;

    FECrsMatrix (const FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = delete;




    // Enum for activity
    enum FEWhichActive
    {
      FE_ACTIVE_OWNED,
      FE_ACTIVE_OWNED_PLUS_SHARED
    };

    // The FECrsGraph from construction time
    Teuchos::RCP<const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node> > feGraph_;

    // This is whichever multivector isn't currently active
    Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > inactiveCrsMatrix_;
    // This is in RCP to make shallow copies of the FECrsMatrix work correctly
    Teuchos::RCP<FEWhichActive> activeCrsMatrix_;

};    // end class FECrsMatrix



}     // end namespace Tpetra



#endif // TPETRA_FECRSMATRIX_DECL_HPP
