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

#ifndef KOKKOS_CRSGRAPH_HPP
#define KOKKOS_CRSGRAPH_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  /// \class CrsGraph
  /// \brief A default compressed-row sparse graph.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Ordinal Same as the LocalOrdinal template parameter of
  ///   Tpetra objects.
  /// \tparam Node Kokkos Node type; same as the Node template
  ///   parameter of Tpetra objects.
  template <class Ordinal,
            class Node>
  class CrsGraph {
  public:
    typedef Ordinal               OrdinalType;
    typedef Node                  NodeType;

    //! @name Constructors/Destructor
    //@{

    //! Default constuctor.
    CrsGraph (size_t numRows, const RCP<Node> &node);

    //! Destructor.
    virtual ~CrsGraph();

    //@}
    //! @name Accessor routines.
    //@{

    //! Node accessor.
    RCP<Node> getNode() const;

    //@}
    //! @name Data entry and accessor methods.
    //@{

    //! Return the number of rows in the graph.
    size_t getNumRows() const;

    //! Return the number of entries in the graph.
    size_t getNumEntries() const;

    //! Indicates that the graph is filled, but empty.
    bool isEmpty() const;

    //! Indicatest that the graph has been finalized.
    bool isFinalized() const;

    /// \brief Indicate that the structure is 1D.
    ///
    /// It will never be the case that both is1DStructure() and
    /// is2DStructure() return true.
    bool is1DStructure() const;

    /// \brief Indicate that the structure is 2D.
    ///
    /// It will never be the case that both is1DStructure() and
    /// is2DStructure() return true.
    bool is2DStructure() const;

    //! Indicate that the stucture is optimized.
    bool isOptimized() const;

    /// \brief Submit the indices and offset for 1D storage.
    ///
    /// \post is1DStructure() == true
    void set1DStructure(ArrayRCP<Ordinal> inds,
                        ArrayRCP<size_t>  rowBegs,
                        ArrayRCP<size_t>  rowEnds);

    /// Submit the indices for 2D storage.
    ///
    /// \post is2DStructure() == true
    void set2DStructure(ArrayRCP<ArrayRCP<Ordinal> > inds,
                        ArrayRCP<size_t> numEntriesPerRow);

    //! Retrieve the structure for 1D storage.
    /**
          If is1DStructure() == false, then
          \post inds == rowBegs == rowEnds == null

          Otherwise,
          \post indices for row \c r are inds[r], where \f$r \in [b,e)\f$, where \f$b = rowBegs[r]\f$ and \f$e = rowEnds[r]\f$
          \post rowBegs has getNumRows()+1 entries; the last entry is inds.size()
          \post rowEnds has getNumRows() entries
     */
    void get1DStructure(ArrayRCP<Ordinal> &inds,
                        ArrayRCP<size_t>  &rowBegs,
                        ArrayRCP<size_t>  &rowEnds);

    /// \brief Finalize storage for the graph.
    ///
    /// Instruct the graph to perform any necessary manipulation,
    /// including (optionally) optimizing the storage of the graph
    /// data.
    ///
    /// \param OptimizeStorage [in] If true, permit the graph to
    ///   reallocate storage on the host in order to provide optimal
    ///   storage and/or performance.
    ///
    /// \post if OptimizeStorage == true, then is2DStructure() == true
    ///   on return.
    void finalize (bool OptimizeStorage);

    /// \brief Finalize storage for the graph with associated matrix values.
    ///
    /// This is typically called from a CrsMatrix.  It performs the
    /// finalize for the graph and the matrix at the same time, so the
    /// matrix doesn't have to.  In that case, the Scalar template
    /// parameter here should be the same as CrsMatrix's Scalar
    /// template parameter.
    ///
    /// \param OptimizeStorage [in] If true, permit the graph to
    ///   reallocate storage on the host in order to provide optimal
    ///   storage and/or performance.
    ///
    /// \param values2D [in/out] 2D-structured matrix values. Required
    ///   to be nonnull on input if is2DStructure() is true. Set to
    ///   null on output if OptimizeStorage is true.
    ///
    /// \param values1D [in/out] 1D-structured matrix values. Required
    ///   to be nonnull on input if is1DStructure() is true. Allocated
    ///   on output if OptimizeStorage is true.
    ///
    /// \post if OptimizeStorage == true or already is2DStructure(),
    ///   then is2DStructure() == true.
    template <class Scalar>
    void finalize (bool OptimizeStorage,
                   ArrayRCP<ArrayRCP<Scalar> >& values2D,
                   ArrayRCP<Scalar>& values1D);

    //! Release data associated with this graph.
    virtual void clear ();

    //@}

  protected:
    //! Copy constructor (protected and not implemented)
    CrsGraph(const CrsGraph& sources);

    RCP<Node> node_;
    size_t numRows_, numEntries_;
    bool isFinalized_, isEmpty_, is1D_, is2D_, isOpt_;

    // 2D storage
    ArrayRCP<ArrayRCP<Ordinal> >  indices2D_;
    ArrayRCP<size_t>              numEntriesPerRow_;
    // 1D storage
    ArrayRCP<Ordinal>             indices1D_;
    ArrayRCP<size_t>              rowBegs_, rowEnds_;
  };


  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraph<Ordinal,Node>::CrsGraph(size_t numRows, const RCP<Node> &node) 
  : node_(node)
  , numRows_(numRows)
  {
    CrsGraph<Ordinal,Node>::clear();
  }

  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraph<Ordinal,Node>::~CrsGraph() {
  }

  // ======= clear ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::clear() {
    isFinalized_   = false;
    isEmpty_       = false;
    is1D_          = false;
    is2D_          = false;
    isOpt_         = false;
    numEntries_    = 0;
    indices2D_        = null;
    numEntriesPerRow_ = null;
    rowBegs_          = null;
    rowEnds_          = null;
    indices1D_        = null;
  }

  // ======= node ===========
  template <class Ordinal, class Node>
  RCP<Node> CrsGraph<Ordinal,Node>::getNode() const {
    return node_;
  }

  // ======= numrows ===========
  template <class Ordinal, class Node>
  size_t CrsGraph<Ordinal,Node>::getNumRows() const {
    return numRows_;
  }

  // ======= numentries ===========
  template <class Ordinal, class Node>
  size_t CrsGraph<Ordinal,Node>::getNumEntries() const {
    return numEntries_;
  }

  // ======= isempty ===========
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::isEmpty() const {
    return isEmpty_;
  }

  // ======= isfinalized ===========
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::isFinalized() const {
    return isFinalized_;
  }

  // ======= is1d ===========
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::is1DStructure() const {
    return is1D_;
  }

  // ======= is2d ===========
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::is2DStructure() const {
    return is2D_;
  }

  // ======= isopt ===========
  template <class Ordinal, class Node>
  bool CrsGraph<Ordinal,Node>::isOptimized() const {
    return isOpt_;
  }

  // ======= get 1d ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::get1DStructure(ArrayRCP<Ordinal> &inds, 
                                                         ArrayRCP<size_t>  &rowBegs,
                                                         ArrayRCP<size_t>  &rowEnds)
  {
    inds = indices1D_;
    rowBegs = rowBegs_;
    rowEnds = rowEnds_;
  }

  // ======= get 2d ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::get2DStructure(ArrayRCP<ArrayRCP<Ordinal> > &inds,
                                                         ArrayRCP<size_t> &numEntriesPerRow) 
  {
    inds = indices2D_;
    numEntriesPerRow = numEntriesPerRow_;
  }

  // ======= set 1d ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::set1DStructure(ArrayRCP<Ordinal> inds, 
                                                         ArrayRCP<size_t>  rowBegs,
                                                         ArrayRCP<size_t>  rowEnds)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t)rowBegs.size() != numRows_+1 || (size_t)rowEnds.size() != numRows_,
      std::runtime_error, Teuchos::typeName(*this) << "::set1DStructure(inds,"
      "rowBegs,rowEnds): rowBegs and rowEnds do not have the right sizes.  "
      "Either rowBegs.size() = " << rowBegs.size() << " != numRows_+1 = " << (numRows_+1)
      << ", or rowEnds.size() = " << rowEnds.size() << " != numRows_ = " << numRows_ << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t)rowBegs[numRows_] > (size_t)inds.size(), std::runtime_error,
      Teuchos::typeName(*this) << "::set1DStructure(inds,rowBegs,rowEnds): "
      "rowBegs[numRows_=" << numRows_ << "] = " << rowBegs[numRows_]
      << " > inds.size() = " << inds.size() << ".");
    this->clear();

    indices1D_ = inds;
    rowBegs_ = rowBegs;
    rowEnds_ = rowEnds;
    if (numRows_ > 0) {
      for (size_t i=0; i < this->getNumRows(); ++i) {
        numEntries_ += (this->rowEnds_[i] - this->rowBegs_[i]);
#ifdef HAVE_KOKKOS_DEBUG
        // row i goes like [ begs[i] , ends[i] )
        // sanity        : begs[i] <= ends[i]
        // ordering      : begs[i] <= begs[i+1]
        // no overlapping: ends[i] <= begs[i+1]
        TEUCHOS_TEST_FOR_EXCEPTION(
          rowBegs_[i+1] < rowBegs_[i] || rowEnds_[i] < rowBegs_[i] || rowEnds_[i] > rowBegs_[i+1],
          std::runtime_error, Teuchos::typeName(*this) << "::set1DStructure("
          "inds,rowBegs,rowEnds): ends and begs are not consistent.");
#endif
      }
    }
    is1D_ = true;
  }

  // ======= set 2d ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::set2DStructure(ArrayRCP<ArrayRCP<Ordinal> > inds,
                                                         ArrayRCP<size_t> numEntriesPerRow)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t)inds.size() != numRows_ || (size_t)numEntriesPerRow.size() != numRows_,
      std::runtime_error, Teuchos::typeName(*this) << "::set2DStructure(inds,"
      "numEntriesPerRow): numEntriesPerRow and inds must have as many entries "
      "as the number of rows given to the constructor.");
    this->clear ();

    indices2D_  = inds;
    if (indices2D_ != null) {
      numEntriesPerRow_ = numEntriesPerRow;
      numEntries_ = std::accumulate (this->numEntriesPerRow_.begin(), this->numEntriesPerRow_.end(), 0);
#ifdef HAVE_KOKKOS_DEBUG
      for (size_t i = 0; i < numRows_; ++i) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          (size_t)inds[i].size() < numEntriesPerRow[i], std::runtime_error,
          Teuchos::typeName(*this) << "::set2DStructure(): inds[" << i
          << "] == " << inds[i] << " is not large enough for the specified "
          "number of entries, " << " numEntriesPerRow[" << i << "] == "
          << numEntriesPerRow[i]);
      }
#endif
    }
    is2D_ = true;
  }

  // ======= finalize ===========
  template <class Ordinal, class Node>
  void CrsGraph<Ordinal,Node>::finalize(bool OptimizeStorage)
  {
    if (isFinalized () && ! (OptimizeStorage && ! isOptimized ())) {
      // If we've already finalized, and if we don't need to optimize
      // the graph's storage, then we don't have to do anything.
      return;
    }
    if ((indices1D_ == null && indices2D_ == null) || (this->getNumEntries() == 0)) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      if (OptimizeStorage) {
        // move into packed 1D storage
        if (is1DStructure() == false) {
          // allocate 1D storage
          // we these are for host use, so we'll forgo the view
          indices1D_ = arcp<Ordinal>(this->getNumEntries());
        }
        ArrayRCP<size_t> offsets = arcp<size_t>(numRows_+1);
        // copy/pack data
        size_t curoffset = 0;
        size_t curnuminds;
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        newinds = indices1D_.begin();
        for (size_t i=0; i < numRows_; ++i) {
          offsets[i] = curoffset;
          if (is1DStructure()) {
            curnuminds = rowEnds_[i] - rowBegs_[i];
            oldinds = indices1D_.begin() + rowBegs_[i];
          }
          else {
            curnuminds = numEntriesPerRow_[i];
            oldinds = indices2D_[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          newinds += curnuminds;
          curoffset += curnuminds;
        }
        offsets[numRows_] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error,
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        // done with the original row beg/end offsets, can point to the new overlapping one
        rowBegs_   = offsets;
        rowEnds_   = offsets.persistingView(1,numRows_);
        isOpt_     = true;
        is1D_      = true;
        // delete 2D storage (if there was any)
        is2D_      = false;
        numEntriesPerRow_ = null;
        indices2D_        = null;
      }
    }
    isFinalized_ = true;
  }


  // ======= finalize ===========
  // finalize() storage for the graph with associated matrix values
  // this is called from a CrsMatrix, and we're doing the finalize the for the graph and matrix at the same time, so the matrix doesn't have to.
  template <class Ordinal, class Node>
  template <class Scalar>
  void CrsGraph<Ordinal,Node>::finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D)
  {
    if (isFinalized () && ! (OptimizeStorage && ! isOptimized ())) {
      // If we've already finalized, and if we don't need to optimize
      // the graph's storage, then we don't have to do anything.
      return;
    }
    if ((indices1D_.is_null () && indices2D_.is_null ()) ||
        (this->getNumEntries () == 0)) {
      isEmpty_ = true; // The graph is empty.
    }
    else {
      isEmpty_ = false;
      // move into packed 1D storage
      if (OptimizeStorage) {
        if (! is1DStructure ()) {
          // allocate 1D storage
          // we know this is a host-based node, so we'll forgo the view of rowBegs_,rowEnds_
          indices1D_ = arcp<Ordinal>(this->getNumEntries());
          values1D   = arcp<Scalar >(this->getNumEntries());
        }
        ArrayRCP<size_t> offsets = arcp<size_t>(numRows_+1);
        // copy/pack data
        size_t curoffset = 0;
        size_t curnuminds;
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        typename ArrayRCP<Scalar >::iterator oldvals, newvals;
        newinds = indices1D_.begin();
        newvals = values1D.begin();
        for (size_t i=0; i < numRows_; ++i) {
          offsets[i] = curoffset;
          if (is1DStructure()) {
            curnuminds = rowEnds_[i] - rowBegs_[i];
            oldinds = indices1D_.begin() + rowBegs_[i];
            oldvals = values1D.begin() + rowBegs_[i];
          }
          else {
            curnuminds = numEntriesPerRow_[i];
            oldinds = indices2D_[i].begin();
            oldvals = values2D[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          std::copy(oldvals, oldvals+curnuminds, newvals);
          newinds += curnuminds;
          newvals += curnuminds;
          curoffset += curnuminds;
        }
        offsets[numRows_] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error,
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        // done with the original row beg/end offsets, can point to the new overlapping one
        rowBegs_   = offsets;
        rowEnds_   = offsets.persistingView(1,numRows_);
        is1D_      = true;
        isOpt_     = true;
        // delete 2D storage (if there was any)
        is2D_      = false;
        numEntriesPerRow_ = null;
        indices2D_        = null;
        values2D          = null;
      }
    }
    isFinalized_ = true;
  }

} // namespace Kokkos

#endif /* KOKKOS_CRSGRAPH_HPP */
