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
#ifndef __Tpetra_MultiVectorFiller_hpp
#define __Tpetra_MultiVectorFiller_hpp

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <iterator>
#include <set>

namespace Tpetra {
namespace Details {
  // \fn sortAndMergeIn
  // \brief Sort and merge newEntries into allEntries, and make unique.
  //
  // \param allEntries [in/out]: Array in which current entries have
  //   already been stored, and in which the new entries are to be
  //   stored.  Passed by reference in case we need to resize it.
  //
  // \param currentEntries [in/out]: Current entries, which we assume
  //   have already been sorted and made unique.  Aliases the
  //   beginning of allEntries.
  //
  // \param newEntries [in/out] New entries, which have not yet been
  //   sorted or made unique.  This does <i>not</i> alias allEntries.
  //
  // Sort and make entries of newEntries unique.  Resize allEntries if
  // necessary to fit the unique entries of newEntries.  Merge
  // newEntries into allEntries and make the results unique.  (This is
  // cheaper than sorting the whole array.)
  //
  // \return A view of all the entries (current and new) in allEntries.
  template<class T>
  Teuchos::ArrayView<T>
  sortAndMergeIn (Teuchos::Array<T>& allEntries,
                  Teuchos::ArrayView<T> currentEntries,
                  Teuchos::ArrayView<T> newEntries)
  {
    using Teuchos::ArrayView;
    using Teuchos::as;
    typedef typename ArrayView<T>::iterator iter_type;
    typedef typename ArrayView<T>::size_type size_type;

    std::sort (newEntries.begin(), newEntries.end());
    iter_type it = std::unique (newEntries.begin(), newEntries.end());
    const size_type numNew = as<size_type> (it - newEntries.begin());
    // View of the sorted, made-unique new entries to merge in.
    ArrayView<T> newEntriesView = newEntries.view (0, numNew);

    const size_type numCur = currentEntries.size();
    if (allEntries.size() < numCur + numNew) {
      allEntries.resize (numCur + numNew);
    }
    ArrayView<T> allView = allEntries.view (0, numCur + numNew);
    ArrayView<T> newView = allEntries.view (numCur, numNew); // target of copy

    std::copy (newEntries.begin(), newEntries.end(), newView.begin());
    std::inplace_merge (allView.begin(), newView.begin(), allView.end());
    iter_type it2 = std::unique (allView.begin(), allView.end());
    const size_type numTotal = as<size_type> (it2 - allView.begin());

    return allEntries.view (0, numTotal);
  }

  /// \class MultiVectorFillerData
  /// \brief Implementation of fill and local assembly for MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// \param map [in] Map, which this object may or may not use as a
    ///   hint to separate local from nonlocal data.  This need not be
    ///   the same Map as that of the multivector output of
    ///   globalAssemble().
    ///
    /// Before using this object, you should call setNumColumns() to
    /// set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of sumIntoGlobalValues()
    /// won't actually do anything.
    MultiVectorFillerData (const Teuchos::RCP<const map_type>& map) :
      map_ (map),
      numCols_ (0)
    {}

    /// \brief Constructor.
    ///
    /// \param map [in] Map, which this object may or may not use as a
    ///   hint to separate local from nonlocal data.  This need not be
    ///   the same Map as that of the multivector output of
    ///   globalAssemble().
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of sumIntoGlobalValues()
    ///   won't do the right thing.
    MultiVectorFillerData (const Teuchos::RCP<const map_type>& map,
                           const size_t numColumns) :
      map_ (map),
      numCols_ (numColumns),
      sourceIndices_ (numColumns),
      sourceValues_ (numColumns)
    {}

    //! Set the number of columns in the output multivector.
    void
    setNumColumns (const size_t newNumColumns)
    {
      const size_t oldNumColumns = getNumColumns();
      if (newNumColumns >= oldNumColumns) {
        for (size_t j = oldNumColumns; j < newNumColumns; ++j) {
          sourceIndices_.push_back (Teuchos::Array<global_ordinal_type> ());
          sourceValues_.push_back (Teuchos::Array<scalar_type> ());
        }
      }
      else {
        // This may not necessarily deallocate any data, but that's OK.
        sourceIndices_.resize (newNumColumns);
        sourceValues_.resize (newNumColumns);
      }
      numCols_ = oldNumColumns;
    }

    void
    sumIntoGlobalValues (const Teuchos::ArrayView<const global_ordinal_type>& rows,
                         const size_t column,
                         const Teuchos::ArrayView<const scalar_type>& values)
    {
      if (column >= getNumColumns()) {
        for (size_t j = column; j < getNumColumns(); ++j) {
          sourceIndices_.push_back (Teuchos::Array<global_ordinal_type> ());
          sourceValues_.push_back (Teuchos::Array<scalar_type> ());
        }
      }
      std::copy (rows.begin(), rows.end(), std::back_inserter (sourceIndices_[column]));
      std::copy (values.begin(), values.end(), std::back_inserter (sourceValues_[column]));
    }

    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (const Teuchos::ArrayView<const global_ordinal_type>& rows,
                         const Teuchos::ArrayView<const scalar_type>& values)
    {
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef typename Teuchos::ArrayView<const GO>::const_iterator GoIter;
      typedef typename Teuchos::ArrayView<const ST>::const_iterator StIter;

      const size_t numColumns = getNumColumns();
      GoIter rowIter = rows.begin();
      StIter valIter = values.begin();
      for (size_t j = 0; j < numColumns; ++j) {
        GoIter rowIterNext = rowIter + numColumns;
        StIter valIterNext = valIter + numColumns;
        std::copy (rowIter, rowIterNext, std::back_inserter (sourceIndices_[j]));
        std::copy (valIter, valIterNext, std::back_inserter (sourceValues_[j]));
        rowIter = rowIterNext;
        valIter = valIterNext;
      }
    }

    /// \brief Locally assemble into X.
    ///
    /// \param X [in/out] Multivector (overlapping source distribution).
    ///
    /// \param f [in/out] Binary function that defines the combine
    ///   mode.  It must define scalar_type operator (const
    ///   scalar_type&, const scalar_type&).  It need not necessarily
    ///   be commutative or even associative, but it should be
    ///   thread-safe in case we decide to parallelize local assembly.
    ///   We call it via X(i,j) = f(X(i,j), Y(i,j)), so write your
    ///   possibly nonassociative or noncommutative operation
    ///   accordingly.
    ///
    /// X is distributed by the source Map (with possible overlap) of
    /// the Export operation.  The source Map of the Export includes
    /// both the elements owned by this object's constructor's input
    /// Map, and the indices inserted by sumIntoGlobalValues().
    ///
    /// Precondition: The set of global indices in X's Map equals the
    /// the union of the entries of nonlocalIndices_[j] for all valid
    /// columns j.
    ///
    /// \note You can get the usual ADD combine mode by supplying f =
    ///   std::plus<scalar_type>.
    template<class BinaryFunction>
    void
    locallyAssemble (MV& X, BinaryFunction& f)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef node_type NT;

      RCP<const ::Tpetra::Map<LO, GO, NT> > map = X.getMap();
      Array<LO> localIndices;
      const size_t numColumns = getNumColumns();
      for (size_t j = 0; j < numColumns; ++j) {
        const typename Array<const GO>::size_type numIndices = sourceIndices_[j].size();
        // Precompute all the local indices before saving to the
        // vector.  Hopefully this gives us a little bit of locality
        // in the global->local conversion, at the expense of a little
        // more storage.
        if (localIndices.size() < numIndices) {
          localIndices.resize (numIndices);
        }
        ArrayView<LO> localIndicesView = localIndices.view (0, numIndices);
        ArrayView<const GO> globalIndicesView = sourceIndices_[j].view (0, numIndices);
        for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
          localIndices[i] = map->getLocalElement (globalIndicesView[i]);
        }

        ArrayRCP<ST> X_j = X.getDataNonConst (j);
        ArrayView<const ST> localValues = sourceValues_[j].view (0, numIndices);
        for (typename ArrayView<const GO>::size_type i = 0; i < numIndices; ++i) {
          const LO localInd = localIndices[i];
          X_j[localInd] = f (X_j[localInd], localValues[i]);
        }
      }
    }

    //! Overload of locallyAssemble() for the usual ADD combine mode.
    void
    locallyAssemble (MV& X)
    {
      std::plus<scalar_type> f;
      locallyAssemble<std::plus<scalar_type> > (X, f);
    }

    //! Clear the contents of the vector, making it implicitly a vector of zeros.
    void clear() {
      Teuchos::Array<Teuchos::Array<global_ordinal_type> > newSourceIndices;
      Teuchos::Array<Teuchos::Array<scalar_type> > newSourceValues;
      // The standard STL idiom for clearing the contents of a vector completely.
      std::swap (sourceIndices_, newSourceIndices);
      std::swap (sourceValues_, newSourceValues);
    }

    //! All source indices (local and nonlocal) of the source Map, sorted and unique.
    Teuchos::Array<global_ordinal_type>
    getSourceIndices () const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      typedef global_ordinal_type GO;
      typedef typename Array<GO>::size_type size_type;

      Array<GO> allInds (0); // will resize below
      const size_t numCols = getNumColumns();

      if (numCols == 1) {
        // Special case for 1 column avoids copying indices twice.
        // Pick the size of the array exactly the first time so there
        // are at most two allocations (the constructor may choose to
        // allocate).
        const size_type numNew = sourceIndices_[0].size();
        allInds.resize (allInds.size() + numNew);
        std::copy (sourceIndices_[0].begin(), sourceIndices_[0].end(),
                   allInds.begin());
        std::sort (allInds.begin(), allInds.end());
        typename Array<GO>::iterator it =
          std::unique (allInds.begin(), allInds.end());
        const size_type numFinal = as<size_type> (it - allInds.begin());
        allInds.resize (numFinal);
      }
      else {
        // Carefully collect all the row indices one column at a time.
        // This ensures that the total allocation size in this routine
        // is independent of the number of columns.  Also, only sort
        // the current column's indices.  Use merges to ensure sorted
        // order in the collected final result.
        ArrayView<GO> curIndsView = allInds.view (0, 0); // will grow
        Array<GO> newInds;
        for (size_t j = 0; j < numCols; ++j) {
          const size_type numNew = sourceIndices_[j].size();
          if (numNew > newInds.size()) {
            newInds.resize (numNew);
          }
          ArrayView<GO> newIndsView = newInds.view (0, numNew);
          std::copy (sourceIndices_[j].begin(), sourceIndices_[j].end(),
                     newIndsView.begin());
          curIndsView = sortAndMergeIn<GO> (allInds, curIndsView, newIndsView);
        }
      }
      return allInds;
    }

  private:
    Teuchos::RCP<const map_type> map_;
    size_t numCols_;
    Teuchos::Array<Teuchos::Array<global_ordinal_type> > sourceIndices_;
    Teuchos::Array<Teuchos::Array<scalar_type> > sourceValues_;

    size_t getNumColumns() const { return numCols_; }
  };

  /// \class MultiVectorFillerData2
  /// \brief Second implementation of fill and local assembly for
  ///   MultiVectorFiller.
  /// \author Mark Hoemmen
  ///
  /// \tparam MV Specialization of Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFillerData2 : public Teuchos::Describable {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;

    typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Default constructor (sets number of columns to zero).
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///   This need not be the same Map as that of the multivector
    ///   output of globalAssemble().
    ///
    /// Before using this object, you should call setNumColumns() to
    /// set the number of columns in the output multivector.
    /// Otherwise, the two-argument version of sumIntoGlobalValues()
    /// won't actually do anything.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map,
                            const Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_DEFAULT,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out=Teuchos::null) :
      map_ (map),
      numCols_ (0),
      verbLevel_ (verbLevel),
      out_ (out)
    {}

    /// \brief Constructor.
    ///
    /// \param map [in] Map over which to distribute the initial fill.
    ///   This need not be the same Map as that of the multivector
    ///   output of globalAssemble(), but the Map must have the
    ///   same communicator as the multivector output of
    ///   globalAssemble().
    ///
    /// \param numColumns [in] The (expected) number of columns in the
    ///   output multivector.  You can always change this later by
    ///   calling setNumColumns().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, you should
    ///   call setNumColumns() first before inserting any data.
    ///   Otherwise, the two-argument version of sumIntoGlobalValues()
    ///   won't do the right thing.
    MultiVectorFillerData2 (const Teuchos::RCP<const map_type>& map,
                            const size_t numColumns,
                            const Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_DEFAULT,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out=Teuchos::null) :
      map_ (map),
      numCols_ (numColumns),
      localVec_ (new MV (map, numColumns)),
#if 0
      nonlocalIndices_ (numColumns),
      nonlocalValues_ (numColumns),
#endif // 0
      verbLevel_ (verbLevel),
      out_ (out)
    {}

    std::string
    description() const
    {
      std::ostringstream oss;
      oss << "Tpetra::MultiVectorFillerData2<"
          << Teuchos::TypeNameTraits<MV>::name () << ">";
      return oss.str();
    }

    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const
    {
      using std::endl;
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;

      // Set default verbosity if applicable.
      const Teuchos::EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

      RCP<const Teuchos::Comm<int> > comm = map_->getComm();
      const int myImageID = comm->getRank();
      const int numImages = comm->getSize();

      if (vl != VERB_NONE) {
        // Don't set the tab level unless we're printing something.
        Teuchos::OSTab tab (out);

        if (myImageID == 0) { // >= VERB_LOW prints description()
          out << this->description() << endl;
        }
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            if (vl != VERB_LOW) {
              // At verbosity > VERB_LOW, each process prints something.
              out << "Process " << myImageID << ":" << endl;

              Teuchos::OSTab procTab (out);
              // >= VERB_MEDIUM: print the local vector length.
              out << "local length=" << localVec_->getLocalLength();
              if (vl != VERB_MEDIUM) {
                // >= VERB_HIGH: print isConstantStride() and getStride()
                if (localVec_->isConstantStride()) {
                  out << ", constant stride=" << localVec_->getStride() << endl;
                }
                else {
                  out << ", not constant stride" << endl;
                }
                if (vl == VERB_EXTREME) {
                  // VERB_EXTREME: print all the values in the multivector.
                  out << "Local values:" << endl;
                  ArrayRCP<ArrayRCP<const scalar_type> > X = localVec_->get2dView();
                  for (size_t i = 0; i < localVec_->getLocalLength(); ++i) {
                    for (size_t j = 0; j < localVec_->getNumVectors(); ++j) {
                      out << X[j][i];
                      if (j < localVec_->getNumVectors() - 1) {
                        out << " ";
                      }
                    } // for each column
                    out << endl;
                  } // for each row

#if 0
                  out << "Nonlocal indices and values:" << endl;
                  for (size_t j = 0; j < (size_t)nonlocalIndices_.size(); ++j) {
                    ArrayView<const global_ordinal_type> inds = nonlocalIndices_[j]();
                    ArrayView<const scalar_type> vals = nonlocalValues_[j]();
                    for (typename ArrayView<const global_ordinal_type>::size_type k = 0;
                         k < inds.size(); ++k) {
                      out << "X(" << inds[k] << "," << j << ") = " << vals[k] << endl;
                    }
                  }
#endif // 0
                } // if vl == VERB_EXTREME
              } // if (vl != VERB_MEDIUM)
              else { // vl == VERB_LOW
                out << endl;
              }
            } // if vl != VERB_LOW
          } // if it is my process' turn to print
        } // for each process in the communicator
      } // if vl != VERB_NONE
    }

    //! All source indices (local and nonlocal) of the source Map, sorted and unique.
    Teuchos::Array<global_ordinal_type>
    getSourceIndices (const Teuchos::Comm<int>& comm,
                      const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null,
                      const bool debug = false) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Tpetra::global_size_t;
      typedef global_ordinal_type GO;
      const char prefix[] = "Tpetra::MultiVectorFiller::getSourceIndices: ";

      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": getSourceIndices" << std::endl;
        *out << os.str ();
      }

      // Get the nonlocal row indices, sorted and made unique.
      // It's fair to assume that these are not contiguous.
      Array<GO> nonlocalIndices = getSortedUniqueNonlocalIndices (comm, out, debug);

      // Get the local row indices, not necessarily sorted or unique.
      ArrayView<const GO> localIndices = getLocalIndices ();

      // Copy the local indices into the full indices array, and sort
      // them there.  We'll merge in the nonlocal indices below.  This
      // can be more efficient than just sorting all the indices, if
      // there are a lot of nonlocal indices.
      Array<GO> indices (localIndices.size () + nonlocalIndices.size ());
      ArrayView<GO> localIndView = indices.view (0, localIndices.size ());
      TEUCHOS_TEST_FOR_EXCEPTION
        (localIndView.size () > indices.size (), std::logic_error,
         prefix << "localIndView.size() = " << localIndView.size ()
         << " > indices.size() = " << indices.size () << ".");
      TEUCHOS_TEST_FOR_EXCEPTION
        (localIndView.size () != localIndices.size (), std::logic_error,
         prefix << "localIndView.size() = " << localIndView.size ()
         << " != localIndices.size() = " << localIndices.size () << ".");

      std::copy (localIndices.begin (), localIndices.end (), localIndView.begin ());
      std::sort (localIndView.begin (), localIndView.end ());

      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": Right after copying and sorting" << std::endl;
        *out << os.str ();
      }

      // Merge the local and nonlocal indices.
      if (nonlocalIndices.size () > 0) {
        ArrayView<GO> nonlclIndView =
          indices.view (localIndices.size (), nonlocalIndices.size ());

        // We need a view of the output array, because
        // std::inplace_merge needs all its iterator inputs to have
        // the same type.  Debug mode builds are pickier than release
        // mode builds, because the iterators in a debug mode build
        // are of a different type that does run-time checking (they
        // aren't just raw pointers).
        ArrayView<GO> indView = indices.view (0, indices.size ());
        if (debug && ! out.is_null ()) {
          std::ostringstream os;
          os << "Right before std::copy" << std::endl;
          *out << os.str ();
        }
        std::copy (nonlocalIndices.begin (),
                   nonlocalIndices.end (),
                   nonlclIndView.begin ());
        if (debug && ! out.is_null ()) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ": Right before std::inplace_merge"
             << std::endl;
          *out << os.str ();
        }
        std::inplace_merge (indView.begin (),
                            indView.begin () + localIndices.size (),
                            indView.end ());
      }

      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": Done with getSourceIndices"
           << std::endl;
        *out << os.str ();
      }
      return indices;
    }

    /// \brief Set the number of columns in the output multivector.
    ///
    /// \param newNumColumns [in] The new number of columns in the
    ///   multivector.  Zero is allowed; it means "clear local
    ///   storage."
    ///
    /// Setting the number of columns to zero effectively clears out
    /// all local storage, but may not necessarily deallocate nonlocal
    /// storage.  Call clear() to clear out all nonlocal storage.
    void
    setNumColumns (const size_t newNumColumns)
    {
      using Teuchos::Array;
      using Teuchos::Range1D;
      using Teuchos::RCP;

      const size_t oldNumColumns = numCols_;
      if (newNumColumns == oldNumColumns) {
        return; // No side effects if no change.
      }

      RCP<MV> newLclVec;
      if (newNumColumns > oldNumColumns) {
        newLclVec = rcp (new MV (map_, newNumColumns));
        // Assign the contents of the old local multivector to the
        // first oldNumColumns columns of the new local multivector,
        // then get rid of the old local multivector.
        RCP<MV> newLclVecView =
          newLclVec->subViewNonConst (Range1D (0, oldNumColumns-1));
        *newLclVecView = *localVec_;
      }
      else {
        if (newNumColumns == 0) {
          // Tpetra::MultiVector doesn't let you construct a
          // multivector with zero columns.
          newLclVec = Teuchos::null;
        }
        else {
          newLclVec = localVec_->subViewNonConst (Range1D (0, newNumColumns-1));
        }
      }

      // Leave most side effects until the end, for exception safety.
      localVec_ = newLclVec;
      numCols_ = newNumColumns;
    }

    /// \brief Set entry (rows[i],columnIndex) to values[i], for all i.
    ///
    /// \param rows [in] Which rows of the MultiVector to modify.
    /// \param columnIndex [in] Which column of the MultiVector to modify.
    /// \param values [in] The new values.
    void
    sumIntoGlobalValues (const Teuchos::ArrayView<const global_ordinal_type>& rows,
                         const size_t columnIndex,
                         const Teuchos::ArrayView<const scalar_type>& values)
    {
      using Teuchos::ArrayView;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef scalar_type ST;
      typedef decltype (rows.size ()) size_type;
      const char prefix[] = "Tpetra::MultiVectorFiller::sumIntoGlobalValues: ";

      if (map_.is_null ()) {
        return; // the calling process is not participating
      }
      const size_type numEnt = rows.size ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (numEnt != values.size (), std::invalid_argument, prefix
         << "rows.size() = " << numEnt << " != values.size() = "
         << values.size () << ".");

      // FIXME (31 Aug 2015) Don't do this.
      if (columnIndex >= getNumColumns()) {
        // Automatically expand the number of columns.  This
        // implicitly ensures that localVec_ is not null.
        setNumColumns (columnIndex + 1);
      }

      // It would be a lot more efficient for users to get the Kokkos
      // view outside of their fill loop, and fill into it directly.
      auto X_j = localVec_->getVector (columnIndex);
      auto X_j_2d = X_j->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X_j_1d = Kokkos::subview (X_j_2d, Kokkos::ALL (), 0);

      const map_type& theMap = *map_;
      for (size_type k = 0; k < numEnt; ++k) {
        const ST val = values[k];
        const GO gblRowInd = rows[k];
        const LO lclRowInd = theMap.getLocalElement (gblRowInd);
        if (lclRowInd == Teuchos::OrdinalTraits<LO>::invalid ()) {
          // The row doesn't belong to the calling process, so stuff
          // its data in nonlocals_.

          auto& innerMap = nonlocals_[columnIndex];
          auto innerIter = innerMap.find (gblRowInd);
          if (innerIter == innerMap.end ()) {
            innerMap.insert (std::make_pair (gblRowInd, values[k]));
          } else {
            innerIter->second += val;
          }
        }
        else {
          // FIXME (mfh 27 Feb 2012) Allow different combine modes.
          // The current combine mode just adds to the current value
          // at that location.
          X_j_1d[lclRowInd] += val;
        }
      }
    }

    /// \brief Set entry (rows[i],j) to values[i], for all i and j.
    ///
    /// Data for each column are stored contiguously in rows and in
    /// values.  Thus, rows and values are in rowwise order, even
    /// though they may be stored in columnwise order in the
    /// multivector.
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    void
    sumIntoGlobalValues (const Teuchos::ArrayView<const global_ordinal_type>& rows,
                         const Teuchos::ArrayView<const scalar_type>& values)
    {
      using Teuchos::ArrayView;
      typedef typename ArrayView<const global_ordinal_type>::size_type size_type;

      const size_t numCols = getNumColumns();
      for (size_t j = 0; j < numCols; ++j) {
        const size_type offset = numCols*j;
        const size_type len = numCols;
        this->sumIntoGlobalValues (rows.view (offset, len), j, values.view (offset, len));
      }
    }

    /// \brief Locally assemble into X, with user-specified combine mode.
    ///
    /// \param X [in/out] Multivector (overlapping source distribution).
    ///
    /// \param f [in/out] Binary function that defines the combine
    ///   mode.  It must define scalar_type operator (const
    ///   scalar_type&, const scalar_type&).  It need not necessarily
    ///   be commutative or even associative, but it should be
    ///   thread-safe in case we decide to parallelize local assembly.
    ///   We call it via X(i,j) = f(X(i,j), Y(i,j)), so write your
    ///   possibly nonassociative or noncommutative operation
    ///   accordingly.
    ///
    /// X is distributed by the source Map (with possible overlap) of
    /// the Export operation.  The source Map of the Export includes
    /// both the elements owned by this object's constructor's input
    /// Map, and the indices inserted by sumIntoGlobalValues().
    ///
    /// Precondition: The set of global indices in X's Map equals the
    /// union of the global indices in map_ and (the union of the
    /// entries of nonlocalIndices_[j] for all valid columns j).
    ///
    /// \note You can get the usual Tpetra::ADD combine mode by
    ///   supplying f = std::plus<scalar_type>.
    template<class BinaryFunction>
    void
    locallyAssemble (MV& X, BinaryFunction& f)
    {
      typedef scalar_type ST;
      typedef local_ordinal_type LO;
      typedef global_ordinal_type GO;
      typedef typename MV::dual_view_type::t_host::memory_space host_memory_space;

      if (X.getMap ().is_null ()) {
        return; // this process is not participating
      }
      const map_type& srcMap = * (X.getMap ());



      for (size_t j = 0; j < X.getNumVectors (); ++j) {
        auto X_j = X.getVector (j);
        auto X_j_2d = X_j->template getLocalView<host_memory_space> ();
        auto X_j_1d = Kokkos::subview (X_j_2d, Kokkos::ALL (), 0);

        // First, assemble in the local components from localVec_.
        if (! localVec_.is_null () && ! localVec_->getMap ().is_null ()
            && j < localVec_->getNumVectors ()) {
          auto lcl_j = localVec_->getVector (j);
          auto lcl_j_2d = lcl_j->template getLocalView<host_memory_space> ();
          auto lcl_j_1d = Kokkos::subview (lcl_j_2d, Kokkos::ALL (), 0);

          // Iterate over all rows of localVec_.  i_lcl is a local
          // index with respect to localVec_'s Map, which may differ
          // from X's Map.
          const map_type& lclMap = * (localVec_->getMap ());
          const LO lclNumRows = static_cast<LO> (lcl_j->getLocalLength ());
          for (LO i_lcl = 0; i_lcl < lclNumRows; ++i_lcl) {
            const GO i_gbl = lclMap.getGlobalElement (i_lcl);
            const LO i_X = srcMap.getLocalElement (i_gbl);
            X_j_1d(i_X) = f (X_j_1d(i_X), lcl_j_1d(i_lcl));
          }
        }

        // Second, assemble in the nonlocal components.
        auto outerIter = nonlocals_.find (j);
        if (outerIter != nonlocals_.end ()) {
          auto beg = outerIter->second.begin ();
          auto end = outerIter->second.end ();
          for (auto innerIter = beg; innerIter != end; ++innerIter) {
            const GO gblRowInd = innerIter->first;
            const LO lclRowInd = srcMap.getLocalElement (gblRowInd);

            if (lclRowInd != Teuchos::OrdinalTraits<LO>::invalid ()) {
              const ST val = innerIter->second;
              X_j_1d(lclRowInd) = f (X_j_1d(lclRowInd), val);
            }
          }
        }
      }
    }

    /// \brief Do the "local" (not MPI) part of assembly
    ///
    /// Call this method for the usual Tpetra::ADD combine mode.  If
    /// you want a different combine mode, call the overload that
    /// takes two arguments, one of which is a binary function.
    void
    locallyAssemble (MV& X)
    {
      std::plus<scalar_type> f;
      locallyAssemble<std::plus<scalar_type> > (X, f);
    }

    /// \brief Clear the contents of the multivector.
    ///
    /// This fills the vector with zeros, and also removes nonlocal
    /// data.  It does <i>not</i> deallocate all storage.  For that,
    /// you need to set the number of columns to zero.
    void clear ()
    {
      std::map<size_t, std::map<global_ordinal_type, scalar_type> > newNonlcls;
      // The standard STL idiom for clearing the contents of a
      // container completely.  Setting the size to zero may not
      // necessarily deallocate data.
      std::swap (nonlocals_, newNonlcls);

      // Don't actually deallocate the multivector of local entries.
      // Just fill it with zero.  This is because the caller hasn't
      // reset the number of columns.
      if (! localVec_.is_null ()) {
        localVec_->putScalar (Teuchos::ScalarTraits<scalar_type>::zero());
      }
    }

  private:
    //! Map that tells us which entries are locally owned.
    Teuchos::RCP<const map_type> map_;

    //! The number of columns in the (output) multivector.
    size_t numCols_;

    //! Multivector of locally stored (i.e., owned by \c map_) entries.
    Teuchos::RCP<MV> localVec_;

    //! Map from column index to (map from row index to value).
    std::map<size_t, std::map<global_ordinal_type, scalar_type> > nonlocals_;

    //! Verbosity level of this object (mainly used for debugging output).
    Teuchos::EVerbosityLevel verbLevel_;

    //! Output stream (mainly used for debugging output).
    Teuchos::RCP<Teuchos::FancyOStream> out_;

    //! The number of columns in the (output) multivector.
    size_t getNumColumns() const { return numCols_; }

    //! The locally owned row indices, not necessarily sorted or unique.
    Teuchos::ArrayView<const global_ordinal_type>
    getLocalIndices() const
    {
      return map_->getNodeElementList ();
    }

    //! The inserted nonlocal row indices, sorted and made unique.
    Teuchos::Array<global_ordinal_type>
    getSortedUniqueNonlocalIndices (const Teuchos::Comm<int>& comm,
                                    const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null,
                                    const bool debug = false) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::as;
      typedef global_ordinal_type GO;

      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": getSortedUniqueNonlocalIndices"
           << std::endl;
        *out << os.str ();
      }

      // FIXME (mfh 30 Aug 2015) The std::set algorithm is harder to
      // thread-parallelize, but it should be easier to make the new
      // test pass with this.
      std::set<GO> indsOut;
      const size_t numCols = getNumColumns ();
      for (size_t j = 0; j < numCols; ++j) {
        auto outerIter = nonlocals_.find (j);
        if (outerIter != nonlocals_.end ()) {
          auto beg = outerIter->second.begin ();
          auto end = outerIter->second.end ();
          for (auto innerIter = beg; innerIter != end; ++innerIter) {
            // *innerIter: (global row index, value) pair.
            indsOut.insert (innerIter->first);
          }
        }
      }

      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": Made nonlocals set" << std::endl;
        *out << os.str ();
      }

      Array<GO> allNonlocals (indsOut.begin (), indsOut.end ());
      if (debug && ! out.is_null ()) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ": Done with "
          "getSortedUniqueNonlocalIndices" << std::endl;
        *out << os.str ();
      }
      return allNonlocals;
    }
  };
} // namespace Details
} // namespace Tpetra

namespace Tpetra {

  /// \class MultiVectorFiller
  /// \brief Adds nonlocal sum-into functionality to Tpetra::MultiVector.
  /// \author Mark Hoemmen
  ///
  /// MultiVector does not allow modifications to rows not owned by
  /// the calling process.  This makes the implementation of
  /// MultiVector simple, fast, and local, but it complicates
  /// operations like assembly of the right-hand side in the finite
  /// element method.  CrsMatrix already supports insertion into
  /// nonowned rows, so that makes assembly not consistent for sparse
  /// matrices and vectors.  MultiVectorFiller remedies this by
  /// specializing for the case of filling a (multi)vector.  It
  /// provides a constructor and sumIntoGlobalValues() methods that
  /// look just like those of MultiVector, so you can swap it in at
  /// the assembly stage.  Then, you can call globalAssemble() to turn
  /// the MultiVectorFiller into a standard MultiVector.
  ///
  /// The globalAssemble() method will redistribute data if necessary.
  /// You can force it to use the Map provided to the constructor (or
  /// the Map provided to the last call to globalAssemble(), which is
  /// a common use case for doing nonlinear solves where the sparse
  /// matrix structure doesn't change).
  ///
  /// \tparam MV Specialization of Tpetra::MultiVector.
  template<class MV>
  class MultiVectorFiller {
  public:
    typedef typename MV::scalar_type scalar_type;
    typedef typename MV::local_ordinal_type local_ordinal_type;
    typedef typename MV::global_ordinal_type global_ordinal_type;
    typedef typename MV::node_type node_type;
    typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    /// \brief Constructor.
    ///
    /// The constructor takes the same arguments as MV's constructor.
    /// We did this on purpose, so that you can construct
    /// MultiVectorFiller as you would a standard MV.
    ///
    /// \param map [in] A Map with the same communicator and Kokkos
    ///   Node as the output multivector of globalAssemble().  This
    ///   need not be the same Map as that output multivector's Map.
    ///
    /// \param numCols [in] Expected number of columns in the output
    ///   multivector of globalAssemble().
    ///
    /// \note If the number of columns given here is not the same as
    ///   the number of columns in the output multivector, the
    ///   two-argument version of sumIntoGlobalValues() won't do
    ///   the right thing.  Use the three-argument version of that
    ///   method if you don't know how many columns there will be in
    ///   advance.
    ///
    /// \note Not providing the output multivector in the constructor
    ///   gives globalAssemble() more flexibility.  For example, its
    ///   output multivector may have any distribution with the same
    ///   global number of rows.  Furthermore, decoupling insertion
    ///   from the output multivector lets this object store
    ///   preassembled data in whatever format it likes.  It doesn't
    ///   have to insert directly into the output multivector until
    ///   assembly time.  (This may improve performance by amortizing
    ///   the cost of global->local index conversions, for example.)
    MultiVectorFiller (const Teuchos::RCP<const map_type>& map,
                       const size_t numCols);

    /// \brief Assemble all the data (local and nonlocal) into X_out.
    ///
    /// You can call this method multiple times with different
    /// multivector output arguments.  If those arguments have the
    /// same Map, this method will attempt to reuse the Export object
    /// each time.  It will only reuse if the new target Map is the
    /// same as (in the sense of Map::isSameAs()) the previous target
    /// Map, unless you force reuse with the second argument (that
    /// saves a few global reductions for the check).
    ///
    /// \param X_out [in/out] MultiVector with a nonoverlapping Map.
    ///   That Map need not be the same as the one with which this
    ///   object was created.
    ///
    /// \param forceReuseMap [in] If true, assume that X_out has the
    ///   same Map (in the sense of Map::isSameAs()) as the target Map
    ///   in the previous call to this method.  If this method was not
    ///   called before, then assume that X_out has the same Map as
    ///   the argument to the constructor of this object.
    void
    globalAssemble (MV& X_out, const bool forceReuseMap = false);

    //! Clear contents, in preparation for another fill cycle.
    void clear () {
      data_.clear ();
    }

    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      data_.describe (out, verbLevel);
    }

    /// \brief Sum data into the multivector.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param column [in] Index of the column in which to insert.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows,
                         size_t column,
                         Teuchos::ArrayView<const scalar_type> values)
    {
      data_.sumIntoGlobalValues (rows, column, values);
    }

    /// \brief Sum data into the multivector.
    ///
    /// In rows and values, the data for each column are stored
    /// contiguously.  Thus, each array is really a matrix in rowwise
    /// order.  (The multivector may use columnwise order internally.)
    ///
    /// Be sure that the number of columns is set correctly before
    /// calling this.
    ///
    /// \param rows [in] Array of global rows for which to insert
    ///   values.  Must have the same length as values.
    ///
    /// \param values [in] Array of values to insert.  Must have the
    ///   same length as rows.  rows[i] is the row in which values[i]
    ///   is to be inserted.
    void
    sumIntoGlobalValues (Teuchos::ArrayView<const global_ordinal_type> rows,
                         Teuchos::ArrayView<const scalar_type> values)
    {
      data_.sumIntoGlobalValues (rows, values);
    }

  private:
    //! Map with which this object was created ("ctor" == "constructor").
    Teuchos::RCP<const map_type> ctorMap_;

    /// \brief Source Map of the last call to globalAssemble().
    ///
    /// A possibly overlapping Map that describes the distribution of
    /// input data, and is the source of the Export.
    Teuchos::RCP<const map_type> sourceMap_;

    /// \brief The target Map of the last call to globalAssemble().
    ///
    /// A nonoverlapping Map which is the target of the Export, and
    /// describes the distribution of the output multivector of
    /// globalAssemble().
    Teuchos::RCP<const map_type> targetMap_;

    /// \brief The source MV of the Export operation in globalAssemble().
    ///
    /// The globalAssemble() method uses this multivector as the
    /// source of the Export operation.  The Export redistributes the
    /// data from a possibly overlapping distribution (reflecting how
    /// elements were inserted) to a nonoverlapping distribution (that
    /// of the output multivector \c X_out of globalAssemble()).  This
    /// is the simplest way to implement redistribution.
    ///
    /// We avoid resizing and reallocating \c sourceVec_ by using a
    /// contiguous subview as the source of the Export, if \c X_out
    /// has fewer columns than \c sourceVec_.
    Teuchos::RCP<MV> sourceVec_;

    /// \brief Storage for inserted indices and values.
    ///
    /// We separate this to facilitate experimentation with different
    /// storage formats.
    Tpetra::Details::MultiVectorFillerData2<MV> data_;

    typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, node_type> export_type;
    Teuchos::RCP<export_type> exporter_;

    /// \brief Assemble the local data into \c X_in.
    ///
    /// This method is called by globalAssemble(), in which \c X_in
    /// is the multivector with the (possibly overlapping) source
    /// distribution.
    void locallyAssemble (MV& X_in) {
      data_.locallyAssemble (X_in);
    }

    //! All source indices (local and nonlocal) of the source Map.
    Teuchos::Array<global_ordinal_type>
    getSourceIndices (const Teuchos::Comm<int>& comm,
                      const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null,
                      const bool debug = false) const
    {
      return data_.getSourceIndices (comm, out, debug);
    }

    /// \brief Compute the source Map (for the Export) from the source indices.
    ///
    /// indexBase, comm, and node are the same as the input arguments
    /// of the Map constructors.
    ///
    /// \note This is a collective operation.
    ///
    /// \return The (possibly overlapping) source Map.
    Teuchos::RCP<const map_type>
    computeSourceMap (const global_ordinal_type indexBase,
                      const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                      const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null,
                      const bool debug = false);
  };

  template<class MV>
  MultiVectorFiller<MV>::MultiVectorFiller (const Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>& map, const size_t numCols)
    : ctorMap_ (map),
      sourceMap_ (Teuchos::null),
      targetMap_ (Teuchos::null),
      data_ (map, numCols),
      exporter_ (Teuchos::null)
  {}

  template<class MV>
  Teuchos::RCP<const typename MultiVectorFiller<MV>::map_type>
  MultiVectorFiller<MV>::
  computeSourceMap (const global_ordinal_type indexBase,
                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const bool debug)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::rcp;
    typedef global_ordinal_type GO;

    if (debug && ! out.is_null ()) {
      out->pushTab ();
      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Proc " << myRank << ": computeSourceMap" << std::endl;
      *out << os.str ();
    }

    Array<GO> indices = getSourceIndices (*comm, out, debug);

    if (debug && ! out.is_null ()) {
      std::ostringstream os;
      const int myRank = comm->getRank ();
      os << "Proc " << myRank << ": computeSourceMap: About to create Map"
         << std::endl;
      *out << os.str ();
      out->popTab ();
    }

    // Passing "invalid" for the numGlobalElements argument of the Map
    // constructor tells the Map to compute the global number of
    // elements itself.
    typedef Tpetra::global_size_t GST;
    const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
    return rcp (new map_type (INV, indices, indexBase, comm));
  }

  template<class MV>
  void
  MultiVectorFiller<MV>::
  globalAssemble (MV& X_out, const bool forceReuseMap)
  {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::rcp;
    RCP<Teuchos::FancyOStream> outPtr;
    const bool debug = false;

    if (debug) {
      outPtr = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
      outPtr->pushTab ();
    }

    const size_t numVecs = X_out.getNumVectors ();
    if (numVecs == 0) {
      // Nothing to do!  Of course, this does not check for whether
      // X_out has the right number of rows.  That's OK, though.
      // Compare to the definition of the BLAS' _AXPY for an input
      // vector containing NaNs, but multiplied by alpha=0.
      return;
    }
    //
    // Get the target Map of the Export.  If X_out's Map is the same
    // as the target Map from last time, then we may be able to
    // recycle the Export from last time, if the source Map is also
    // the same.
    //
    RCP<const map_type> targetMap;
    bool assumeSameTargetMap = false;
    if (targetMap_.is_null()) {
      targetMap_ = X_out.getMap();
      targetMap = targetMap_;
      assumeSameTargetMap = false;
    }
    else {
      if (forceReuseMap) {
        targetMap = targetMap_;
        assumeSameTargetMap = true;
      }
      else  {
        // If X_out's Map is the same as targetMap_, we may be able to
        // reuse the Export.  Constructing the Export may be more
        // expensive than calling isSameAs() (which requires just a
        // few reductions and reading through the lists of owned
        // global indices), so it's worth checking.
        if (targetMap_->isSameAs (*(X_out.getMap()))) {
          assumeSameTargetMap = true;
          targetMap = targetMap_;
        }
      }
    }

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right before getting source Map" << std::endl;
      *outPtr << os.str ();
    }

    //
    // Get the source Map of the Export.  If the source Map of the
    // Export is the same as last time, then we may be able to recycle
    // the Export from last time, if the target Map is also the same.
    //
    RCP<const map_type> sourceMap;
    bool computedSourceMap = false;
    {
      if (forceReuseMap && ! sourceMap_.is_null()) {
        sourceMap = sourceMap_;
      }
      else {
        sourceMap = computeSourceMap (ctorMap_->getIndexBase (),
                                      ctorMap_->getComm (),
                                      outPtr, debug);
        computedSourceMap = true;
      }
    }
    if (computedSourceMap) {
      sourceMap_ = sourceMap;
    }

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right before computing Export" << std::endl;
      *outPtr << os.str ();
    }

    //
    // Now that we have the source and target Maps of the Export, we
    // can check whether we can recycle the Export from last time.
    //
    const bool mustComputeExport =
      (exporter_.is_null() || (assumeSameTargetMap && ! computedSourceMap));
    if (mustComputeExport) {
      exporter_ = rcp (new export_type (sourceMap_, targetMap_));
    }

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right after computing Export" << std::endl;
      *outPtr << os.str ();
    }

    // Source multivector for the Export.
    RCP<MV> X_in;
    const bool mustReallocateVec = sourceVec_.is_null() ||
      sourceVec_->getNumVectors() < numVecs || ! assumeSameTargetMap;

    if (mustReallocateVec) {
      // Reallocating nonlocalVec_ ensures that it has the right Map.
      sourceVec_ = rcp (new MV (sourceMap_, numVecs));
      X_in = sourceVec_;
    } else {
      if (sourceVec_->getNumVectors() == numVecs) {
        X_in = sourceVec_;
      } else { // sourceVec_ has more vectors than needed.
        X_in = sourceVec_->subViewNonConst (Range1D (0, numVecs-1));
      }
    }

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right before locallyAssemble" << std::endl;
      *outPtr << os.str ();
    }

    // "Locally assemble" the data into X_in by summing together
    // entries with the same indices.
    locallyAssemble (*X_in);

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right after locallyAssemble" << std::endl;
      *outPtr << os.str ();
    }

    // Do the Export.
    const Tpetra::CombineMode combineMode = Tpetra::ADD;
    X_out.doExport (*X_in, *exporter_, combineMode);

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Right after Export" << std::endl;
      *outPtr << os.str ();
    }

    // Clean out the locally owned data for next time.
    X_in->putScalar (Teuchos::ScalarTraits<scalar_type>::zero ());

    // FIXME (mfh 27 Aug 2015) Clean out the remote data for next
    // time.  See Bug 6398.

    if (debug && ! outPtr.is_null ()) {
      std::ostringstream os;
      const Teuchos::Comm<int>& comm = * (X_out.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      os << "Proc " << myRank << ": Done with globalAssemble" << std::endl;
      *outPtr << os.str ();
      outPtr->pushTab ();
    }
  }

  namespace Test {

    /// \class MultiVectorFillerTester
    /// \brief Tests for MultiVectorFiller
    /// \author Mark Hoemmen
    ///
    /// \tparam MV A specialization of Tpetra::MultiVector.
    template<class MV>
    class MultiVectorFillerTester {
    public:
      typedef typename MV::scalar_type scalar_type;
      typedef typename MV::local_ordinal_type local_ordinal_type;
      typedef typename MV::global_ordinal_type global_ordinal_type;
      typedef typename MV::node_type node_type;
      typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

      /// \brief Test global assembly when constructor Map = target Map.
      ///
      /// Constructor Map = target Map is a common case for finite
      /// element assembly.  This method current only tests the
      /// version of MultiVectorFiller::sumIntoGlobalValues() that
      /// works on one column at a time.
      ///
      /// If any test fails, this method throws an exception.
      static void
      testSameMap (const Teuchos::RCP<const map_type>& targetMap,
                   const global_ordinal_type eltSize, // Must be odd
                   const size_t numCols,
                   const Teuchos::RCP<Teuchos::FancyOStream>& outStream=Teuchos::null,
                   const Teuchos::EVerbosityLevel verbosityLevel=Teuchos::VERB_DEFAULT)
      {
        using Teuchos::Array;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::as;
        using Teuchos::Comm;
        using Teuchos::FancyOStream;
        using Teuchos::getFancyOStream;
        using Teuchos::oblackholestream;
        using Teuchos::ptr;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcpFromRef;
        using Teuchos::REDUCE_SUM;
        using Teuchos::reduceAll;
        using std::cerr;
        using std::endl;

        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef scalar_type ST;
        typedef Teuchos::ScalarTraits<ST> STS;

        TEUCHOS_TEST_FOR_EXCEPTION(eltSize % 2 == 0, std::invalid_argument,
          "Element size (eltSize) argument must be odd.");
        TEUCHOS_TEST_FOR_EXCEPTION(numCols == 0, std::invalid_argument,
          "Number of columns (numCols) argument must be nonzero.");
        // Default behavior is to print nothing out.
        RCP<FancyOStream> out = outStream.is_null() ?
          getFancyOStream (rcp (new oblackholestream)) : outStream;
        const Teuchos::EVerbosityLevel verbLevel =
          (verbosityLevel == Teuchos::VERB_DEFAULT) ?
          Teuchos::VERB_NONE : verbosityLevel;

        //RCP<MV> X = rcp (new MV (targetMap, numCols));

        const GO indexBase = targetMap->getIndexBase();
        Array<GO> rows (eltSize);
        Array<ST> values (eltSize);
        std::fill (values.begin(), values.end(), STS::one());

        // Make this a pointer so we can free its contents, in case
        // those contents depend on the input to globalAssemble().
        RCP<MultiVectorFiller<MV> > filler =
          rcp (new MultiVectorFiller<MV> (targetMap, numCols));

        TEUCHOS_TEST_FOR_EXCEPTION(! targetMap->isContiguous(),
          std::invalid_argument, "MultiVectorFiller test currently only works "
          "for contiguous Maps.");

        const GO minGlobalIndex = targetMap->getMinGlobalIndex();
        const GO maxGlobalIndex = targetMap->getMaxGlobalIndex();
        const GO minAllGlobalIndex = targetMap->getMinAllGlobalIndex();
        const GO maxAllGlobalIndex = targetMap->getMaxAllGlobalIndex();
        for (size_t j = 0; j < numCols; ++j) {
          for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
            // Overlap over processes, without running out of bounds.
            const GO start = std::max (i - eltSize/2, minAllGlobalIndex);
            const GO end = std::min (i + eltSize/2, maxAllGlobalIndex);
            const GO len = end - start + 1;

            TEUCHOS_TEST_FOR_EXCEPTION(len > eltSize, std::logic_error,
              "At start,end = " << start << "," << end << ", len = " << len
              << " > eltSize = " << eltSize << ".");

            for (GO k = 0; k < len; ++k) {
              rows[k] = start + k;
            }
            if (verbLevel == Teuchos::VERB_EXTREME) {
              *out << "Inserting: "
                   << Teuchos::toString (rows.view(0,len)) << ", "
                   << Teuchos::toString (values.view(0, len)) << std::endl;
            }
            filler->sumIntoGlobalValues (rows.view(0, len), j, values.view(0, len));
          }
        }

        if (verbLevel == Teuchos::VERB_EXTREME) {
          *out << "Filler:" << std::endl;
          filler->describe (*out, verbLevel);
          *out << std::endl;
        }

        MV X_out (targetMap, numCols);
        filler->globalAssemble (X_out);
        filler = Teuchos::null;

        if (verbLevel == Teuchos::VERB_EXTREME) {
          *out << "X_out:" << std::endl;
          X_out.describe (*out, verbLevel);
          *out << std::endl;
        }

        // Create multivector for comparison.
        MV X_expected (targetMap, numCols);
        const scalar_type two = STS::one() + STS::one();
        for (size_t j = 0; j < numCols; ++j) {
          ArrayRCP<ST> X_j = X_expected.getDataNonConst (j);

          // Each entry of the column should have the value eltSize,
          // except for the first and last few entries in the whole
          // column (globally, not locally).
          for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
            const LO localIndex = targetMap->getLocalElement (i);
            TEUCHOS_TEST_FOR_EXCEPTION(i == Teuchos::OrdinalTraits<LO>::invalid(),
              std::logic_error, "Global index " << i << " is not in the "
              "multivector's Map.");

            if (i <= minAllGlobalIndex + eltSize/2) {
              X_j[localIndex] = two + as<ST>(i) - as<ST>(indexBase);
            }
            else if (i >= maxAllGlobalIndex - eltSize/2) {
              X_j[localIndex] = two + as<ST>(maxAllGlobalIndex) - as<ST>(i);
            }
            else {
              X_j[localIndex] = as<ST>(eltSize);
            }
          } // for each global index which my process owns
        } // for each column of the multivector

        if (verbLevel == Teuchos::VERB_EXTREME) {
          *out << "X_expected:" << std::endl;
          X_expected.describe (*out, verbLevel);
          *out << std::endl;
        }

        Array<GO> errorLocations;
        for (size_t j = 0; j < numCols; ++j) {
          ArrayRCP<const ST> X_out_j = X_out.getData (j);
          ArrayRCP<const ST> X_expected_j = X_expected.getData (j);

          // Each entry of the column should have the value eltSize,
          // except for the first and last few entries in the whole
          // column (globally, not locally).
          for (GO i = minGlobalIndex; i <= maxGlobalIndex; ++i) {
            const LO localIndex = targetMap->getLocalElement (i);
            TEUCHOS_TEST_FOR_EXCEPTION(i == Teuchos::OrdinalTraits<LO>::invalid(),
              std::logic_error, "Global index " << i << " is not in the "
              "multivector's Map.");

            // The floating-point additions should be exact in this
            // case, except for very large values of eltSize.
            if (X_out_j[localIndex] != X_expected_j[localIndex]) {
              errorLocations.push_back (i);
            }
          } // for each global index which my process owns

          const typename Array<GO>::size_type localNumErrors = errorLocations.size();
          typename Array<GO>::size_type globalNumErrors = 0;
          RCP<const Comm<int> > comm = targetMap->getComm();
          reduceAll (*comm, REDUCE_SUM, localNumErrors, ptr (&globalNumErrors));

          if (globalNumErrors != 0) {
            std::ostringstream os;
            os << "Proc " << comm->getRank() << ": " << localNumErrors
               << " incorrect value" << (localNumErrors != 1 ? "s" : "")
               << ".  Error locations: [ ";
            std::copy (errorLocations.begin(), errorLocations.end(),
                       std::ostream_iterator<GO> (os, " "));
            os << "].";
            // Iterate through all processes in the communicator,
            // printing out each process' local errors.
            for (int p = 0; p < comm->getSize(); ++p) {
              if (p == comm->getRank()) {
                cerr << os.str() << endl;
              }
              // Barriers to let output finish.
              comm->barrier();
              comm->barrier();
              comm->barrier();
            }
            TEUCHOS_TEST_FOR_EXCEPTION(globalNumErrors != 0, std::logic_error,
              "Over all procs: " << globalNumErrors << " total error"
              << (globalNumErrors != 1 ? "s" : "") << ".");
          } // if there were any errors in column j
        } // for each column j
      }
    };

    //! Instantiate a MultiVectorFillerTester and run the test.
    template<class ScalarType,
             class LocalOrdinalType,
             class GlobalOrdinalType,
             class NodeType>
    void
    testMultiVectorFiller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                           const size_t unknownsPerNode,
                           const GlobalOrdinalType unknownsPerElt,
                           const size_t numCols,
                           const Teuchos::RCP<Teuchos::FancyOStream>& outStream,
                           const Teuchos::EVerbosityLevel verbLevel)
    {
      using Tpetra::createContigMapWithNode;
      using Teuchos::FancyOStream;
      using Teuchos::getFancyOStream;
      using Teuchos::oblackholestream;
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using std::cerr;
      using std::endl;

      typedef ScalarType ST;
      typedef LocalOrdinalType LO;
      typedef GlobalOrdinalType GO;
      typedef NodeType NT;
      typedef ::Tpetra::Map<LO, GO, NT> map_type;
      typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;
      typedef Tpetra::global_size_t GST;

      RCP<FancyOStream> out = outStream.is_null() ?
        getFancyOStream (rcp (new oblackholestream)) : outStream;
      const GST INV = Teuchos::OrdinalTraits<GST>::invalid ();
      const GO indexBase = 0;
      RCP<const map_type> targetMap =
        rcp (new map_type (INV, unknownsPerNode, indexBase, comm));

      std::ostringstream os;
      try {
        MultiVectorFillerTester<MV>::testSameMap (targetMap, unknownsPerElt, numCols, out, verbLevel);
      } catch (std::exception& e) {
        os << e.what();
      }

      for (int p = 0; p < comm->getSize(); ++p) {
        if (p == comm->getRank()) {
          *out << "On process " << comm->getRank() << ": " << os.str() << endl;
        }
        comm->barrier();
        comm->barrier();
        comm->barrier();
      }
    }

    /// \brief Test the sortAndMergeIn() utility function.
    ///
    /// If any test fails, this function throws std::logic_error with
    /// an informative message.  If all tests pass, this function
    /// returns with no externally visible side effects.
    void
    testSortAndMergeIn ();

  } // namespace Test
} // namespace Tpetra


#endif // __Tpetra_MultiVectorFiller_hpp
