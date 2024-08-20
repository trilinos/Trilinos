// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_COOMATRIX_HPP
#define TPETRA_DETAILS_COOMATRIX_HPP

/// \file Tpetra_Details_CooMatrix.hpp
/// \brief Declaration and definition of Tpetra::Details::CooMatrix,
///   an implementation detail of Tpetra::CrsMatrix (sparse matrix)
///   file input and output.

#include "TpetraCore_config.h"
#include "Tpetra_Details_PackTriples.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>


namespace Tpetra {
namespace Details {

// Implementation details of Tpetra::Details.
// So, users REALLY should not use anything in here.
namespace Impl {

/// \brief Type of each (row index, column index) pair in the
///   Tpetra::Details::CooMatrix (see below).
template<class IndexType>
struct CooGraphEntry {
  IndexType row;
  IndexType col;
};

/// \brief Function comparing two CooGraphEntry structs,
///   lexicographically, first by row index, then by column index.
///
/// CooMatrix (see below) will use this.
template<class IndexType>
struct CompareCooGraphEntries {
  bool
  operator () (const CooGraphEntry<IndexType>& x,
               const CooGraphEntry<IndexType>& y) const
  {
    return (x.row < y.row) || (x.row == y.row && x.col < y.col);
  }
};

/// \brief Implementation detail of Tpetra::Details::CooMatrix (which
///   see below).
template<class SC, class GO>
class CooMatrixImpl {
private:
  /// \brief Type of the set of this matrix's locally stored entries.
  ///
  /// The set sorts its entries lexicographically, first by row index,
  /// then by column index.  It uses global indices, because it's
  /// easier to do Matrix Market input and output that way.  (How else
  /// could you process globally indexed matrix triples when you don't
  /// have Maps yet?)
  typedef std::map<CooGraphEntry<GO>, SC,
                   CompareCooGraphEntries<GO> > entries_type;

  /// \brief Set of this matrix's locally stored entries, sorted
  ///   lexicographically, first by row index, then by column index.
  entries_type entries_;

public:
  //! Type for packing and unpacking data.
  typedef char packet_type;

  /// \brief Default constructor.
  CooMatrixImpl () = default;

  /// \brief Insert one entry locally into the sparse matrix, if it
  ///   does not exist there yet.  If it does exist, sum the values.
  ///
  /// \param gblRowInd [in] Global row index of the entry to insert.
  /// \param gblColInd [in] Global column index of the entry to insert.
  /// \param val [in] Value of the matrix entry to insert / sum.
  void
  sumIntoGlobalValue (const GO gblRowInd,
                      const GO gblColInd,
                      const SC& val)
  {
    // There's no sense in worrying about the insertion hint, since
    // indices may be all over the place.  Users make no promises of
    // sorting or locality of input.
    CooGraphEntry<GO> ij {gblRowInd, gblColInd};
    auto result = this->entries_.insert ({ij, val});
    if (! result.second) { // already in the map
      result.first->second += val; // sum in the new value
    }
  }

  /// \brief Insert multiple entries locally into the sparse matrix.
  ///
  /// This works like multiple calls to sumIntoGlobalValue.
  ///
  /// \param gblRowInd [in] Global row indices of the entries to insert.
  /// \param gblColInd [in] Global column indices of the entries to insert.
  /// \param val [in] Values of the matrix entries to insert / sum.
  /// \param numEnt [in] Number of entries to insert.
  void
  sumIntoGlobalValues (const GO gblRowInds[],
                       const GO gblColInds[],
                       const SC vals[],
                       const std::size_t numEnt)
  {
    for (std::size_t k = 0; k < numEnt; ++k) {
      // There's no sense in worrying about the insertion hint, since
      // indices may be all over the place.  Users make no promises of
      // sorting or locality of input.
      CooGraphEntry<GO> ij {gblRowInds[k], gblColInds[k]};
      const SC val = vals[k];
      auto result = this->entries_.insert ({ij, val});
      if (! result.second) { // already in the map
        result.first->second += val; // sum in the new value
      }
    }
  }

  //! Number of entries in the sparse matrix on the calling process.
  std::size_t
  getLclNumEntries () const
  {
    return this->entries_.size ();
  }

  /// \brief Execute the given function for all entries of the sparse
  ///   matrix, sequentially (no thread parallelism).  Do not modify
  ///   entries.
  void
  forAllEntries (std::function<void (const GO, const GO, const SC&)> f) const
  {
    for (auto iter = this->entries_.begin ();
         iter != this->entries_.end (); ++iter) {
      f (iter->first.row, iter->first.col, iter->second);
    }
  }

  /// \brief Into global row tgtGblRow of \c *this, merge global row
  ///   srcGblRow of \c src.
  ///
  /// For matrix entries with the same row <i>and</i> column indices,
  /// sum the values into the entry in tgtEntries.
  void
  mergeIntoRow (const GO tgtGblRow,
                const CooMatrixImpl<SC, GO>& src,
                const GO srcGblRow)
  {
    // const char prefix[] =
    //   "Tpetra::Details::Impl::CooMatrixImpl::mergeIntoRow: ";

    entries_type& tgtEntries = this->entries_;
    const entries_type& srcEntries = src.entries_;

    // Find all entries with the given global row index.  The min GO
    // value is guaranteed to be the least possible column index, so
    // beg1 is a lower bound for all (row index, column index) pairs.
    // Lower bound is inclusive; upper bound is exclusive.
    auto srcBeg = srcEntries.lower_bound ({srcGblRow, std::numeric_limits<GO>::min ()});
    auto srcEnd = srcEntries.upper_bound ({srcGblRow+1, std::numeric_limits<GO>::min ()});
    auto tgtBeg = tgtEntries.lower_bound ({tgtGblRow, std::numeric_limits<GO>::min ()});
    auto tgtEnd = tgtEntries.upper_bound ({tgtGblRow+1, std::numeric_limits<GO>::min ()});

    // Don't waste time iterating over the current row of *this, if
    // the current row of src is empty.
    if (srcBeg != srcEnd) {
      for (auto tgtCur = tgtBeg; tgtCur != tgtEnd; ++tgtCur) {
        auto srcCur = srcBeg;
        while (srcCur != srcEnd && srcCur->first.col < tgtCur->first.col) {
          ++srcCur;
        }
        // At this point, one of the following is true:
        //
        // 1. srcCur == srcEnd, thus nothing more to insert.
        // 2. srcCur->first.col > tgtCur->first.col, thus the current
        //    row of src has no entry matching tgtCur->first, so we
        //    have to insert it.  Insertion does not invalidate
        //    tgtEntries iterators, and neither srcEntries nor
        //    tgtEntries have duplicates, so this is safe.
        // 3. srcCur->first.col == tgtCur->first.col, so add the two entries.
        if (srcCur != srcEnd) {
          if (srcCur->first.col == tgtCur->first.col) {
            tgtCur->second += srcCur->second;
          }
          else {
            // tgtCur is (optimally) right before where we want to put
            // the new entry, since srcCur points to the first entry
            // in srcEntries whose column index is greater than that
            // of the entry to which tgtCur points.
            (void) tgtEntries.insert (tgtCur, *srcCur);
          }
        } // if srcCur != srcEnd
      } // for each entry in the current row (lclRow) of *this
    } // if srcBeg != srcEnd
  }

  /// \brief Count the number of packets (bytes, in this case) needed
  ///   to pack the given row of the matrix.
  ///
  /// \param numPackets [out] Number of packets (bytes, in this case)
  ///   needed for the row.  Must be an int for MPI's sake.
  /// \param gblRow [in] Global index of the row to pack.
  /// \param comm [in] Communicator for packing.
  ///
  /// \return Error code; MPI_SUCESSS (0) if no error.
  int
  countPackRow (int& numPackets,
                const GO gblRow,
                const ::Teuchos::Comm<int>& comm,
                std::ostream* errStrm = NULL) const
  {
    using ::Tpetra::Details::countPackTriples;
    using ::Tpetra::Details::countPackTriplesCount;
    using std::endl;
    const char prefix[] = "Tpetra::Details::Impl::CooMatrixImpl::countPackRow: ";
#ifdef HAVE_TPETRACORE_MPI
    int errCode = MPI_SUCCESS;
#else
    int errCode = 0;
#endif // HAVE_TPETRACORE_MPI

    // Count the number of entries in the given row.
    const GO minGO = std::numeric_limits<GO>::min ();
    auto beg = this->entries_.lower_bound ({gblRow, minGO});
    auto end = this->entries_.upper_bound ({gblRow+1, minGO});
    int numEnt = 0;
    for (auto iter = beg; iter != end; ++iter) {
      ++numEnt;
      if (numEnt == std::numeric_limits<int>::max ()) {
        if (errStrm != NULL) {
          *errStrm << prefix << "In (global) row " << gblRow << ", the number "
            "of entries thus far, numEnt = " << numEnt << ", has reached the "
            "maximum int value.  We need to store this count as int for MPI's "
            "sake." << endl;
        }
        return -1;
      }
    }

    int numPacketsForCount = 0; // output argument of countPackTriplesCount
    {
      errCode =
        countPackTriplesCount (comm, numPacketsForCount, errStrm);
      if (errCode != 0) {
        if (errStrm != NULL) {
          *errStrm << prefix << "countPackTriplesCount "
            "returned errCode = " << errCode << " != 0." << endl;
        }
        return errCode;
      }
      if (numPacketsForCount < 0) {
        if (errStrm != NULL) {
          *errStrm << prefix << "countPackTriplesCount returned "
            "numPacketsForCount = " << numPacketsForCount << " < 0." << endl;
        }
        return -1;
      }
    }

    int numPacketsForTriples = 0; // output argument of countPackTriples
    {
      errCode = countPackTriples<SC, GO> (numEnt, comm, numPacketsForTriples);
      TEUCHOS_TEST_FOR_EXCEPTION
        (errCode != 0, std::runtime_error, prefix << "countPackTriples "
         "returned errCode = " << errCode << " != 0.");
      TEUCHOS_TEST_FOR_EXCEPTION
        (numPacketsForTriples < 0, std::logic_error, prefix << "countPackTriples "
         "returned numPacketsForTriples = " << numPacketsForTriples << " < 0.");
    }

    numPackets = numPacketsForCount + numPacketsForTriples;
    return errCode;
  }

  /// \brief Pack the given row of the matrix.
  ///
  /// \param outBuf [out] Output pack buffer.
  /// \param outBufSize [out] Total output buffer size in bytes.
  /// \param outBufCurPos [in/out] Current position from which to
  ///   start writing to the output buffer.  This corresponds to the
  ///   'position' in/out argument of MPI_Pack.
  /// \param comm [in] The communicator (MPI wants this).
  ///
  /// \param gblRowInds [in/out] Temporary space for row indices.
  /// \param gblColInds [in/out] Temporary space for column indices.
  /// \param vals [in/out] Temporary space for matrix values.
  ///
  /// \param gblRow [in] Global index of the row to pack.
  void
  packRow (packet_type outBuf[],
           const int outBufSize,
           int& outBufCurPos, // in/out argument
           const ::Teuchos::Comm<int>& comm,
           std::vector<GO>& gblRowInds,
           std::vector<GO>& gblColInds,
           std::vector<SC>& vals,
           const GO gblRow) const
  {
    using ::Tpetra::Details::packTriples;
    using ::Tpetra::Details::packTriplesCount;
    const char prefix[] = "Tpetra::Details::Impl::CooMatrixImpl::packRow: ";

    const GO minGO = std::numeric_limits<GO>::min ();
    auto beg = this->entries_.lower_bound ({gblRow, minGO});
    auto end = this->entries_.upper_bound ({gblRow+1, minGO});

    // This doesn't actually deallocate.  Only swapping with an empty
    // std::vector does that.
    gblRowInds.resize (0);
    gblColInds.resize (0);
    vals.resize (0);

    int numEnt = 0;
    for (auto iter = beg; iter != end; ++iter) {
      gblRowInds.push_back (iter->first.row);
      gblColInds.push_back (iter->first.col);
      vals.push_back (iter->second);
      ++numEnt;
      TEUCHOS_TEST_FOR_EXCEPTION
        (numEnt == std::numeric_limits<int>::max (), std::runtime_error, prefix
         << "In (global) row " << gblRow << ", the number of entries thus far, "
         "numEnt = " << numEnt << ", has reached the maximum int value.  "
         "We need to store this count as int for MPI's sake.");
    }

    {
      const int errCode =
        packTriplesCount (numEnt, outBuf, outBufSize, outBufCurPos, comm);
      TEUCHOS_TEST_FOR_EXCEPTION
        (errCode != 0, std::runtime_error, prefix
         << "In (global) row " << gblRow << ", packTriplesCount returned "
         "errCode = " << errCode << " != 0.");
    }
    {
      const int errCode =
        packTriples (gblRowInds.data (),
                     gblColInds.data (),
                     vals.data (),
                     numEnt,
                     outBuf,
                     outBufSize,
                     outBufCurPos, // in/out argument
                     comm);
      TEUCHOS_TEST_FOR_EXCEPTION
        (errCode != 0, std::runtime_error, prefix << "In (global) row "
         << gblRow << ", packTriples returned errCode = " << errCode
         << " != 0.");
    }
  }

  /// \brief Get the global row indices on this process, sorted and
  ///   made unique, and return the minimum global row index on this
  ///   process.
  ///
  /// \param rowInds [out] The global row indices on this process.
  ///
  /// \return Minimum global row index on this process.
  GO
  getMyGlobalRowIndices (std::vector<GO>& rowInds) const
  {
    rowInds.clear ();

    GO lclMinRowInd = std::numeric_limits<GO>::max (); // compute local min
    for (typename entries_type::const_iterator iter = this->entries_.begin ();
         iter != this->entries_.end (); ++iter) {
      const GO rowInd = iter->first.row;
      if (rowInd < lclMinRowInd) {
        lclMinRowInd = rowInd;
      }
      if (rowInds.empty () || rowInds.back () != rowInd) {
        rowInds.push_back (rowInd); // don't insert duplicates
      }
    }
    return lclMinRowInd;
  }

  template<class OffsetType>
  void
  buildCrs (std::vector<OffsetType>& rowOffsets,
            GO gblColInds[],
            SC vals[]) const
  {
    static_assert (std::is_integral<OffsetType>::value,
                   "OffsetType must be a built-in integer type.");

    // clear() doesn't generally free storage; it just resets the
    // length.  Thus, we reuse any existing storage here.
    rowOffsets.clear ();

    const std::size_t numEnt = this->getLclNumEntries ();
    if (numEnt == 0) {
      rowOffsets.push_back (0);
    }
    else {
      typename entries_type::const_iterator iter = this->entries_.begin ();
      GO prevGblRowInd = iter->first.row;

      OffsetType k = 0;
      for ( ; iter != this->entries_.end (); ++iter, ++k) {
        const GO gblRowInd = iter->first.row;
        if (k == 0 || gblRowInd != prevGblRowInd) {
          // The row offsets array always has at least one entry.  The
          // first entry is always zero, and the last entry is always
          // the number of matrix entries.
          rowOffsets.push_back (k);
          prevGblRowInd = gblRowInd;
        }
        gblColInds[k] = iter->first.col;

        static_assert (std::is_same<typename std::decay<decltype (iter->second)>::type, SC>::value,
                       "The type of iter->second != SC.");
        vals[k] = iter->second;
      }
      rowOffsets.push_back (static_cast<OffsetType> (numEnt));
    }
  }

  /// \brief Build a locally indexed version of CRS storage.
  ///
  /// Build a locally indexed version of compressed row sparse (CRS,
  /// also known as "compressed sparse row," CSR) storage.
  ///
  /// \param rowOffsets [out] Row offsets.
  /// \param lclColInds [out] The matrix's local column indices; must
  ///   have at least getLclNumEntries() entries.
  /// \param vals [out] The matrix's values; must have at least
  ///   getLclNumEntries() entries.
  /// \param gblToLcl [in] Closure that can convert a global
  ///   <i>column</i> index to a local column index.
  template<class OffsetType, class LO>
  void
  buildLocallyIndexedCrs (std::vector<OffsetType>& rowOffsets,
                          LO lclColInds[],
                          SC vals[],
                          std::function<LO (const GO)> gblToLcl) const
  {
    static_assert (std::is_integral<OffsetType>::value,
                   "OffsetType must be a built-in integer type.");
    static_assert (std::is_integral<LO>::value,
                   "LO must be a built-in integer type.");

    // clear() doesn't generally free storage; it just resets the
    // length.  Thus, we reuse any existing storage here.
    rowOffsets.clear ();

    const std::size_t numEnt = this->getLclNumEntries ();
    if (numEnt == 0) {
      rowOffsets.push_back (0);
    }
    else {
      typename entries_type::const_iterator iter = this->entries_.begin ();
      GO prevGblRowInd = iter->first.row;

      OffsetType k = 0;
      for ( ; iter != this->entries_.end (); ++iter, ++k) {
        const GO gblRowInd = iter->first.row;
        if (k == 0 || gblRowInd != prevGblRowInd) {
          // The row offsets array always has at least one entry.  The
          // first entry is always zero, and the last entry is always
          // the number of matrix entries.
          rowOffsets.push_back (k);
          prevGblRowInd = gblRowInd;
        }
        lclColInds[k] = gblToLcl (iter->first.col);
        vals[k] = iter->second;
      }
      rowOffsets.push_back (static_cast<OffsetType> (numEnt));
    }
  }
};

} // namespace Impl

/// \brief Sparse matrix used only for file input / output.
///
/// This class stores a sparse matrix in coordinate format.  It is
/// meant only to help file input and output.  Thus, it implements
/// DistObject, but does NOT implement RowMatrix or even Operator.
///
/// Unlike other DistObject subclasses in Tpetra, this class'
/// constructor need not necessarily take a Map.  If the class does
/// NOT have a Map, it builds its Map at fillComplete(), after
/// construction, as a function of the input indices.
///
/// Users are only allowed to insert matrix entries if the class does
/// NOT have a Map.  Insertion is local to each process.  Users insert
/// entries by calling sumIntoGlobalValue().  Each entry is a "triple"
/// consisting of a global row index, a global column index, and a
/// matrix value.  Users call fillComplete() when they are done
/// inserting entries.  At that point, the object builds its Map.
/// (Unlike CrsMatrix, this class does NOT allow multiple resumeFill()
/// / fillComplete() cycles.)
///
/// Once the class has a Map, users may apply DistObject methods like
/// doImport() and doExport() to redistribute the data.  The target of
/// an Import or Export must have a Map, and users may not have
/// inserted entries into it.
///
/// Here is an example of how to use this class:
///
/// \code
/// Tpetra::Details::CooMatrix<double, int, long long> A_in;
/// for (size_t k = 0; k < numEntriesToInsert; ++k) {
///   long long gblRowInd, gblColInd;
///   double val;
///   // We don't implement this here.  You have to do it.
///   readEntryFromFile (&gblRowInd, &gblColInd, &val);
///
///   A_in.sumIntoGlobalValue (gblRowInd, gblColInd, val);
/// }
///
/// // You are responsible for supplying the output matrix's
/// // communicator and Map.
/// A_in.fillComplete (comm);
/// Teuchos::RCP<const Tpetra::Map<int, long long> > outMap = ...;
///
/// Tpetra::Details::CooMatrix<double, int, long long> A_out (rowMap);
/// Tpetra::Export<int, long long> exporter (A_in.getMap (), outMap);
///
/// A_out.doExport (A_in, exporter, Tpetra::ADD);
/// \endcode
template<class SC,
         class LO = ::Tpetra::DistObject<char>::local_ordinal_type,
         class GO = ::Tpetra::DistObject<char>::global_ordinal_type,
         class NT = ::Tpetra::DistObject<char>::node_type>
class CooMatrix : public ::Tpetra::DistObject<char, LO, GO, NT> {
public:
  //! This class transfers data as bytes (MPI_BYTE).
  typedef char packet_type;
  //! Type of each entry (value) in the sparse matrix.
  typedef SC scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef NT node_type;
  typedef typename NT::device_type device_type;
  //! Type of the Map specialization to give to the constructor.
  typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

private:
  //! Type of the base class of CooMatrix.
  typedef ::Tpetra::DistObject<packet_type, local_ordinal_type,
                               global_ordinal_type, node_type> base_type;

  //! Implementation (internal data structure).
  Impl::CooMatrixImpl<SC, GO> impl_;

public:
  /// \brief Default constructor.
  ///
  /// This creates the object with a null Map.  Users may insert
  /// entries into this object, until they call fillComplete.  This
  /// object may NOT be the target of an Import or Export operation.
  CooMatrix () :
    base_type (::Teuchos::null),
    localError_ (new bool (false)),
    errs_ (new std::shared_ptr<std::ostringstream> ()) // ptr to a null ptr
  {}

  /// \brief Constructor that takes a Map.
  ///
  /// \param map [in] Input Map.
  ///
  /// If \c map is nonnull, then users may use this object ONLY as the
  /// target of an Import or Export operation.  They may NOT insert
  /// entries into it.
  CooMatrix (const ::Teuchos::RCP<const map_type>& map) :
    base_type (map),
    localError_ (new bool (false)),
    errs_ (new std::shared_ptr<std::ostringstream> ()) // ptr to a null ptr
  {}

  //! Destructor (virtual for memory safety of derived classes).
  virtual ~CooMatrix () {}

  /// \brief Insert one entry locally into the sparse matrix, if it
  ///   does not exist there yet.  If it does exist, sum the values.
  ///
  /// \param gblRowInd [in] Global row index of the entry to insert.
  /// \param gblColInd [in] Global column index of the entry to insert.
  /// \param val [in] Value of the matrix entry to insert / sum.
  void
  sumIntoGlobalValue (const GO gblRowInd,
                      const GO gblColInd,
                      const SC& val)
  {
    this->impl_.sumIntoGlobalValue (gblRowInd, gblColInd, val);
  }

  /// \brief Insert multiple entries locally into the sparse matrix.
  ///
  /// This works like multiple calls to sumIntoGlobalValue.
  ///
  /// \param gblRowInd [in] Global row indices of the entries to insert.
  /// \param gblColInd [in] Global column indices of the entries to insert.
  /// \param val [in] Values of the matrix entries to insert / sum.
  /// \param numEnt [in] Number of entries to insert.
  void
  sumIntoGlobalValues (const GO gblRowInds[],
                       const GO gblColInds[],
                       const SC vals[],
                       const std::size_t numEnt)
  {
    this->impl_.sumIntoGlobalValues (gblRowInds, gblColInds, vals, numEnt);
  }

  /// \brief Initializer-list overload of the above method (which see).
  void
  sumIntoGlobalValues (std::initializer_list<GO> gblRowInds,
                       std::initializer_list<GO> gblColInds,
                       std::initializer_list<SC> vals,
                       const std::size_t numEnt)
  {
    this->impl_.sumIntoGlobalValues (gblRowInds.begin (), gblColInds.begin (),
                                     vals.begin (), numEnt);
  }

  //! Number of entries in the sparse matrix on the calling process.
  std::size_t
  getLclNumEntries () const
  {
    return this->impl_.getLclNumEntries ();
  }

  template<class OffsetType>
  void
  buildCrs (::Kokkos::View<OffsetType*, device_type>& rowOffsets,
            ::Kokkos::View<GO*, device_type>& gblColInds,
            ::Kokkos::View<typename ::Kokkos::ArithTraits<SC>::val_type*, device_type>& vals)
  {
    static_assert (std::is_integral<OffsetType>::value,
                   "OffsetType must be a built-in integer type.");
    using ::Kokkos::create_mirror_view;
    using ::Kokkos::deep_copy;
    using ::Kokkos::View;
    typedef typename ::Kokkos::ArithTraits<SC>::val_type ISC;

    const std::size_t numEnt = this->getLclNumEntries ();

    gblColInds = View<GO*, device_type> ("gblColInds", numEnt);
    vals = View<ISC*, device_type> ("vals", numEnt);
    auto gblColInds_h = create_mirror_view (gblColInds);
    auto vals_h = create_mirror_view (vals);

    std::vector<std::size_t> rowOffsetsSV;
    this->impl_.buildCrs (rowOffsetsSV,
                          gblColInds_h.data (),
                          vals_h.data ());
    rowOffsets =
      View<OffsetType*, device_type> ("rowOffsets", rowOffsetsSV.size ());
    typename View<OffsetType*, device_type>::HostMirror
      rowOffsets_h (rowOffsetsSV.data (), rowOffsetsSV.size ());
    deep_copy (rowOffsets, rowOffsets_h);

    deep_copy (gblColInds, gblColInds_h);
    deep_copy (vals, vals_h);
  }

  /// \brief Tell the matrix that you are done inserting entries
  ///   locally, and that the matrix should build its Map now.
  ///
  /// This is the preferred version of fillComplete().
  ///
  /// \pre The matrix does not yet have a Map.
  /// \pre fillComplete() has never been called before on this object.
  /// \post The matrix has a Map.
  ///
  /// \param comm [in] Input communicator to use for the Map.
  void
  fillComplete (const ::Teuchos::RCP<const ::Teuchos::Comm<int> >& comm)
  {
    if (comm.is_null ()) {
      this->map_ = ::Teuchos::null;
    }
    else {
      this->map_ = this->buildMap (comm);
    }
  }

  /// \brief Special version of fillComplete that assumes that the
  ///   matrix already has a Map, and reuses its communicator to
  ///   create a new Map.
  ///
  /// DO NOT call this method unless you know what you are doing.
  ///
  /// \pre The matrix DOES have a (nonnull) Map with a nonnull
  ///   communicator.
  void
  fillComplete ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->getMap ().is_null (), std::runtime_error, "Tpetra::Details::"
       "CooMatrix::fillComplete: This object does not yet have a Map.  "
       "You must call the version of fillComplete "
       "that takes an input communicator.");
    this->fillComplete (this->getMap ()->getComm ());
  }

  /// \brief Whether this object had an error on the calling process.
  ///
  /// Import and Export operations using this object as the target of
  /// the Import or Export may incur local errors, for example due to
  /// invalid rows or incorrectly sized data.  On local error
  /// detection, we don't want to throw an exception right away,
  /// because not all processes may throw an exception; this
  /// inconsistency across processes can result in deadlock or put
  /// Tpetra in an incorrect state.  Instead, we set a local error
  /// flag on the affected process, and ignore the incorrect data.  If
  /// you want to check whether any process experienced an error, you
  /// must do a reduction or all-reduce over this flag.  Every time
  /// you initiate a new Import or Export with this object as the
  /// target, we clear this flag.
  bool localError () const {
    return *localError_;
  }

  /// \brief The current stream of error messages.
  ///
  /// This is only nonempty on the calling process if localError()
  /// returns true.  In that case, it stores a stream of
  /// human-readable, endline-separated error messages encountered
  /// during an Import or Export cycle.  Every time you initiate a new
  /// Import or Export with this object as the target, we clear this
  /// stream.
  ///
  /// If you want to print this, you are responsible for ensuring that
  /// it is valid for the calling MPI process to print to whatever
  /// output stream you use.  On some MPI implementations, you may
  /// need to send the string to Process 0 in MPI_COMM_WORLD for
  /// printing.
  ///
  /// Note to developers: we clear the stream at the beginning of
  /// checkSizes().
  std::string errorMessages () const {
    return ((*errs_).get () == NULL) ? std::string ("") : (*errs_)->str ();
  }

private:
  /// \brief Whether this object on the calling process is in an error
  ///   state.
  ///
  /// See the documentation of localError() for details.
  ///
  /// This pointer is always nonnull.  Using a pointer rather than a
  /// \c bool value here ensures that all views of this object have
  /// access to the error state, because all views have the same
  /// (nonnull at construction) pointer.
  ///
  /// Note to developers: we clear this flag at the beginning of
  /// checkSizes().
  std::shared_ptr<bool> localError_;

  /// \brief Stream of error messages.
  ///
  /// The outer pointer is always nonnull, but the inner pointer is
  /// only nonnull if localError_ is true.  Using a pointer to a
  /// pointer ensures that all views of this object have access to the
  /// error stream, because all views have the same (nonnull at
  /// construction) outer pointer.
  std::shared_ptr<std::shared_ptr<std::ostringstream> > errs_;

  //! Mark that a local error occurred, and get a stream for reporting it.
  std::ostream&
  markLocalErrorAndGetStream ()
  {
    * (this->localError_) = true;
    if ((*errs_).get () == NULL) {
      *errs_ = std::shared_ptr<std::ostringstream> (new std::ostringstream ());
    }
    return **errs_;
  }

public:
  /// \brief One-line descriptiion of this object; overrides
  ///   Teuchos::Describable method.
  virtual std::string description () const {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "\"Tpetra::Details::CooMatrix\": {"
       << "SC: " << TypeNameTraits<SC>::name ()
       << ", LO: " << TypeNameTraits<LO>::name ()
       << ", GO: " << TypeNameTraits<GO>::name ()
       << ", NT: " << TypeNameTraits<NT>::name ();
    if (this->getObjectLabel () != "") {
      os << ", Label: \"" << this->getObjectLabel () << "\"";
    }
    os << ", Has Map: " << (this->map_.is_null () ? "false" : "true")
       << "}";
    return os.str ();
  }

  /// \brief Print a descriptiion of this object to the given output
  ///   stream; overrides Teuchos::Describable method.
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const
  {
    using ::Tpetra::Details::gathervPrint;
    using ::Teuchos::EVerbosityLevel;
    using ::Teuchos::OSTab;
    using ::Teuchos::TypeNameTraits;
    using ::Teuchos::VERB_DEFAULT;
    using ::Teuchos::VERB_LOW;
    using ::Teuchos::VERB_MEDIUM;
    using std::endl;

    const auto vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    auto comm = this->getMap ().is_null () ?
      Teuchos::null : this->getMap ()->getComm ();
    const int myRank = comm.is_null () ? 0 : comm->getRank ();
    // const int numProcs = comm.is_null () ? 1 : comm->getSize ();

    if (vl != Teuchos::VERB_NONE) {
      // Convention is to start describe() implementations with a tab.
      OSTab tab0 (out);
      if (myRank == 0) {
        out << "\"Tpetra::Details::CooMatrix\":" << endl;
      }
      OSTab tab1 (out);
      if (myRank == 0) {
        out << "Template parameters:" << endl;
        {
          OSTab tab2 (out);
          out << "SC: " << TypeNameTraits<SC>::name () << endl
              << "LO: " << TypeNameTraits<LO>::name () << endl
              << "GO: " << TypeNameTraits<GO>::name () << endl
              << "NT: " << TypeNameTraits<NT>::name () << endl;
        }
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
        out << "Has Map: " << (this->map_.is_null () ? "false" : "true") << endl;
      } // if myRank == 0

      // Describe the Map, if it exists.
      if (! this->map_.is_null ()) {
        if (myRank == 0) {
          out << "Map:" << endl;
        }
        OSTab tab2 (out);
        this->map_->describe (out, vl);
      }

      // At verbosity > VERB_LOW, each process prints something.
      if (vl > VERB_LOW) {
        std::ostringstream os;
        os << "Process " << myRank << ":" << endl;
        //OSTab tab3 (os);

        const std::size_t numEnt = this->impl_.getLclNumEntries ();
        os << " Local number of entries: " << numEnt << endl;
        if (vl > VERB_MEDIUM) {
          os << " Local entries:" << endl;
          //OSTab tab4 (os);
          this->impl_.forAllEntries ([&os] (const GO row, const GO col, const SC& val) {
              os << "  {"
                 << "row: " << row
                 << ", col: " << col
                 << ", val: " << val
                 << "}" << endl;
            });
        }
        gathervPrint (out, os.str (), *comm);
      }
    } // vl != Teuchos::VERB_NONE
  }

private:
  /// \brief Build the Map from the local sparse matrix entries.
  ///
  /// This method has collective semantics over the input communicator.
  ///
  /// \param comm [in] Communicator to use for the Map.
  ///
  /// \return The Map, whose global indices on the calling process are
  ///   exactly those in \c entries.
  Teuchos::RCP<const map_type>
  buildMap (const ::Teuchos::RCP<const ::Teuchos::Comm<int> >& comm)
  {
    using ::Teuchos::outArg;
    using ::Teuchos::rcp;
    using ::Teuchos::REDUCE_MIN;
    using ::Teuchos::reduceAll;
    typedef ::Tpetra::global_size_t GST;
    //const char prefix[] = "Tpetra::Details::CooMatrix::buildMap: ";

    // Processes where comm is null don't participate in the Map.
    if (comm.is_null ()) {
      return ::Teuchos::null;
    }

    // mfh 17 Jan 2017: We just happen to use row indices, because
    // that's what Tpetra::CrsMatrix currently uses.  That's probably
    // not the best thing to use, but it's not bad for commonly
    // encountered matrices.  A better more general approach might be
    // to hash (row index, column index) pairs to a global index.  One
    // could make that unique by doing a parallel scan at map
    // construction time.

    std::vector<GO> rowInds;
    const GO lclMinGblRowInd = this->impl_.getMyGlobalRowIndices (rowInds);

    // Compute the globally min row index for the "index base."
    GO gblMinGblRowInd = 0; // output argument
    reduceAll<int, GO> (*comm, REDUCE_MIN, lclMinGblRowInd,
                        outArg (gblMinGblRowInd));
    const GO indexBase = gblMinGblRowInd;
    const GST INV = Tpetra::Details::OrdinalTraits<GST>::invalid ();
    return rcp (new map_type (INV, rowInds.data (), rowInds.size (),
                              indexBase, comm));
  }

protected:
  /// \brief By returning 0, tell DistObject that this class may not
  ///   necessarily have a constant number of "packets" per local
  ///   index.
  virtual size_t constantNumberOfPackets () const {
    return static_cast<size_t> (0);
  }

  /// \brief Compare the source and target (\e this) objects for compatibility.
  ///
  /// \return True if they are compatible, else false.
  virtual bool
  checkSizes (const ::Tpetra::SrcDistObject& source)
  {
    using std::endl;
    typedef CooMatrix<SC, LO, GO, NT> this_COO_type;
    const char prefix[] = "Tpetra::Details::CooMatrix::checkSizes: ";

    const this_COO_type* src = dynamic_cast<const this_COO_type* > (&source);

    if (src == NULL) {
      std::ostream& err = markLocalErrorAndGetStream ();
      err << prefix << "The source object of the Import or Export "
        "must be a CooMatrix with the same template parameters as the "
        "target object." << endl;
    }
    else if (this->map_.is_null ()) {
      std::ostream& err = markLocalErrorAndGetStream ();
      err << prefix << "The target object of the Import or Export "
        "must be a CooMatrix with a nonnull Map." << endl;
    }
    return ! (* (this->localError_));
  }

  //! Kokkos::Device specialization for DistObject communication buffers.
  using buffer_device_type =
    typename ::Tpetra::DistObject<char, LO, GO, NT>::buffer_device_type;

  virtual void
  copyAndPermute
  (const ::Tpetra::SrcDistObject& sourceObject,
   const size_t numSameIDs,
   const Kokkos::DualView<const LO*,
     buffer_device_type>& permuteToLIDs,
   const Kokkos::DualView<const LO*,
     buffer_device_type>& permuteFromLIDs,
   const CombineMode /* CM */)
  {
    using std::endl;
    using this_COO_type = CooMatrix<SC, LO, GO, NT>;
    const char prefix[] = "Tpetra::Details::CooMatrix::copyAndPermute: ";

    // There's no communication in this method, so it's safe just to
    // return on error.

    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "The target object of the Import or Export is "
        "already in an error state." << endl;
      return;
    }

    const this_COO_type* src = dynamic_cast<const this_COO_type*> (&sourceObject);
    if (src == nullptr) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "Input argument 'sourceObject' is not a CooMatrix."
          << endl;
      return;
    }

    const size_t numPermuteIDs =
      static_cast<size_t> (permuteToLIDs.extent (0));
    if (numPermuteIDs != static_cast<size_t> (permuteFromLIDs.extent (0))) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "permuteToLIDs.extent(0) = "
          << numPermuteIDs << " != permuteFromLIDs.extent(0) = "
          << permuteFromLIDs.extent (0) << "." << endl;
      return;
    }
    if (sizeof (int) <= sizeof (size_t) &&
        numPermuteIDs > static_cast<size_t> (std::numeric_limits<int>::max ())) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "numPermuteIDs = " << numPermuteIDs
          << ", a size_t, overflows int." << endl;
      return;
    }

    // Even though this is an std::set, we can start with numSameIDs
    // just by iterating through the first entries of the set.

    if (sizeof (int) <= sizeof (size_t) &&
        numSameIDs > static_cast<size_t> (std::numeric_limits<int>::max ())) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "numSameIDs = " << numSameIDs
          << ", a size_t, overflows int." << endl;
      return;
    }

    //
    // Copy in entries from any initial rows with the same global row indices.
    //
    const LO numSame = static_cast<int> (numSameIDs);
    // Count of local row indices encountered here with invalid global
    // row indices.  If nonzero, something went wrong.  If something
    // did go wrong, we'll defer responding until the end of this
    // method, so we can print as much useful info as possible.
    LO numInvalidSameRows = 0;
    for (LO lclRow = 0; lclRow < numSame; ++lclRow) {
      // All numSame initial rows have the same global index in both
      // source and target Maps, so we only need to convert to global
      // once.
      const GO gblRow = this->map_->getGlobalElement (lclRow);
      if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
        ++numInvalidSameRows;
        continue;
      }
      else {
        this->impl_.mergeIntoRow (gblRow, src->impl_, gblRow);
      }
    }

    //
    // Copy in entries from remaining rows that are permutations, that
    // is, that live in both the source and target Maps, but aren't
    // included in the "same" list (see above).
    //
    const LO numPermute = static_cast<int> (numPermuteIDs);
    // Count of local "from" row indices encountered here with invalid
    // global row indices.  If nonzero, something went wrong.  If
    // something did go wrong, we'll defer responding until the end of
    // this method, so we can print as much useful info as possible.
    LO numInvalidRowsFrom = 0;
    // Count of local "to" row indices encountered here with invalid
    // global row indices.  If nonzero, something went wrong.  If
    // something did go wrong, we'll defer responding until the end of
    // this method, so we can print as much useful info as possible.
    LO numInvalidRowsTo = 0;

    TEUCHOS_ASSERT( ! permuteFromLIDs.need_sync_host () );
    TEUCHOS_ASSERT( ! permuteToLIDs.need_sync_host () );
    auto permuteFromLIDs_h = permuteFromLIDs.view_host ();
    auto permuteToLIDs_h = permuteToLIDs.view_host ();

    for (LO k = 0; k < numPermute; ++k) {
      const LO lclRowFrom = permuteFromLIDs_h[k];
      const LO lclRowTo = permuteToLIDs_h[k];
      const GO gblRowFrom = src->map_->getGlobalElement (lclRowFrom);
      const GO gblRowTo = this->map_->getGlobalElement (lclRowTo);

      bool bothConversionsValid = true;
      if (gblRowFrom == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
        ++numInvalidRowsFrom;
        bothConversionsValid = false;
      }
      if (gblRowTo == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
        ++numInvalidRowsTo;
        bothConversionsValid = false;
      }
      if (bothConversionsValid) {
        this->impl_.mergeIntoRow (gblRowTo, src->impl_, gblRowFrom);
      }
    }

    // Print info if any errors occurred.
    if (numInvalidSameRows != 0 || numInvalidRowsFrom != 0 ||
        numInvalidRowsTo != 0) {
      // Collect and print all the invalid input row indices, for the
      // "same," "from," and "to" lists.
      std::vector<std::pair<LO, GO> > invalidSameRows;
      invalidSameRows.reserve (numInvalidSameRows);
      std::vector<std::pair<LO, GO> > invalidRowsFrom;
      invalidRowsFrom.reserve (numInvalidRowsFrom);
      std::vector<std::pair<LO, GO> > invalidRowsTo;
      invalidRowsTo.reserve (numInvalidRowsTo);

      for (LO lclRow = 0; lclRow < numSame; ++lclRow) {
        // All numSame initial rows have the same global index in both
        // source and target Maps, so we only need to convert to global
        // once.
        const GO gblRow = this->map_->getGlobalElement (lclRow);
        if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          invalidSameRows.push_back ({lclRow, gblRow});
        }
      }

      for (LO k = 0; k < numPermute; ++k) {
        const LO lclRowFrom = permuteFromLIDs_h[k];
        const LO lclRowTo = permuteToLIDs_h[k];
        const GO gblRowFrom = src->map_->getGlobalElement (lclRowFrom);
        const GO gblRowTo = this->map_->getGlobalElement (lclRowTo);

        if (gblRowFrom == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          invalidRowsFrom.push_back ({lclRowFrom, gblRowFrom});
        }
        if (gblRowTo == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          invalidRowsTo.push_back ({lclRowTo, gblRowTo});
        }
      }

      std::ostringstream os;
      if (numInvalidSameRows != 0) {
        os << "Invalid permute \"same\" (local, global) index pairs: [";
        for (std::size_t k = 0; k < invalidSameRows.size (); ++k) {
          const auto& p = invalidSameRows[k];
          os << "(" << p.first << "," << p.second << ")";
          if (k + 1 < invalidSameRows.size ()) {
            os << ", ";
          }
        }
        os << "]" << std::endl;
      }
      if (numInvalidRowsFrom != 0) {
        os << "Invalid permute \"from\" (local, global) index pairs: [";
        for (std::size_t k = 0; k < invalidRowsFrom.size (); ++k) {
          const auto& p = invalidRowsFrom[k];
          os << "(" << p.first << "," << p.second << ")";
          if (k + 1 < invalidRowsFrom.size ()) {
            os << ", ";
          }
        }
        os << "]" << std::endl;
      }
      if (numInvalidRowsTo != 0) {
        os << "Invalid permute \"to\" (local, global) index pairs: [";
        for (std::size_t k = 0; k < invalidRowsTo.size (); ++k) {
          const auto& p = invalidRowsTo[k];
          os << "(" << p.first << "," << p.second << ")";
          if (k + 1 < invalidRowsTo.size ()) {
            os << ", ";
          }
        }
        os << "]" << std::endl;
      }

      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << os.str ();
      return;
    }
  }

  virtual void
  packAndPrepare
  (const ::Tpetra::SrcDistObject& sourceObject,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& exportLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type>& exports,
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID,
   size_t& constantNumPackets)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using std::endl;
    using this_COO_type = CooMatrix<SC, LO, GO, NT>;
    const char prefix[] = "Tpetra::Details::CooMatrix::packAndPrepare: ";
    const char suffix[] = "  This should never happen.  "
      "Please report this bug to the Tpetra developers.";

    // Tell the caller that different rows may have different numbers
    // of matrix entries.
    constantNumPackets = 0;

    const this_COO_type* src = dynamic_cast<const this_COO_type*> (&sourceObject);
    if (src == nullptr) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "Input argument 'sourceObject' is not a CooMatrix."
          << endl;
    }
    else if (* (src->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "The source (input) object of the Import or Export "
        "is already in an error state on this process."
          << endl;
    }
    else if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "The target (output, \"this\") object of the Import "
        "or Export is already in an error state on this process." << endl;
    }
    // Respond to detected error(s) by resizing 'exports' to zero (so
    // we won't be tempted to read it later), and filling
    // numPacketsPerLID with zeros.
    if (* (this->localError_)) {
      // Resize 'exports' to zero, so we won't be tempted to read it.
      Details::reallocDualViewIfNeeded (exports, 0, "CooMatrix exports");
      // Trick to get around const DualView& being const.
      {
        auto numPacketsPerLID_tmp = numPacketsPerLID;
        numPacketsPerLID_tmp.sync_host ();
        numPacketsPerLID_tmp.modify_host ();
      }
      // Fill numPacketsPerLID with zeros.
      Kokkos::deep_copy (numPacketsPerLID.h_view, static_cast<size_t> (0));
      return;
    }

    const size_t numExports = exportLIDs.extent (0);
    if (numExports == 0) {
      Details::reallocDualViewIfNeeded (exports, 0, exports.h_view.label ());
      return; // nothing to send
    }
    RCP<const Comm<int> > comm = src->getMap ().is_null () ?
      Teuchos::null : src->getMap ()->getComm ();
    if (comm.is_null () || comm->getSize () == 1) {
      if (numExports != static_cast<size_t> (0)) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "The input communicator is either null or "
          "has only one process, but numExports = " << numExports << " != 0."
            << suffix << endl;
        return;
      }
      // Don't go into the rest of this method unless there are
      // actually processes other than the calling process.  This is
      // because the pack and unpack functions only have nonstub
      // implementations if building with MPI.
      return;
    }

    numPacketsPerLID.sync_host ();
    numPacketsPerLID.modify_host ();

    TEUCHOS_ASSERT( ! exportLIDs.need_sync_host () );
    auto exportLIDs_h = exportLIDs.view_host ();

    int totalNumPackets = 0;
    size_t numInvalidRowInds = 0;
    std::ostringstream errStrm; // current loop iteration's error messages
    for (size_t k = 0; k < numExports; ++k) {
      const LO lclRow = exportLIDs_h[k];
      // We're packing the source object's data, so we need to use the
      // source object's Map to convert from local to global indices.
      const GO gblRow = src->map_->getGlobalElement (lclRow);
      if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
        // Mark the error later; just count for now.
        ++numInvalidRowInds;
        numPacketsPerLID.h_view[k] = 0;
        continue;
      }

      // Count the number of bytes needed to pack the current row of
      // the source object.
      int numPackets = 0;
      const int errCode =
        src->impl_.countPackRow (numPackets, gblRow, *comm, &errStrm);
      if (errCode != 0) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << errStrm.str () << endl;
        numPacketsPerLID.h_view[k] = 0;
        continue;
      }

      // Make sure that the total number of packets fits in int.
      // MPI requires this.
      const long long newTotalNumPackets =
        static_cast<long long> (totalNumPackets) +
        static_cast<long long> (numPackets);
      if (newTotalNumPackets >
          static_cast<long long> (std::numeric_limits<int>::max ())) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "The new total number of packets "
            << newTotalNumPackets << " does not fit in int." << endl;
        // At this point, we definitely cannot continue.  In order to
        // leave the output arguments in a rational state, we zero out
        // all remaining entries of numPacketsPerLID before returning.
        for (size_t k2 = k; k2 < numExports; ++k2) {
          numPacketsPerLID.h_view[k2] = 0;
        }
        return;
      }
      numPacketsPerLID.h_view[k] = static_cast<size_t> (numPackets);
      totalNumPackets = static_cast<int> (newTotalNumPackets);
    }

    // If we found invalid row indices in exportLIDs, go back,
    // collect, and report them.
    if (numInvalidRowInds != 0) {
      std::vector<std::pair<LO, GO> > invalidRowInds;
      for (size_t k = 0; k < numExports; ++k) {
        const LO lclRow = exportLIDs_h[k];
        // We're packing the source object's data, so we need to use
        // the source object's Map to convert from local to global
        // indices.
        const GO gblRow = src->map_->getGlobalElement (lclRow);
        if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          invalidRowInds.push_back ({lclRow, gblRow});
        }
      }
      std::ostringstream os;
      os << prefix << "We found " << numInvalidRowInds << " invalid row ind"
         << (numInvalidRowInds != static_cast<size_t> (1) ? "ices" : "ex")
         << " out of " << numExports << " in exportLIDs.  Here is the list "
         << " of invalid row indices: [";
      for (size_t k = 0; k < invalidRowInds.size (); ++k) {
        os << "(LID: " << invalidRowInds[k].first << ", GID: "
           << invalidRowInds[k].second << ")";
        if (k + 1 < invalidRowInds.size ()) {
          os << ", ";
        }
      }
      os << "].";

      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << os.str () << std::endl;
      return;
    }

    {
      const bool reallocated =
        Details::reallocDualViewIfNeeded (exports, totalNumPackets,
                                          "CooMatrix exports");
      if (reallocated) {
        exports.sync_host (); // make sure alloc'd on host
      }
    }
    exports.modify_host ();

    // FIXME (mfh 17 Jan 2017) packTriples wants three arrays, not a
    // single array of structs.  For now, we just copy.
    std::vector<GO> gblRowInds;
    std::vector<GO> gblColInds;
    std::vector<SC> vals;

    int outBufCurPos = 0;
    packet_type* outBuf = exports.h_view.data ();
    for (size_t k = 0; k < numExports; ++k) {
      const LO lclRow = exportLIDs.h_view[k];
      // We're packing the source object's data, so we need to use the
      // source object's Map to convert from local to global indices.
      const GO gblRow = src->map_->getGlobalElement (lclRow);
      // Pack the current row of the source object.
      src->impl_.packRow (outBuf, totalNumPackets, outBufCurPos, *comm,
                          gblRowInds, gblColInds, vals, gblRow);
    }
  }

  virtual void
  unpackAndCombine
  (const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& importLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type> imports,
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID,
   const size_t /* constantNumPackets */,
   const ::Tpetra::CombineMode /* combineMode */)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using std::endl;
    const char prefix[] = "Tpetra::Details::CooMatrix::unpackAndCombine: ";
    const char suffix[] = "  This should never happen.  "
      "Please report this bug to the Tpetra developers.";

    TEUCHOS_ASSERT( ! importLIDs.need_sync_host () );
    auto importLIDs_h = importLIDs.view_host ();

    const std::size_t numImports = importLIDs.extent (0);
    if (numImports == 0) {
      return; // nothing to receive
    }
    else if (imports.extent (0) == 0) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "importLIDs.extent(0) = " << numImports << " != 0, "
          << "but imports.extent(0) = 0.  This doesn't make sense, because "
          << "for every incoming LID, CooMatrix packs at least the count of "
          << "triples associated with that LID, even if the count is zero.  "
          << "importLIDs = [";
      for (std::size_t k = 0; k < numImports; ++k) {
        err << importLIDs_h[k];
        if (k + 1 < numImports) {
          err << ", ";
        }
      }
      err << "].  " << suffix << endl;
      return;
    }

    RCP<const Comm<int> > comm = this->getMap ().is_null () ?
      Teuchos::null : this->getMap ()->getComm ();
    if (comm.is_null () || comm->getSize () == 1) {
      if (numImports != static_cast<size_t> (0)) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "The input communicator is either null or "
          "has only one process, but numImports = " << numImports << " != 0."
            << suffix << endl;
        return;
      }
      // Don't go into the rest of this method unless there are
      // actually processes other than the calling process.  This is
      // because the pack and unpack functions only have nonstub
      // implementations if building with MPI.
      return;
    }

    // Make sure that the length of 'imports' fits in int.
    // This is ultimately an MPI restriction.
    if (static_cast<size_t> (imports.extent (0)) >
        static_cast<size_t> (std::numeric_limits<int>::max ())) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "imports.extent(0) = "
          << imports.extent (0) << " does not fit in int." << endl;
      return;
    }
    const int inBufSize = static_cast<int> (imports.extent (0));

    if (imports.need_sync_host ()) {
      imports.sync_host ();
    }
    if (numPacketsPerLID.need_sync_host ()) {
      numPacketsPerLID.sync_host ();
    }
    auto imports_h = imports.view_host ();
    auto numPacketsPerLID_h = numPacketsPerLID.view_host ();

    // FIXME (mfh 17,24 Jan 2017) packTriples wants three arrays, not a
    // single array of structs.  For now, we just copy.
    std::vector<GO> gblRowInds;
    std::vector<GO> gblColInds;
    std::vector<SC> vals;

    const packet_type* inBuf = imports_h.data ();
    int inBufCurPos = 0;
    size_t numInvalidRowInds = 0;
    int errCode = 0;
    std::ostringstream errStrm; // for unpack* error output.
    for (size_t k = 0; k < numImports; ++k) {
      const LO lclRow = importLIDs_h(k);
      const GO gblRow = this->map_->getGlobalElement (lclRow);
      if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
        ++numInvalidRowInds;
        continue;
      }

      // Remember where we were, so we don't overrun the buffer
      // length.  inBufCurPos is an in/out argument of unpackTriples*.
      const int origInBufCurPos = inBufCurPos;

      int numEnt = 0; // output argument of unpackTriplesCount
      errCode = unpackTriplesCount (inBuf, inBufSize, inBufCurPos,
                                    numEnt, *comm, &errStrm);
      if (errCode != 0 || numEnt < 0 || inBufCurPos < origInBufCurPos) {
        std::ostream& err = this->markLocalErrorAndGetStream ();

        err << prefix << "In unpack loop, k=" << k << ": ";
        if (errCode != 0) {
          err << "  unpackTriplesCount returned errCode = " << errCode
              << " != 0." << endl;
        }
        if (numEnt < 0) {
          err << "  unpackTriplesCount returned errCode = 0, but numEnt = "
              << numEnt << " < 0." << endl;
        }
        if (inBufCurPos < origInBufCurPos) {
          err << "  After unpackTriplesCount, inBufCurPos = " << inBufCurPos
              << " < origInBufCurPos = " << origInBufCurPos << "." << endl;
        }
        err << "  unpackTriplesCount report: " << errStrm.str () << endl;
        err << suffix << endl;

        // We only continue in a debug build, because the above error
        // messages could consume too much memory and cause an
        // out-of-memory error, without actually printing.  Printing
        // everything is useful in a debug build, but not so much in a
        // release build.
#ifdef HAVE_TPETRA_DEBUG
        // Clear out the current error stream, so we don't accumulate
        // over loop iterations.
        errStrm.str ("");
        continue;
#else
        return;
#endif // HAVE_TPETRA_DEBUG
      }

      // FIXME (mfh 17,24 Jan 2017) packTriples wants three arrays,
      // not a single array of structs.  For now, we just copy.
      gblRowInds.resize (numEnt);
      gblColInds.resize (numEnt);
      vals.resize (numEnt);

      errCode = unpackTriples (inBuf, inBufSize, inBufCurPos,
                               gblRowInds.data (), gblColInds.data (),
                               vals.data (), numEnt, *comm, &errStrm);
      if (errCode != 0) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "unpackTriples returned errCode = "
            << errCode << " != 0.  It reports: " << errStrm.str ()
            << endl;
        // We only continue in a debug build, because the above error
        // messages could consume too much memory and cause an
        // out-of-memory error, without actually printing.  Printing
        // everything is useful in a debug build, but not so much in a
        // release build.
#ifdef HAVE_TPETRA_DEBUG
        // Clear out the current error stream, so we don't accumulate
        // over loop iterations.
        errStrm.str ("");
        continue;
#else
        return;
#endif // HAVE_TPETRA_DEBUG
      }
      this->sumIntoGlobalValues (gblRowInds.data (), gblColInds.data (),
                                 vals.data (), numEnt);
    }

    // If we found invalid row indices in exportLIDs, go back,
    // collect, and report them.
    if (numInvalidRowInds != 0) {
      // Mark the error now, before we possibly run out of memory.
      // The latter could raise an exception (e.g., std::bad_alloc),
      // but at least we would get the error state right.
      std::ostream& err = this->markLocalErrorAndGetStream ();

      std::vector<std::pair<LO, GO> > invalidRowInds;
      for (size_t k = 0; k < numImports; ++k) {
        const LO lclRow = importLIDs_h(k);
        const GO gblRow = this->map_->getGlobalElement (lclRow);
        if (gblRow == ::Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          invalidRowInds.push_back ({lclRow, gblRow});
        }
      }

      err << prefix << "We found " << numInvalidRowInds << " invalid row ind"
          << (numInvalidRowInds != static_cast<size_t> (1) ? "ices" : "ex")
          << " out of " << numImports << " in importLIDs.  Here is the list "
          << " of invalid row indices: [";
      for (size_t k = 0; k < invalidRowInds.size (); ++k) {
        err << "(LID: " << invalidRowInds[k].first << ", GID: "
            << invalidRowInds[k].second << ")";
        if (k + 1 < invalidRowInds.size ()) {
          err << ", ";
        }
      }
      err << "].";
      return;
    }
  }
};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COOMATRIX_HPP
