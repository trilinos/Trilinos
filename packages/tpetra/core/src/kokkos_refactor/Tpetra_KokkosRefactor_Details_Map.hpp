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

#ifndef TPETRA_KOKKOSREFACTOR_DETAILS_MAP_HPP
#define TPETRA_KOKKOSREFACTOR_DETAILS_MAP_HPP

#include <Teuchos_Describable.hpp>
#include <Kokkos_UnorderedMap.hpp>

namespace Tpetra {
  namespace Details {

    /// \class Iota
    /// \brief Functor for the iota function.
    ///
    /// This functor fills the given View with integers 0, 1, ...,
    /// x.dimension_0() - 1.  It implements the iota() function below.
    /// Users should use the iota() function directly.
    template<class IntType, class DeviceType>
    class Iota {
    public:
      typedef DeviceType execution_space;

      Iota (const Kokkos::View<IntType*, DeviceType>& x,
            const IntType first) : x_ (x), first_ (first) {}

      KOKKOS_INLINE_FUNCTION void
      operator () (const typename DeviceType::size_type i) const {
        x_(i) = first_ + static_cast<IntType> (i);
      }

    private:
      Kokkos::View<IntType*, DeviceType> x_;
      const IntType first_;
    };

    //! Fill x with integers first, first + 1, ..., first + x.dimension_0() - 1.
    template<class IntType, class DeviceType>
    void
    iota (const Kokkos::View<IntType*, DeviceType>& x,
          const IntType first)
    {
      Kokkos::parallel_for (x.dimension_0 (), Iota<IntType, DeviceType> (x, first));
    }

    /// \struct MapData
    /// \brief Used by GlobalToLocalTableFiller.
    ///
    /// \note Structs given to functors just get copied over directly.
    template<class LO, class GO, class DeviceType>
    struct MapData {
      GO minMyGID; //!< My process' minimum GID
      GO maxMyGID; //!< My process' maximum GID
      /// \brief Count of failed inserts into GID->LID table.
      ///
      /// This should always be zero.  If not, there is a bug in the
      /// code that called fillGlobalToLocalTable().
      LO failedCount;
      /// \brief Count of duplicate inserts into GID->LID table.
      ///
      /// This should only ever be nonzero if the GID list on the
      /// calling process contains duplicates.  Map never supported
      /// that case well, in either Epetra or Tpetra.  For example, it
      /// would give those duplicate GIDs different LIDs.  If we ever
      /// want to fix this, for example by detecting and filtering out
      /// duplicate GIDs, we could use duplicateCount to detect that
      /// case.  If detected, we could then reassign LIDs to GIDs by
      /// merging duplicate GIDs, and rebuild the GID->LID table.  We
      /// would have to do this _before_ any MPI_Allreduce or
      /// MPI_*scan for computing the total GID count or offsets,
      /// since it would change the number of GIDs on the affected
      /// process(es).
      LO duplicateCount;

      //! Default constructor (to make this a proper value type).
      MapData () :
        minMyGID (0), maxMyGID (0), failedCount (0), duplicateCount (0)
      {}
    };

    /// \class GlobalToLocalTableFiller
    /// \brief Kokkos reduce functor for filling GID->LID table and
    ///   computing min and max GID
    ///
    /// The functor is a reduce functor, not a for functor, because it
    /// also conveniently computes the min and max locally owned GID.
    template<class LO, class GO, class DeviceType>
    class GlobalToLocalTableFiller {
    public:
      typedef DeviceType execution_space;
      typedef typename DeviceType::size_type size_type;
      typedef MapData<LO, GO, DeviceType> value_type;

      /// \brief Constructor
      ///
      /// \param glMap [out] GID->LID table to fill; must already be
      ///   preallocated with sufficient space.
      /// \param entries [in] GIDs to put into GID->LID table.  These
      ///   are only the zero or more "noncontiguous GIDs" that follow
      ///   the initial sequence of one or more contiguous GIDs
      ///   [firstContiguousGID, lastContiguousGID].
      /// \param firstContiguousGID [in] The first contiguous GID.
      /// \param lastContiguousGID [in] The last contiguous GID.
      GlobalToLocalTableFiller (const Kokkos::UnorderedMap<GO, LO, DeviceType>& glMap,
                                const Kokkos::View<const GO*, DeviceType>& entries,
                                const GO firstContiguousGID,
                                const GO lastContiguousGID) :
        glMap_ (glMap),
        entries_ (entries),
        firstContiguousGID_ (firstContiguousGID),
        lastContiguousGID_ (lastContiguousGID)
      {}

      /// \brief Set the initial value of the reduction.
      ///
      /// Pre-assign to minMyGID the first seen contiguous GID, and to
      /// maxMyGID the last seen contiguous GID.  [first,last] form an
      /// inclusive range.  Map's constructor (that takes a GID list)
      /// should not call this functor if the calling process owns no
      /// GIDs.  Even if there is only one GID, that will be both the
      /// first and the last contiguous GID.
      KOKKOS_INLINE_FUNCTION void
      init (value_type& dst) const
      {
        dst.minMyGID = firstContiguousGID_;
        dst.maxMyGID = lastContiguousGID_;
        dst.failedCount = 0;
        dst.duplicateCount = 0;
      }

      /// \brief Combine two intermediate reduction results.
      ///
      /// This sets both the min and max GID, if necessary.
      KOKKOS_INLINE_FUNCTION void
      join (volatile value_type& dst ,
            const volatile value_type& src) const
      {
        if (src.maxMyGID > dst.maxMyGID) {
          dst.maxMyGID = src.maxMyGID;
        }
        if (src.minMyGID < dst.minMyGID) {
          dst.minMyGID = src.minMyGID;
        }
        dst.failedCount += src.failedCount;
        dst.duplicateCount += src.duplicateCount;
      }

      /// \brief Do this for every element of entries.
      ///
      /// Add (entries_(i), i) to the GID->LID lookup table, and
      /// update the min and max GID seen thus far in entries.
      KOKKOS_INLINE_FUNCTION void
      operator () (const size_type i, value_type& dst) const
      {
        // [firstContiguousGID_, lastContiguousGID_] as GIDs map to
        // [0, lastContiguousGID_ - firstContiguousGID_ + 1] as LIDs.
        const LO lid =
          static_cast<LO> (lastContiguousGID_ - firstContiguousGID_ + 1) +
          static_cast<LO> (i);
        const GO gid = entries_(i);

        if (gid > dst.maxMyGID) {
          dst.maxMyGID = gid;
        }
        if (gid < dst.minMyGID) {
          dst.minMyGID = gid;
        }
        // The table should not run out of space, but if it does, the
        // caller will try again.
        const Kokkos::UnorderedMapInsertResult result = glMap_.insert (gid, lid);

        if (result.failed ()) {
          dst.failedCount++;
        }
        else if (result.existing ()) {
          dst.duplicateCount++;
        }
      }

    private:
      //! The GID->LID table to fill.
      Kokkos::UnorderedMap<GO, LO, DeviceType> glMap_;
      //! The input list of GIDs (not including initial contiguous GIDs).
      Kokkos::View<const GO*, DeviceType> entries_;
      //! First contiguous GID in initial contiguous GIDs.
      const GO firstContiguousGID_;
      //! Last contiguous GID (inclusive) in initial contiguous GIDs.
      const GO lastContiguousGID_;
    };

    /// \brief Fill the GID->LID lookup table, and compute the local
    ///   min and max GID.  Return them through the returned struct.
    ///
    /// The "noncontiguous" (GID list) constructor of
    /// Tpetra::Details::Map calls this function.  It is an
    /// implementation detail of Tpetra and not meant for users.
    ///
    /// \param glMap [out] GID->LID table to fill; must already be
    ///   preallocated with sufficient space.
    /// \param entries [in] GIDs to put into GID->LID table.  These
    ///   are only the zero or more "noncontiguous GIDs" that follow
    ///   the initial sequence of one or more contiguous GIDs
    ///   [firstContiguousGID, lastContiguousGID].
    /// \param firstContiguousGID [in] The first contiguous GID.
    /// \param lastContiguousGID [in] The last contiguous GID.
    ///
    /// \return Struct containing the min and max GID on the calling
    ///   process, and the counts of failed and duplicate inserts into
    ///   the GID->LID table.  The count of failed inserts should
    ///   always be zero.  The count of duplicate inserts may be
    ///   nonzero, if the GID list on the calling process contains
    ///   duplicates.  This is generally a bad idea, however.
    template<class LO, class GO, class DeviceType>
    MapData<LO, GO, DeviceType>
    fillGlobalToLocalTable (Kokkos::UnorderedMap<GO, LO, DeviceType>& glMap,
                            const Kokkos::View<const GO*, DeviceType>& entries,
                            const GO firstContiguousGID,
                            const GO lastContiguousGID)
    {
      typedef typename DeviceType::size_type size_type;
      typedef GlobalToLocalTableFiller<LO, GO, DeviceType> functor_type;

      const size_type numEnt = entries.dimension_0 ();
      const size_type mapCap = glMap.capacity ();
      if (mapCap < numEnt) {
        Teuchos::ArrayView<const GO> lgMapAv (entries.ptr_on_device (), numEnt);
        TEUCHOS_TEST_FOR_EXCEPTION(
          mapCap < numEnt, std::logic_error, "fillGlobalToLocalTable: Input "
          "GID->LID table to fill does not have sufficient capacity for the "
          "given list of GIDs.  The table has capacity " << mapCap << ", but the "
          "input list of GIDs has " << numEnt << " entr" <<
          (numEnt != 1 ? "ies" : "y") << ".");
      }

      functor_type filler (glMap, entries, firstContiguousGID, lastContiguousGID);

      MapData<LO, GO, DeviceType> result;
      Kokkos::parallel_reduce (numEnt, filler, result);
      return result;
    }

    /// \class Map
    /// \brief Device implementation of Tpetra::Map
    /// \tparam LO The local ordinal type
    /// \tparam GO The global ordinal type
    /// \tparam DeviceType The Kokkos device type
    ///
    /// All methods marked \c const may be called in Kokkos parallel
    /// kernels.  No method that is not marked \c const may be called
    /// in a Kokkos parallel kernel.
    ///
    /// This class must <i>not</i> inherit from Teuchos::Describable,
    /// because that class contains data (e.g., the object label) that
    /// cannot live in CUDA device memory.  Every datum in this device
    /// must be able to live on a CUDA device, without requiring UVM.
    template<class LO, class GO, class DeviceType>
    class Map {
      template<class LO2, class GO2, class D2> friend class Map;
    public:
      //! \name Typedefs
      //@{

      //! The type of local indices.
      typedef LO local_ordinal_type;
      //! The type of global indices.
      typedef GO global_ordinal_type;
      //! The type of the Kokkos Device.
      typedef DeviceType execution_space;

      //@}
      //! \name Constructors
      //@{

      //! Empty constructor is helpful for clone(), etc.
      Map () :
        invalidGlobalIndex_ (Teuchos::OrdinalTraits<GO>::invalid ()), // final
        invalidLocalIndex_ (Teuchos::OrdinalTraits<LO>::invalid ()), // final
        globalNumIndices_ (0),
        myNumIndices_ (0),
        indexBase_ (0),
        firstContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        lastContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        minMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        maxMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        minAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        maxAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()),
        contiguous_ (false),
        distributed_ (true), // could be either false or true
        uniform_ (false),
        initialized_ (false)
      {}

      //! Contiguous uniform constructor.
      Map (const GO globalNumIndices,
           const GO indexBase,
           const Teuchos::Comm<int>& comm,
           const LocalGlobal lOrG) :
        invalidGlobalIndex_ (Teuchos::OrdinalTraits<GO>::invalid ()), // final
        invalidLocalIndex_ (Teuchos::OrdinalTraits<LO>::invalid ()), // final
        globalNumIndices_ (globalNumIndices), // final
        myNumIndices_ (0), // set below
        indexBase_ (indexBase), // final
        minMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // set below
        maxMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // set below
        minAllGID_ (globalNumIndices == 0 ?
                    invalidGlobalIndex_ :
                    indexBase), // final
        maxAllGID_ (globalNumIndices == 0 ?
                    invalidGlobalIndex_ :
                    indexBase + globalNumIndices - 1), // final
        contiguous_ (true), // final
        distributed_ (true), // this may be changed below
        uniform_ (true), // final
        initialized_ (true) // final
      {
        using Teuchos::as;
        using Teuchos::broadcast;
        using Teuchos::outArg;
        using Teuchos::reduceAll;
        using Teuchos::REDUCE_MIN;
        using Teuchos::REDUCE_MAX;
        using Teuchos::typeName;

#ifdef HAVE_TPETRA_DEBUG
        // In debug mode only, check whether globalNumIndices and
        // indexBase are the same over all processes in comm.
        {
          GO proc0NumGlobalIndices = globalNumIndices;
          broadcast<int, GO> (comm, 0, outArg (proc0NumGlobalIndices));
          GO minGlobalNumIndices = globalNumIndices;
          GO maxGlobalNumIndices = globalNumIndices;
          reduceAll<int, GO> (comm, REDUCE_MIN, globalNumIndices, outArg (minGlobalNumIndices));
          reduceAll<int, GO> (comm, REDUCE_MAX, globalNumIndices, outArg (maxGlobalNumIndices));
          TEUCHOS_TEST_FOR_EXCEPTION(
             minGlobalNumIndices != maxGlobalNumIndices || globalNumIndices != minGlobalNumIndices,
             std::invalid_argument,
             "Tpetra::Map constructor: All processes must provide the same number "
             "of global indices.  Process 0 set globalNumIndidces = "
             << proc0NumGlobalIndices << ".  The calling process "
             << comm.getRank () << " set globalNumIndices = " << globalNumIndices
             << ".  The min and max values over all processes are "
             << minGlobalNumIndices << " resp. " << maxGlobalNumIndices << ".");

          GO proc0IndexBase = indexBase;
          broadcast<int, GO> (comm, 0, outArg (proc0IndexBase));
          GO minIndexBase = indexBase;
          GO maxIndexBase = indexBase;
          reduceAll<int, GO> (comm, REDUCE_MIN, indexBase, outArg (minIndexBase));
          reduceAll<int, GO> (comm, REDUCE_MAX, indexBase, outArg (maxIndexBase));
          TEUCHOS_TEST_FOR_EXCEPTION(
            minIndexBase != maxIndexBase || indexBase != minIndexBase,
          std::invalid_argument,
          "Tpetra::Map constructor: "
          "All processes must provide the same indexBase argument.  "
          "Process 0 set indexBase = " << proc0IndexBase << ".  The calling "
          "process " << comm.getRank () << " set indexBase = " << indexBase
          << ".  The min and max values over all processes are "
          << minIndexBase << " resp. " << maxIndexBase << ".");
        }
#endif // HAVE_TPETRA_DEBUG

        // Distribute the GIDs (global indices) across the processes in
        // the given communicator so that they are
        //
        //   - Nonoverlapping (only one process owns each GID)
        //   - Contiguous (the sequence of GIDs is nondecreasing, and no
        //     two adjacent GIDs differ by more than one)
        //   - As evenly distributed as possible (the numbers of GIDs on
        //     two different processes do not differ by more than one)

        // All processes have the same globalNumIndices, but we still need
        // to check that it is valid.  globalNumIndices must be positive
        // and not the "invalid" value.
        //
        // This comparison looks funny, but it avoids compiler warnings
        // for comparing unsigned integers (GO might possibly be unsigned).
        TEUCHOS_TEST_FOR_EXCEPTION(
          (globalNumIndices < 1 && globalNumIndices != 0), std::invalid_argument,
          "Tpetra::Map constructor: globalNumIndices (= " << globalNumIndices
          << ") must be nonnegative.");

        TEUCHOS_TEST_FOR_EXCEPTION(
          globalNumIndices == static_cast<GO> (getInvalidGlobalIndex ()), std::invalid_argument,
          "Tpetra::Map constructor: You provided globalNumIndices = Teuchos::"
          "OrdinalTraits<GO>::invalid().  This version of the constructor "
          "requires a valid value of globalNumIndices.  You probably mistook "
          "this constructor for the \"contiguous nonuniform\" constructor, "
          "which can compute the global number of indices for you if you set "
          "globalNumIndices to that value.");

        if (lOrG == GloballyDistributed) {
          if (globalNumIndices == 0) {
            myNumIndices_ = 0;
            minMyGID_ = getInvalidGlobalIndex ();
            maxMyGID_ = getInvalidGlobalIndex ();
          }
          else { // globalNumIndices != 0, should be > 0
            // Compute myNumIndices:
            //
            // If globalNumIndices == numProcs * B + remainder,
            // then Proc r gets B+1 elements if r < remainder,
            // and B elements if r >= remainder.
            //
            // This strategy is valid for any value of globalNumIndices and
            // numProcs, including the following border cases:
            //   - numProcs == 1
            //   - myNumIndices < numProcs
            //
            // In the former case, remainder == 0 && globalNumIndices ==
            // myNumIndices.  In the latter case, remainder ==
            // globalNumIndices && myNumIndices is either 0 or 1.
            const GO numProcs = static_cast<GO> (comm.getSize ());
            const GO myRank = static_cast<GO> (comm.getRank ());
            const GO quotient  = globalNumIndices / numProcs;
            const GO remainder = globalNumIndices - quotient * numProcs;

            GO startIndex;
            if (myRank < remainder) {
              myNumIndices_ = 1 + quotient;
              // myRank was originally an int, so it should never overflow
              // reasonable GO types.
              startIndex = static_cast<GO> (myRank) * static_cast<GO> (myNumIndices_);
            } else {
              myNumIndices_ = static_cast<LO> (quotient);
              startIndex = static_cast<GO> (myRank) * static_cast<GO> (myNumIndices_) +
                static_cast<GO> (remainder);
            }
            minMyGID_ = indexBase + startIndex;
            maxMyGID_ = indexBase + startIndex + myNumIndices_ - 1;
          }
          distributed_ = (comm.getSize () > 1);
        }
        else {  // lOrG == LocallyReplicated
          myNumIndices_ = static_cast<LO> (globalNumIndices);
          if (globalNumIndices == 0) {
            minMyGID_ = getInvalidGlobalIndex ();
            maxMyGID_ = getInvalidGlobalIndex ();
          }
          else {
            minMyGID_ = indexBase;
            maxMyGID_ = indexBase + globalNumIndices - 1;
          }
          distributed_ = false;
        }

        firstContiguousGID_ = minMyGID_;
        lastContiguousGID_ = maxMyGID_;
      }

      //! Contiguous but possibly nonuniform constructor.
      Map (const GO globalNumIndices,
           const LO myNumIndices,
           const GO indexBase,
           const Teuchos::Comm<int>& comm) :
        invalidGlobalIndex_ (Teuchos::OrdinalTraits<GO>::invalid ()), // final
        invalidLocalIndex_ (Teuchos::OrdinalTraits<LO>::invalid ()), // final
        globalNumIndices_ (globalNumIndices), // provisional, if invalid()
        myNumIndices_ (myNumIndices), // final
        indexBase_ (indexBase), // final
        firstContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        lastContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        minMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // set below
        maxMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // set below
        minAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        maxAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        contiguous_ (true),  // final
        distributed_ (true), // set below
        uniform_ (false),    // final (conservative; we could try to detect this)
        initialized_ (true)  // final
      {
        using Teuchos::broadcast;
        using Teuchos::outArg;
        using Teuchos::reduceAll;
        using Teuchos::REDUCE_MIN;
        using Teuchos::REDUCE_MAX;
        using Teuchos::REDUCE_SUM;
        using Teuchos::scan;

#ifdef HAVE_TPETRA_DEBUG
        // Keep this for later debug checks.
        GO debugGlobalSum = 0; // Will be global sum of myNumIndices
        reduceAll<int, GO> (comm, REDUCE_SUM, static_cast<GO> (myNumIndices),
                            outArg (debugGlobalSum));
        // In debug mode only, check whether globalNumIndices and
        // indexBase are the same over all processes in comm.
        {
          GO proc0GlobalNumIndices = globalNumIndices;
          broadcast<int, GO> (comm, 0, outArg (proc0GlobalNumIndices));
          GO minGlobalNumIndices = globalNumIndices;
          GO maxGlobalNumIndices = globalNumIndices;
          reduceAll<int, GO> (comm, REDUCE_MIN, globalNumIndices, outArg (minGlobalNumIndices));
          reduceAll<int, GO> (comm, REDUCE_MAX, globalNumIndices, outArg (maxGlobalNumIndices));
          TEUCHOS_TEST_FOR_EXCEPTION(
            minGlobalNumIndices != maxGlobalNumIndices || globalNumIndices != minGlobalNumIndices,
            std::invalid_argument,
            "Tpetra::Map constructor: All processes must provide the same number "
            "of global indices.  This is true even if that argument is Teuchos::"
            "OrdinalTraits<global_size_t>::invalid() to signal that the Map should "
            "compute the global number of elements.  Process 0 set globalNumIndices"
            " = " << proc0GlobalNumIndices << ".  The calling process "
            << comm.getRank () << " set globalNumIndices = " << globalNumIndices
            << ".  The min and max values over all processes are "
            << minGlobalNumIndices << " resp. " << maxGlobalNumIndices << ".");

          GO proc0IndexBase = indexBase;
          broadcast<int, GO> (comm, 0, outArg (proc0IndexBase));
          GO minIndexBase = indexBase;
          GO maxIndexBase = indexBase;
          reduceAll<int, GO> (comm, REDUCE_MIN, indexBase, outArg (minIndexBase));
          reduceAll<int, GO> (comm, REDUCE_MAX, indexBase, outArg (maxIndexBase));
          TEUCHOS_TEST_FOR_EXCEPTION(
            minIndexBase != maxIndexBase || indexBase != minIndexBase,
            std::invalid_argument,
            "Tpetra::Map constructor: "
            "All processes must provide the same indexBase argument.  "
            "Process 0 set indexBase = " << proc0IndexBase << ".  The calling "
            "process " << comm.getRank () << " set indexBase = " << indexBase
            << ".  The min and max values over all processes are "
            << minIndexBase << " resp. " << maxIndexBase << ".");

          // Make sure that the sum of myNumIndices over all processes
          // equals globalNumIndices.
          TEUCHOS_TEST_FOR_EXCEPTION(
            globalNumIndices != getInvalidGlobalIndex () && debugGlobalSum != globalNumIndices,
            std::invalid_argument,
            "Tpetra::Map constructor: The sum of myNumIndices over all "
            "processes = " << debugGlobalSum << " != globalNumIndices = "
            << globalNumIndices << ".  If you would like this constructor to "
            "compute globalNumIndices for you, you may set globalNumIndices = "
            "Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid() on input.");
        }
#endif // HAVE_TPETRA_DEBUG

        // Distribute the GIDs (global indices) across the processes in
        // the given communicator so that they are
        //
        //   - Nonoverlapping (only one process owns each GID)
        //   - Contiguous (the sequence of GIDs is nondecreasing, and no
        //     two adjacent GIDs differ by more than one)
        //
        // This differs from the first Map constructor (that only takes a
        // global number of GIDs) in that the user has specified the
        // number of LIDs (local indices), so that the GIDs are not
        // (necessarily) evenly distributed over the processes.

        // Compute my local offset.  This is an inclusive scan, so to get
        // the final offset, we subtract off the input.
        GO scanResult = 0;
        scan<int, GO> (comm, REDUCE_SUM, myNumIndices, outArg (scanResult));
        const GO myOffset = scanResult - myNumIndices;

        if (globalNumIndices != static_cast<GO> (getInvalidGlobalIndex ())) {
          globalNumIndices_ = globalNumIndices; // Use the user's value.
        }
        else {
          // Inclusive scan means that the last process has the final sum.
          // Rather than doing a reduceAll to get the sum of
          // myNumIndices, we can just have the last process broadcast
          // its result.  That saves us a round of log(numProcs) messages.
          const int numProcs = comm.getSize ();
          GO globalSum = scanResult;
          if (numProcs > 1) {
            broadcast<int, GO> (comm, numProcs - 1, outArg (globalSum));
          }
          globalNumIndices_ = globalSum;

#ifdef HAVE_TPETRA_DEBUG
          // No need for an all-reduce here; both come from collectives.
          TEUCHOS_TEST_FOR_EXCEPTION(
            globalSum != debugGlobalSum, std::logic_error,
            "Tpetra::Map constructor (contiguous nonuniform): "
            "globalSum = " << globalSum << " != debugGlobalSum = " << debugGlobalSum
            << ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
        }
        // Now globalNumIndices_ is set.

        if (globalNumIndices_ != 0) {
          minMyGID_ = indexBase + myOffset;
          maxMyGID_ = indexBase + myOffset + myNumIndices - 1;
          firstContiguousGID_ = minMyGID_;
          lastContiguousGID_ = maxMyGID_;
          minAllGID_ = indexBase;
          maxAllGID_ = indexBase + globalNumIndices_ - 1;
        }

        distributed_ = checkIsDist (comm);
      }

      //! Possibly noncontiguous constructor.
      Map (const GO globalNumIndices,
           const Kokkos::View<const GO*, DeviceType>& myGlobalIndices,
           const GO indexBase,
           const Teuchos::Comm<int>& comm) :
        invalidGlobalIndex_ (Teuchos::OrdinalTraits<GO>::invalid ()), // final
        invalidLocalIndex_ (Teuchos::OrdinalTraits<LO>::invalid ()), // final
        globalNumIndices_ (globalNumIndices), // provisional, if invalid()
        myNumIndices_ (myGlobalIndices.dimension_0 ()), // final
        indexBase_ (indexBase), // final
        firstContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        lastContiguousGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        minMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        maxMyGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        minAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        maxAllGID_ (Teuchos::OrdinalTraits<GO>::invalid ()), // initialize below
        lgMap_ (myGlobalIndices), // final (input assumed correct)
        contiguous_ (false), // final (conservative; we could try to detect this)
        distributed_ (true), // set below
        uniform_ (false),    // final (conservative; we could try to detect this)
        initialized_ (true)  // final
      {
        using Teuchos::outArg;
        using Teuchos::reduceAll;

        // FIXME (mfh 20 Feb 2013, 05 Feb 2014) The global reduction
        // is redundant, since the directory Map will have to do the
        // same thing.  We should actually do the scan and broadcast
        // for the directory Map here, and give the computed offsets
        // to the directory Map's constructor.
        if (globalNumIndices_ == invalidGlobalIndex_) {
          // The user wants us to compute the sum.
          reduceAll<int, GO> (comm, Teuchos::REDUCE_SUM,
                              static_cast<GO> (myNumIndices_),
                              outArg (globalNumIndices_));
        } // now globalNumIndices_ is final.

        Kokkos::View<const GO*, DeviceType> nonContigEntries;
        if (myNumIndices_ > 0) {
          // Find contiguous GID range, with the restriction that the
          // beginning of the range starts with the first entry.

          // NOTE (mfh 05 Feb 2014) This assumes UVM, in that the
          // host here is reading from a device View.
          firstContiguousGID_ = lgMap_(0);
          //lastContiguousGID_ = firstContiguousGID_ + 1;

          // This is always true, as long as there is at least one
          // GID; that is, there is always at least one entry in the
          // sequence of contiguous GIDs, namely the first entry.
          lastContiguousGID_ = firstContiguousGID_;

          // By choosing LO so that it can count the local number of
          // elements, the user asserts that it can fit any values in
          // [0, myNumIndices_-1].  Thus, it's OK for i to be LO.
          GO i = 1;
          for ( ; i < myNumIndices_; ++i) {
            const GO curGid = lgMap_(i);
            const GO newLastContigGID = lastContiguousGID_ + 1;
            if (newLastContigGID != curGid) break;
            lastContiguousGID_ = newLastContigGID;
          }
          //--lastContiguousGID_;

          // [firstContiguousGID, lastContigousGID] is the initial
          // sequence of contiguous GIDs.  The sequence always
          // contains at least the first GID.

          const size_t numContig =
            static_cast<size_t> (lastContiguousGID_ - firstContiguousGID_ + 1);
          const size_t numNonContig =
            static_cast<size_t> (myNumIndices_) - numContig;

          // Make a view of the given GIDs, _not_ including the
          // initial sequence of contiguous GIDs.  The GID -> LID
          // table does not store this initial sequence.

          std::pair<size_t, size_t> offsetInfo (numContig, numContig + numNonContig);
          nonContigEntries =
            Kokkos::subview<Kokkos::View<const GO*, DeviceType> > (lgMap_, offsetInfo);

          TEUCHOS_TEST_FOR_EXCEPTION(
            static_cast<size_t> (nonContigEntries.dimension_0 ()) != numNonContig,
            std::logic_error, "Tpetra::Details::Map: On Process " << comm.getRank ()
            << ", computed offset into noncontiguous entries incorrectly.  "
            "nonContigEntries.dimension_0() = " << nonContigEntries.dimension_0 ()
            << ", but numNonContig = " << numNonContig << ".");

          TEUCHOS_TEST_FOR_EXCEPTION(
            numNonContig != 0 && nonContigEntries(0) != lgMap_(numContig),
            std::logic_error, "Tpetra::Details::Map: On Process " << comm.getRank ()
            << ", computed view of noncontiguous entries incorrectly.  "
            "nonContigEntries(0) = " << nonContigEntries(0) << " != "
            "lgMap_(numNonContig = " << numNonContig << ") = " <<
            lgMap_(numNonContig) << ".");

          //
          // Fill the GID -> LID table, and compute local min and max
          // GIDs, both at the same time in a single parallel kernel.
          //

          // Allocate space for GID -> LID table.  Leave extra space
          // to avoid excessive collisions.
          typedef typename DeviceType::size_type size_type;
          const size_type tableSize =
            static_cast<size_type> (1.25 * static_cast<double> (numNonContig));
          Kokkos::UnorderedMap<GO, LO, DeviceType> glMap (tableSize);
          MapData<LO, GO, DeviceType> result =
            fillGlobalToLocalTable (glMap, nonContigEntries,
                                    firstContiguousGID_,
                                    lastContiguousGID_);

          // There should be no failed inserts, since we made the
          // UnorderedMap more than big enough to hold all the GIDs.
          TEUCHOS_TEST_FOR_EXCEPTION(
            result.failedCount != 0, std::logic_error, "Tpetra::Details::Map: "
            << result.failedCount << " insert(s) into GID->LID table failed "
            "out of " << numNonContig << " total GIDs to insert.");

          minMyGID_ = result.minMyGID;
          maxMyGID_ = result.maxMyGID;

          // This is a shallow copy, that casts to const.  After
          // creating the UnorderedMap, we don't need to modify it, so
          // we store it as const.  This improves look-up performance
          // on NVIDIA GPUs, since lookups can go through (read-only)
          // texture cache.
          glMap_ = glMap;
        }
        else {
          // mfh 10 Feb 2014: A common case for users of Map is to
          // execute the following loop:
          //
          // for (GO g = map.getMinGlobalIndex (); g <= map.getMaxGlobalIndex (); ++g) { ... }
          //
          // Unfortunately, Tpetra chose to make this an inclusive
          // range.  This means that if the calling process owns no
          // GIDs, we need to make the above loop have zero
          // iterations.  We can't do this if both the min and max GID
          // on the calling process are the invalid GID.  We can if we
          // set the min GID to indexBase+1 and the max GID to
          // indexBase.
          minMyGID_ = indexBase + 1;
          maxMyGID_ = indexBase;

          // This insures tests for GIDs in the range [firstContiguousGID,
          // lastContiguousGID] fail for processes with no local elements.
          firstContiguousGID_ = indexBase + 1;
          lastContiguousGID_ = indexBase;
        }

        // mfh 20 Feb 2014: Commenting out the exception test below
        // fixes Bug 5401.  In particular, it needs to be OK to supply
        // -1 as the min GID, even though Teuchos::OrdinalTraits<GO>
        // for signed GO is -1 (which is actually annoying -- the most
        // negative integer (e.g., INT_MIN for GO = int) would be
        // better).
        //
        // TEUCHOS_TEST_FOR_EXCEPTION(
        //   minMyGID_ == getInvalidGlobalIndex () ||
        //   maxMyGID_ == getInvalidGlobalIndex (),
        //   std::logic_error, "Tpetra::Details::Map noncontig ctor: The calling "
        //   "process " << comm.getRank () << " owns " << myNumIndices_ << " (!= "
        //   "0) global indices, but was unable to find the min or max global "
        //   "index on this process.  Please report this bug to the Tpetra "
        //   "developers.  Tell them that minMyGID_ == " << minMyGID_ << " and "
        //   "maxMyGID_ == " << maxMyGID_ << ".");

        if (globalNumIndices_ == 0) {
          // mfh 10 Feb 2014: Users might want to execute the
          // following loop:
          //
          // for (GO g = map.getMinAllGlobalIndex (); g <= map.getMaxAllGlobalIndex (); ++g) { ... }
          //
          // Unfortunately, Tpetra chose to make this an inclusive
          // range.  This means that if the Map is empty (has no GIDs
          // on any process), we need to make the above loop have zero
          // iterations.  We can't do this if both the global min and
          // max GIDs are the invalid GID.  We can if we set the
          // global min GID to indexBase+1 and the global max GID to
          // indexBase.
          minAllGID_ = indexBase + 1;
          maxAllGID_ = indexBase;
          distributed_ = (comm.getSize () > 1);
        }
        else {
          // Compute the min and max of all processes' GIDs.  If
          // myNumIndices_ == 0 on this process, that makes things a bit
          // tricky.  Since Tpetra requires that indexBase be the min
          // GID over all processes, that makes things a bit easier.  If
          // a process owns no GIDs, we can therefore set its min and
          // max GID to indexBase in the all-reduce.  This will not
          // influence the final result.
          //
          // mfh 10 Feb 2014: If Tpetra decides to change to allow
          // indexBase not to be the min GID, that would break the above
          // approach.  One of the Tpetra tests already assumes this: it
          // uses an index base different than the actual global min
          // GID.  In that case, if GO is signed, and a process owns no
          // GIDs, we can do the following on that process:
          //
          //   - Set its min GID (in the all-reduce, not the value of
          //     minMyGID_) to std::numeric_limits<GO>::max(), and its
          //     max GID to std::numeric_limits<GO>::min().
          //   - While we're at it, use the same all-reduce to figure
          //     out if the Map is distributed.  "Distributed" means
          //     that there is at least one process with a number of
          //     local elements less than the number of global elements.
          //
          // This works if GO is signed because in that case,
          // std::numeric_limits<GO>::min() is the biggest negative
          // value, which users are unlikely to want to use as a GID.
          // If GO is unsigned, however, this won't work.  In that case,
          // std::numeric_limits<GO>::min() is zero, which is a typical
          // smallest GID.  There are two approaches we could take:
          //
          //   1. Find the max GID first, in a separate all-reduce.
          //      (On processes that own no GIDs, use 0 (the least
          //      unsigned integer) as the input of the all-reduce.)
          //      Then, on processes that own no GIDs, use that result
          //      (the global max GID) in place of
          //      std::numeric_limits<GO>::max(), when doing the
          //      all-reduce to find the global min GID.
          //
          //   2. Convert from GO to a signed type of the same size,
          //      and use the above approach for signed GO.
          //
          // Option 1 is better because Option 2 reduces the range of
          // allowed GID values by a factor of two.
          if (std::numeric_limits<GO>::is_signed) {
            // Compute the min and max of all processes' GIDs using a
            // single MAX all-reduce, using min(x,y) == -max(-x,-y).
            // This only works if x and y are signed.  If the calling
            // process owns no GIDs, make its min GID input to the
            // all-reduce -std::numeric_limits<GO>::min(), which should
            // be no greater than std::numeric_limits<GO>::max().
            //
            // If each process sets localDist=1 if its number of local
            // elements is strictly less than the number of global
            // elements, and localDist=0 otherwise, then a MAX
            // all-reduce on localDist tells us if the Map is
            // distributed (1 if yes, 0 if no).  Thus, we can append
            // localDist onto the end of the data and get the global
            // result from the all-reduce.

            // Does my process NOT own all the elements?
            const GO localDist =
              (static_cast<GO> (myNumIndices_) < globalNumIndices_) ? 1 : 0;

            GO minMaxInput[3];
            minMaxInput[0] = (myNumIndices_ == 0) ?
              -std::numeric_limits<GO>::min () : -minMyGID_;
            minMaxInput[1] = (myNumIndices_ == 0) ?
              std::numeric_limits<GO>::min () : maxMyGID_;
            minMaxInput[2] = localDist;

            GO minMaxOutput[3];
            minMaxOutput[0] = 0; // arbitrary initialization
            minMaxOutput[1] = 0; // arbitrary initialization
            minMaxOutput[2] = 0; // arbitrary initialization
            reduceAll<int, GO> (comm, Teuchos::REDUCE_MAX, 3,
                                minMaxInput, minMaxOutput);
            minAllGID_ = -minMaxOutput[0];
            maxAllGID_ = minMaxOutput[1];
            const GO globalDist = minMaxOutput[2];
            distributed_ = (comm.getSize () > 1 && globalDist == 1);
          }
          else { // GO is unsigned
            const GO maxInput = (myNumIndices_ == 0) ? 0 : maxMyGID_;
            reduceAll<int, GO> (comm, Teuchos::REDUCE_MAX, maxInput,
                                outArg (maxAllGID_));
            const GO minInput = (myNumIndices_ == 0) ? maxAllGID_ : minMyGID_;
            reduceAll<int, GO> (comm, Teuchos::REDUCE_MIN, minInput,
                                outArg (minAllGID_));

            // FIXME (mfh 10 Feb 2014) We could combine the MIN
            // all-reduce above with the MIN all-reduce below.
            if (comm.getSize () > 1) {
              // The communicator has more than one process, but that
              // doesn't necessarily mean the Map is distributed.
              int localRep = 0;
              if (getGlobalNumIndices () == static_cast<GO> (getMyNumIndices ())) {
                // The number of local elements on this process equals
                // the number of global elements.
                //
                // NOTE (mfh 22 Nov 2011) Does this still work if there
                // were duplicates in the global ID list on input (the
                // third Map constructor), so that the number of local
                // elements (which are not duplicated) on this process
                // could be less than the number of global elements,
                // even if this process owns all the elements?
                localRep = 1;
              }
              int allLocalRep;
              reduceAll<int, int> (comm, Teuchos::REDUCE_MIN, localRep,
                                   outArg (allLocalRep));
              if (allLocalRep != 1) {
                // At least one process does not own all the elements.
                // This makes the Map a distributed Map.
                distributed_ = true;
              }
            }
            else {
              // If the communicator has only one process, then the Map
              // is not distributed.
              distributed_ = false;
            }
          } // if GO is signed
        }

        // mfh 10 Feb 2014: tpetra/test/Map/Map_UnitTests.cpp,
        // indexBaseAndAllMin test uses a different indexBase input
        // (0) than the actual min GID (1).  Other tests seem to
        // depend on this as well.  That's why I commented out the
        // test below.

        // TEUCHOS_TEST_FOR_EXCEPTION(
        //   globalNumIndices_ != 0 && minAllGID_ != indexBase_,
        //   std::invalid_argument,
        //   "Tpetra::Details::Map constructor (noncontiguous): "
        //   "The Map has " << globalNumIndices_ << " > 0 global indices, but the "
        //   "min global index " << minAllGID_ << " over all process(es) does not "
        //   "equal the given indexBase " << indexBase_ << ".");
      }

      /// \brief Make this Map a copy of the input Map, but for a
      ///   (possibly) different device.
      template<class InDeviceType>
      void
      create_copy_view (const Map<LO, GO, InDeviceType>& map)
      {
        if (! map.initialized ()) {
          invalidGlobalIndex_ = Teuchos::OrdinalTraits<GO>::invalid ();
          invalidLocalIndex_ = Teuchos::OrdinalTraits<LO>::invalid ();
          globalNumIndices_ = 0;
          myNumIndices_ = 0;
          indexBase_ = 0;
          firstContiguousGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          lastContiguousGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          minMyGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          maxMyGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          minAllGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          maxAllGID_ = Teuchos::OrdinalTraits<GO>::invalid ();
          contiguous_ = false;
          distributed_ = true;
          uniform_ = false;
          initialized_ = false;
        }
        else {
          invalidGlobalIndex_ = map.invalidGlobalIndex_;
          invalidLocalIndex_ = map.invalidLocalIndex_;
          globalNumIndices_ = map.globalNumIndices_;
          myNumIndices_ = map.myNumIndices_;
          indexBase_ = map.indexBase_;
          firstContiguousGID_ = map.firstContiguousGID_;
          lastContiguousGID_ = map.lastContiguousGID_;
          minMyGID_ = map.minMyGID_;
          maxMyGID_ = map.maxMyGID_;
          minAllGID_ = map.minAllGID_;
          maxAllGID_ = map.maxAllGID_;

          if (! map.contiguous_) {
            if (map.myNumIndices_ != 0 && map.glMap_.size () != 0) {
              // The calling process owns one or more GIDs, and some of
              // these GIDs are not contiguous.
              Kokkos::UnorderedMap<GO, LO, DeviceType> glMap;
              glMap.create_copy_view (map.glMap_);
              glMap_ = glMap;
            }

            // It's OK for this to be dimension 0; that means it hasn't been initialized yet.
            Kokkos::View<GO*, DeviceType> lgMap ("LID->GID", map.lgMap_.dimension_0 ());
            if (map.lgMap_.dimension_0 () != 0) {
              Kokkos::deep_copy (lgMap, map.lgMap_);
            }
            lgMap_ = lgMap; // shallow copy and cast to const
          }

          contiguous_ = map.contiguous_;
          distributed_ = map.distributed_;
          uniform_ = map.uniform_;
          initialized_ = map.initialized_;
        }
      }

      /// \brief Set the \c distributed_ state.
      ///
      /// This method exists to avoid making Tpetra::Map a friend of this class.
      ///
      /// \warning This method only sets the host datum.  Any "mirror"
      ///   that is resident in device memory won't get this setting
      ///   automatically.
      void setDistributed (const bool distributed) {
        distributed_ = distributed;
      }

      /// \brief Whether this Map has been initialized.
      ///
      /// The default Map constructor creates an uninitialized Map.
      /// All other Map constructors create an initialized Map.  The
      /// create_copy_view method set the output Map (the
      /// <tt>*this</tt> argument) as initialized.
      ///
      /// It is an error to use an uninitialized Map.  Results are
      /// undefined, but the resulting Map will likely behave as if it
      /// were empty.
      bool initialized () const {
        return initialized_;
      }

      //@}
      //! \name Methods safe to call in a Kokkos parallel kernel.
      //@{

      //! The number of indices owned by all processes in the Map's communicator.
      KOKKOS_INLINE_FUNCTION GO getGlobalNumIndices () const {
        return globalNumIndices_;
      }

      //! The number of indices owned by the calling process.
      KOKKOS_INLINE_FUNCTION LO getMyNumIndices () const {
        return myNumIndices_;
      }

      //! The index base for this Map.
      KOKKOS_INLINE_FUNCTION GO getIndexBase () const {
        return indexBase_;
      }

      /// \brief An invalid local index.
      ///
      /// Map's methods that return a local index use this index to
      /// indicate an error condition.  For example, if
      /// getLocalIndex() gets a global index that is not owned by the
      /// calling process, it returns an invalid local index, which
      /// equals the return value of this method.
      KOKKOS_INLINE_FUNCTION LO getInvalidLocalIndex () const {
        return invalidLocalIndex_;
      }

      //! The minimum local index (on the calling process).
      KOKKOS_INLINE_FUNCTION LO getMinLocalIndex () const {
        return static_cast<LO> (0);
      }

      /// \brief The maximum local index on the calling process.
      ///
      /// If this process owns no indices, that is, if
      /// <tt>getMyNumIndices() == 0</tt>, then this method returns
      /// <tt>Teuchos::OrdinalTraits<LO>::invalid()</tt>.
      KOKKOS_INLINE_FUNCTION LO getMaxLocalIndex () const {
        const LO myNumIndices = getMyNumIndices ();
        // Local indices are always zero-based.
        return (myNumIndices == 0) ? getInvalidLocalIndex () : (myNumIndices - 1);
      }

      /// \brief An invalid global index.
      ///
      /// Map's methods that return a global index use this index to
      /// indicate an error condition.  For example, if
      /// getGlobalIndex() gets a local index that is not owned by the
      /// calling process, it returns an invalid global index, which
      /// equals the return value of this method.
      KOKKOS_INLINE_FUNCTION LO getInvalidGlobalIndex () const {
        return invalidGlobalIndex_;
      }

      //! The minimum global index owned by the calling process.
      KOKKOS_INLINE_FUNCTION GO getMinGlobalIndex () const {
        return minMyGID_;
      }

      //! The maximum global index owned by the calling process.
      KOKKOS_INLINE_FUNCTION GO getMaxGlobalIndex () const {
        return maxMyGID_;
      }

      //! The minimum global index over all processes in the communicator.
      KOKKOS_INLINE_FUNCTION GO getMinAllGlobalIndex () const {
        return minAllGID_;
      }

      //! The maximum global index over all processes in the communicator.
      KOKKOS_INLINE_FUNCTION GO getMaxAllGlobalIndex () const {
        return maxAllGID_;
      }

      /// \brief The local index corresponding to the given global index.
      ///
      /// If the given global index is not owned by this process, return
      /// Teuchos::OrdinalTraits<LO>::invalid().
      KOKKOS_INLINE_FUNCTION LO getLocalIndex (const GO globalIndex) const {
        if (getMyNumIndices() == 0) return getInvalidLocalIndex();
        if (isContiguous ()) {
          if (globalIndex >= getMinGlobalIndex () &&
              globalIndex <= getMaxGlobalIndex ()) {
            return static_cast<LO> (globalIndex - getMinGlobalIndex ());
          }
          else {
            return getInvalidLocalIndex ();
          }
        }
        else if (globalIndex >= firstContiguousGID_ &&
                 globalIndex <= lastContiguousGID_) {
          return static_cast<LO> (globalIndex - firstContiguousGID_);
        }
        else {
          // This also works if this process owns no indices.
          // In that case, glMap_ will be empty, so this branch
          // will always return the invalid LID.
          const typename global_to_local_table_type::size_type i =
            glMap_.find (globalIndex);
          return glMap_.valid_at (i) ? // if the GID is in the map, ...
            glMap_.value_at (i) : // ... return its LID, ...
            getInvalidLocalIndex (); // ... else return the invalid LID.
        }
      }

      /// \brief The global index corresponding to the given local index.
      ///
      /// If the given local index is not valid on the calling process,
      /// return Teuchos::OrdinalTraits<GO>::invalid().
      KOKKOS_INLINE_FUNCTION GO getGlobalIndex (const LO localIndex) const {
        if (getMyNumIndices() == 0) return getInvalidLocalIndex();
        if (localIndex < getMinLocalIndex () ||
            localIndex > getMaxLocalIndex ()) {
          // This process will always take this branch if it owns no indices.
          return getInvalidGlobalIndex ();
        }
        else if (isContiguous ()) {
          return getMinGlobalIndex () + localIndex;
        }
        else {
          return lgMap_(localIndex);
        }
      }

      //! Whether the given local index is owned by (valid on) the calling process.
      KOKKOS_INLINE_FUNCTION bool
      isOwnedLocalIndex (const LO localIndex) const {
        return localIndex >= getMinLocalIndex () &&
          localIndex <= getMaxLocalIndex ();
      }

      //! Whether the given global index is owned by the calling process.
      KOKKOS_INLINE_FUNCTION bool
      isOwnedGlobalIndex (const GO globalIndex) const {
        return this->getLocalIndex (globalIndex) != getInvalidLocalIndex ();
      }

      /// \brief Whether the range of global indices is uniform.
      ///
      /// This is a conservative quantity.  It need only be true if
      /// the Map was constructed using the first (uniform contiguous)
      /// constructor or a nonmember constructor that calls it.  We
      /// reserve the right to do more work to check this in the
      /// future.
      KOKKOS_INLINE_FUNCTION bool isUniform () const {
        return uniform_;
      }

      /// \brief True if this Map is distributed contiguously, else false.
      ///
      /// Currently, creating this Map using the constructor for a
      /// user-defined arbitrary distribution (that takes a list of
      /// global elements owned on each process) means that this
      /// method always returns false.  We currently make no effort to
      /// test whether the user-provided global indices are actually
      /// contiguous on all the processes.  Many operations may be
      /// faster for contiguous Maps.  Thus, if you know the indices
      /// are contiguous on all processes, you should consider using
      /// one of the constructors for contiguous elements.
      KOKKOS_INLINE_FUNCTION bool isContiguous () const {
        return contiguous_;
      }

      /// \brief Whether this Map is globally distributed or locally replicated.
      ///
      /// \return True if this Map is globally distributed, else false.
      ///
      /// "Globally distributed" means that <i>all</i> of the following
      /// are true:
      /// <ol>
      /// <li> The map's communicator has more than one process.</li>
      /// <li> There is at least one process in the map's
      ///      communicator, whose local number of indices does not
      ///      equal the number of global indices.  (That is, not all
      ///      the indices are replicated over all the processes.)
      /// </li>
      /// </ol>
      ///
      /// If at least one of the above are not true, then the map is
      /// "locally replicated."  (The two are mutually exclusive.)
      ///
      /// Calling this method requires no communication or
      /// computation, because the result is precomputed in Map's
      /// constructors.
      KOKKOS_INLINE_FUNCTION bool isDistributed () const {
        return distributed_;
      }

      //@}
      //! \name Methods <i>not</i> safe to call in a Kokkos parallel kernel
      //@{

      /// \brief Return a view of the global indices owned by this process.
      ///
      /// If you call this method on a contiguous Map, it will create
      /// and cache the list of global indices for later use.  Beware
      /// of calling this if the calling process owns a very large
      /// number of global indices.
      Kokkos::View<const GO*, DeviceType> getMyGlobalIndices () {
        const GO myNumIndices = getMyNumIndices ();

        if (myNumIndices != 0 && lgMap_.dimension_0 () == 0) {
          Kokkos::View<GO*, DeviceType> lgMap ("LID->GID", myNumIndices);
          // fill with [getMinGlobalIndex(), getMinGlobalIndex()+myNumIndices-1].
          iota<GO, DeviceType> (lgMap, this->getMinGlobalIndex ());
          lgMap_ = lgMap; // shallow copy; cast to const
        }
        return lgMap_;
      }

      //! A one-line description of this object.
      std::string description () {
        return "\"Tpetra::Details::Map\"";
      }

      //! Print this object with the given verbosity level to the given Teuchos::FancyOStream.
      void
      describe (Teuchos::FancyOStream &out,
                const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default)
      {
        using std::endl;
        const Teuchos::EVerbosityLevel vl =
          (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

        if (vl != Teuchos::VERB_NONE) {
          // By convention, describe() always starts with a tab.
          Teuchos::OSTab tab0 (out);
          out << description () << endl; // TODO (mfh 05 Feb 2014) Print more info.
        }
      }

    private:
      void validateLocally (const int myRank) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          invalidGlobalIndex_ != Teuchos::OrdinalTraits<GO>::invalid (),
          std::logic_error,
          "Tpetra::Details::Map: On Process " << myRank << ", "
          "invalidGlobalIndex_ was not set.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          invalidLocalIndex_ != Teuchos::OrdinalTraits<LO>::invalid (),
          std::logic_error,
          "Tpetra::Details::Map: On Process " << myRank << ", "
          "invalidLocalIndex_ was not set.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! contiguous_ && myNumIndices_ != lgMap_.dimension_0 (),
          std::logic_error, "Tpetra::Details::Map: On Process " << myRank <<
          ", this noncontiguous Map has myNumIndices_ = " << myNumIndices_ <<
          " != lgMap_.dimension_0 () = " << lgMap_.dimension_0 () << ".");

        // Test whether the Map contains all the noncontiguous indices
        // it claims to contain.
        if (! contiguous_ && myNumIndices_ > 0) {
          bool bad = false;
          Teuchos::Array<GO> badGlobals;
          Teuchos::Array<LO> badLocalsIn;
          Teuchos::Array<LO> badLocalsOut;

          for (LO localIndexIn = 0; localIndexIn < myNumIndices_; ++localIndexIn) {
            const GO globalIndex = this->getGlobalIndex (localIndexIn);
            const LO localIndexOut = this->getLocalIndex (globalIndex);

            if (localIndexIn != localIndexOut) {
              bad = true;
              badGlobals.push_back (globalIndex);
              badLocalsIn.push_back (localIndexIn);
              badLocalsOut.push_back (localIndexOut);
            }
          }

          // Get the GID->LID table's keys and values as Teuchos::Array.
          typedef typename DeviceType::size_type size_type;
          const size_type mapSize = glMap_.size ();
          Teuchos::Array<GO> keys;
          Teuchos::Array<LO> vals;
          keys.reserve (mapSize);
          vals.reserve (mapSize);
          for (size_type k = 0; k < mapSize; ++k) {
            keys.push_back (glMap_.key_at (k));
            vals.push_back (glMap_.value_at (k));
          }
          // Get the LID->GID table as an ArrayView for easy printing.
          Teuchos::ArrayView<const GO> lgMapAv (lgMap_.ptr_on_device (),
                                                lgMap_.dimension_0 ());
          TEUCHOS_TEST_FOR_EXCEPTION(
            bad, std::logic_error, "Tpetra::Details::Map: On Process " <<
            myRank << ", this noncontiguous Map's GID->LID lookup is broken.  "
            "The following GIDs are in this process' GID list, but attempting "
            "to look up their LIDs gives the wrong answer: " << badGlobals () <<
            ".  The corresponding LID list should have been " << badLocalsIn ()
            << ", but was " << badLocalsOut << " instead.  The Map's GID list "
            "is " << lgMapAv << ".  The keys in the GID->LID table are " <<
            keys () << " and their values are " << vals () << ".  The GID->LID "
            "table reports its size as " << mapSize << " and its capacity as "
            << glMap_.capacity () << ".  The initial sequence of contiguous "
            "GIDs on this process is [" << firstContiguousGID_ << "," <<
            lastContiguousGID_ << "].");
        }
      }

      /// \brief Whether this Map is globally distributed or locally replicated.
      ///
      /// \param comm [in] This Map's communicator.
      ///
      /// \return True if this Map is globally distributed, else false.
      ///
      /// This is a collective operation over the input communicator.
      /// See the documentation of \c isDistributed() for definitions
      /// of "globally distributed" and "locally replicated."
      ///
      /// Map invokes this method in its constructors if necessary, to
      /// set the distributed_ flag (and thus the return value of
      /// isDistributed()).  Map doesn't need to call checkIsDist()
      /// when using the uniform contiguous constructor with lg =
      /// GloballyDistributed, since then checking the number of
      /// processes in the communicator suffices.
      bool checkIsDist (const Teuchos::Comm<int>& comm) {
        using Teuchos::outArg;
        using Teuchos::REDUCE_MIN;
        using Teuchos::reduceAll;

        bool global = false;
        if (comm.getSize () > 1) {
          // The communicator has more than one process, but that
          // doesn't necessarily mean the Map is distributed.
          int localRep = 0;
          if (getGlobalNumIndices () == static_cast<GO> (getMyNumIndices ())) {
            // The number of local indices on this process equals the
            // number of global indices.
            //
            // NOTE (mfh 22 Nov 2011) Does this still work if there
            // were duplicates in the GID list on input (the third Map
            // constructor), so that the number of LIDs (which are not
            // duplicated) on this process could be less than the
            // number of GIDs, even if this process owns all the IDs?
            localRep = 1;
          }
          int allLocalRep;
          reduceAll<int, int> (comm, REDUCE_MIN, localRep, outArg (allLocalRep));
          if (allLocalRep != 1) {
            // At least one process does not own all the elements.
            // This makes the Map a distributed Map.
            global = true;
          }
        }
        // If the communicator has only one process, then the Map is
        // not distributed.
        return global;
      }

      /// \brief Same as Teuchos::OrdinalTraits<GlobalOrdinal>::invalid().
      ///
      /// This exists so that Teuchos::OrdinalTraits::invalid() does
      /// not need to work on the device.
      GO invalidGlobalIndex_;

      /// \brief Same as Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
      ///
      /// This exists so that Teuchos::OrdinalTraits::invalid() does
      /// not need to work on the device.
      LO invalidLocalIndex_;

      //! Number of GIDs in this Map over all processes in the Map's communicator.
      GO globalNumIndices_;

      //! Number of GIDs owned by the calling process.  This is always nonnegative.
      GO myNumIndices_;

      //! The index base for global IDs in this Map.
      GO indexBase_;

      //! First GID in the initial range of zero or more contiguous GIDs on the calling process.
      GO firstContiguousGID_;

      //! Last GID in the initial range of zero or more contiguous GIDs on the calling process.
      GO lastContiguousGID_;

      /// \brief Minimum GID on the calling process.
      ///
      /// If the calling process owns no GIDs, this is ???.
      GO minMyGID_;

      /// \brief Maximum GID on the calling process.
      ///
      /// If the calling process owns no GIDs, this is ???.
      GO maxMyGID_;

      /// \brief Minimum global ID in this Map over all processes in the
      ///   Map's communicator.
      GO minAllGID_;

      /// \brief Maximum global ID in this Map over all processes in the
      ///   Map's communicator.
      GO maxAllGID_;

      /// \brief Type of the table that maps global IDs to local IDs.
      ///
      /// This is a const table.  Neither keys nor values may be
      /// changed.  This means that when Map builds this table, it
      /// first has to build a nonconst table, and then assign the
      /// result to this table.
      typedef Kokkos::UnorderedMap<const GO, const LO, DeviceType> global_to_local_table_type;

      /// \brief A mapping from global IDs to local IDs.
      ///
      /// This is a local mapping.  Directory implements the global
      /// mapping for all global indices (both remote and locally
      /// owned).  This object corresponds roughly to
      /// Epetra_BlockMapData's LIDHash_ hash table (which also maps
      /// from global to local indices).
      ///
      /// This mapping is built only for a noncontiguous map, by the
      /// noncontiguous map constructor.  For noncontiguous maps, the
      /// getLocalElement() and isNodeGlobalElement() methods use this
      /// mapping.
      ///
      /// TODO (mfh 03 Feb 2014) By default, UnorderedMap uses the
      /// Murmur hash.  Epetra's hash function is faster for common
      /// cases.
      global_to_local_table_type glMap_;

      /// \brief A mapping from local IDs to global IDs.
      ///
      /// By definition, this mapping is local; it only contains GIDs
      /// owned by this process.  This mapping is created in two
      /// cases:
      /// <ol>
      /// <li> It is always created for a noncontiguous Map, in the
      ///    noncontiguous version of the Map constructor.</li>
      /// <li> In getNodeElementList(), on demand (if it wasn't created
      ///    before).</li>
      /// </ol>
      Kokkos::View<const GO*, DeviceType> lgMap_;

      /// \brief Whether the Map was constructed using one of the
      ///   contiguous constructors.
      ///
      /// The contiguous constructors define the global range of
      /// global indices to be contiguous and ordered.
      bool contiguous_;

      /// \brief Whether this Map's GIDs are non-identically
      ///   distributed among different processes.
      bool distributed_;

      /// \brief Whether the range of global indices is uniform.
      ///
      /// This is only true if the Map was constructed using the first
      /// (uniform contiguous) constructor or a nonmember constructor
      /// that calls it.
      bool uniform_;

      /// \brief Whether this Map has been initialized.
      ///
      /// Please see the documentation of initialized() for details.
      bool initialized_;
    };

    /// \class MapMirrorer
    /// \brief Return a mirror of either deviceMap or hostMap for a
    ///   possibly different output device type.
    ///
    /// \tparam OutMapType A specialization of Map; the type of the
    ///   return value of mirror().
    /// \tparam DeviceMapType A specialization of Map; the type of the
    ///   \c deviceMap input.
    /// \tparam HostMapType A specialization of Map; the type of the
    ///   \c hostMap input.
    ///
    /// MapMirrorer's mirror() method returns a "mirror copy" of
    /// either deviceMap or hostMap.  It assumes that both Maps
    /// represents the same distribution, just stored on possibly
    /// different devices.  A "mirror copy" is a shallow copy if
    /// possible, else a deep copy.  The method will favor the
    /// initialized Map, if one of deviceMap or hostMap is not
    /// initialized.
    ///
    /// The mirror() method will only compile if all three Map types
    /// have the same first two template parameters (LO and GO).
    ///
    /// If mirror() had to do work (a deep copy) to create the mirror,
    /// then we cache that work.  This is why the input arguments are
    /// passed in as nonconst references.
    template<class OutMapType, class DeviceMapType, class HostMapType>
    struct MapMirrorer {
      typedef OutMapType output_map_type;
      typedef DeviceMapType device_map_type;
      typedef HostMapType host_map_type;

      //! Return a mirror copy of either deviceMap or hostMap.
      static output_map_type
      mirror (device_map_type& deviceMap, host_map_type& hostMap);
    };

    // Partial specialization for when all three Maps have different
    // device types.
    template<class LO, class GO, class OutDeviceType, class DeviceType, class HostDeviceType>
    struct MapMirrorer<Map<LO, GO, OutDeviceType>, Map<LO, GO, DeviceType>, Map<LO, GO, HostDeviceType> > {
      typedef Map<LO, GO, OutDeviceType> output_map_type;
      typedef Map<LO, GO, DeviceType> device_map_type;
      typedef Map<LO, GO, HostDeviceType> host_map_type;

      static output_map_type
      mirror (device_map_type& deviceMap, host_map_type& hostMap) {
        output_map_type outputMap;

        // Favor the host Map, since if OutDeviceType != DeviceType,
        // then OutDeviceType is likely not to be a CUDA device (e.g.,
        // DeviceType = Cuda, HostDeviceType = OpenMP, and
        // OutDeviceType = Serial).
        if (hostMap.initialized ()) {
          outputMap.template create_copy_view<HostDeviceType> (hostMap);
        } else {
          outputMap.template create_copy_view<DeviceType> (deviceMap);
        }
        return outputMap;
      }
    };

    // Partial specialization for when the device and host Maps have
    // different device types, and the output Map has the same device
    // type as the device Map.
    template<class LO, class GO, class DeviceType, class HostDeviceType>
    struct MapMirrorer<Map<LO, GO, DeviceType>, Map<LO, GO, DeviceType>, Map<LO, GO, HostDeviceType> > {
      typedef Map<LO, GO, DeviceType> output_map_type;
      typedef Map<LO, GO, DeviceType> device_map_type;
      typedef Map<LO, GO, HostDeviceType> host_map_type;

      static output_map_type
      mirror (device_map_type& deviceMap, host_map_type& hostMap) {
        if (deviceMap.initialized ()) {
          return deviceMap;
        }
        else {
          output_map_type outMap;
          outMap.template create_copy_view<HostDeviceType> (hostMap);
          if (outMap.initialized () && ! deviceMap.initialized ()) {
            deviceMap = outMap; // cache the deep copy
          }
          return outMap;
        }
      }
    };

    // Partial specialization for when the device and host Maps have
    // different device types, and the output Map has the same device
    // type as the host Map.
    template<class LO, class GO, class DeviceType, class HostDeviceType>
    struct MapMirrorer<Map<LO, GO, HostDeviceType>, Map<LO, GO, DeviceType>, Map<LO, GO, HostDeviceType> > {
      typedef Map<LO, GO, HostDeviceType> output_map_type;
      typedef Map<LO, GO, DeviceType> device_map_type;
      typedef Map<LO, GO, HostDeviceType> host_map_type;

      static output_map_type
      mirror (device_map_type& deviceMap, host_map_type& hostMap) {
        if (hostMap.initialized ()) {
          return hostMap;
        }
        else {
          output_map_type outMap;
          outMap.template create_copy_view<DeviceType> (deviceMap);
          if (outMap.initialized ()) {
            hostMap = outMap; // cache the deep copy
          }
          return outMap;
        }
      }
    };

    // Partial specialization for when the device, host, and output
    // Maps have the same device types.
    template<class LO, class GO, class DeviceType>
    struct MapMirrorer<Map<LO, GO, DeviceType>, Map<LO, GO, DeviceType>, Map<LO, GO, DeviceType> > {
      typedef Map<LO, GO, DeviceType> output_map_type;
      typedef Map<LO, GO, DeviceType> device_map_type;
      typedef Map<LO, GO, DeviceType> host_map_type;

      static output_map_type
      mirror (device_map_type& deviceMap, const host_map_type& /* hostMap */ ) {
        return deviceMap;
      }
    };

  }
}

#endif // TPETRA_KOKKOSREFACTOR_DETAILS_MAP_HPP
