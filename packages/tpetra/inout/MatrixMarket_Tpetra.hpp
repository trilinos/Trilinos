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

#ifndef __MatrixMarket_Tpetra_hpp
#define __MatrixMarket_Tpetra_hpp

/// \file MatrixMarket_Tpetra.hpp
/// \brief Matrix Market file readers and writers for Tpetra objects.
/// \author Mark Hoemmen
///
/// This header file implements Matrix Market file readers and writers
/// for both sparse and dense matrices (as \c Tpetra::CrsMatrix
/// resp. \c Tpetra::MultiVector).  The Matrix Market (see their <a
/// href="http://math.nist.gov/MatrixMarket"> web site </a> for
/// details) defines a human-readable ASCII text file format ("Matrix
/// Market format") for interchange of sparse and dense matrices.
///
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_MatrixMarket_Raw_Adder.hpp"
#include "Teuchos_MatrixMarket_SymmetrizingAdder.hpp"
#include "Teuchos_MatrixMarket_assignScalar.hpp"
#include "Teuchos_MatrixMarket_Banner.hpp"
#include "Teuchos_MatrixMarket_CoordDataReader.hpp"
#include "Teuchos_MatrixMarket_SetScientific.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>

namespace Tpetra {
  ///
  /// \namespace MatrixMarket
  /// \brief Matrix Market file readers and writers for sparse and
  ///   dense matrices (as \c CrsMatrix resp. \c MultiVector).
  /// \author Mark Hoemmen
  ///
  /// The Matrix Market (see their <a
  /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
  /// details) defines a human-readable ASCII text file format
  /// ("Matrix Market format") for interchange of sparse and dense
  /// matrices.  This namespace defines classes for reading and
  /// writing sparse or dense matrices from a Matrix Market file or
  /// input stream.
  ///
  /// Matrix Market files are designed for easy reading and writing of
  /// test matrices by both humans and computers.  They are not
  /// intended for high-performance or parallel file input and output.
  /// You should use a true parallel file format if you want to do
  /// parallel input and output of sparse or dense matrices.  Since
  /// the Matrix Market format is not optimized for performance or
  /// parallelism, our readers and writers assume that the entire
  /// matrix can fit in a single MPI process.  We do all the file
  /// input or output on the MPI process with rank 0 ("Process 0").
  /// Distributed input matrices are gathered from all MPI processes
  /// in the participating communicator, and distributed output
  /// matrices are broadcast from Process 0 to all MPI processes in
  /// the participating communicator.
  ///
  namespace MatrixMarket {

    namespace {
#ifdef HAVE_MPI
      // We're communicating integers of type IntType.  Figure
      // out the right MPI_Datatype for IntType.  Usually it
      // is int or long, so these are good enough for now.
      template<class IntType> MPI_Datatype getMpiDatatype () {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented for IntType != int or long");
      }
      template<> MPI_Datatype getMpiDatatype<int> () { return MPI_INT; }
      template<> MPI_Datatype getMpiDatatype<long> () { return MPI_LONG; }
#endif // HAVE_MPI

      template<class IntType>
      void
      reduceSum (const IntType sendbuf[], 
		 IntType recvbuf[],
		 const int count,
		 const int root,
		 const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
#ifdef HAVE_MPI
        using Teuchos::RCP;
        using Teuchos::rcp_dynamic_cast;
        using Teuchos::MpiComm;

        // Get the raw MPI communicator, so we don't have to wrap MPI_Reduce in Teuchos.
        RCP<const MpiComm<int> > mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm);
        if (! mpiComm.is_null ()) {
          MPI_Datatype rawMpiType = getMpiDatatype<IntType> ();
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err = 
	    MPI_Reduce (reinterpret_cast<void*> (const_cast<IntType*> (sendbuf)),
			reinterpret_cast<void*> (const_cast<IntType*> (recvbuf)),
			count, rawMpiType, MPI_SUM, root, rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Reduce failed");
          return;
        }
#endif // HAVE_MPI
        // Serial communicator case: just copy.  count is the
        // amount to receive, so it's the amount to copy.
        std::copy (sendbuf, sendbuf + count, recvbuf);
      }		 

      template<class IntType>
      void
      gather (const IntType sendbuf[], const int sendcnt,
              IntType recvbuf[], const int recvcnt,
              const int root, const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
        using Teuchos::RCP;
        using Teuchos::rcp_dynamic_cast;
#ifdef HAVE_MPI
        using Teuchos::MpiComm;

        // Get the raw MPI communicator, so we don't have to wrap MPI_Gather in Teuchos.
        RCP<const MpiComm<int> > mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm);
        if (! mpiComm.is_null ()) {
          MPI_Datatype sendtype = getMpiDatatype<IntType> ();
          MPI_Datatype recvtype = sendtype;
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err = MPI_Gather (reinterpret_cast<void*> (const_cast<IntType*> (sendbuf)),
                                      sendcnt,
                                      sendtype,
                                      reinterpret_cast<void*> (recvbuf),
                                      recvcnt,
                                      recvtype,
                                      root,
                                      rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(
            err != MPI_SUCCESS, std::runtime_error, "MPI_Gather failed");
          return;
        }
#endif // HAVE_MPI

        TEUCHOS_TEST_FOR_EXCEPTION(
          recvcnt > sendcnt, std::invalid_argument, 
	  "gather: If the input communicator contains only one process, "
	  "then you cannot receive more entries than you send.  "
	  "You aim to receive " << recvcnt << " entries, "
	  "but to send " << sendcnt << ".");
        // Serial communicator case: just copy.  recvcnt is the
        // amount to receive, so it's the amount to copy.
        std::copy (sendbuf, sendbuf + recvcnt, recvbuf);
      }

      template<class IntType>
      void
      gatherv (const IntType sendbuf[], const int sendcnt,
               IntType recvbuf[], const int recvcnts[], const int displs[],
               const int root, const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
        using Teuchos::RCP;
        using Teuchos::rcp_dynamic_cast;
#ifdef HAVE_MPI
        using Teuchos::MpiComm;

        // Get the raw MPI communicator, so we don't have to wrap
        // MPI_Gather in Teuchos.
        RCP<const MpiComm<int> > mpiComm = 
	  rcp_dynamic_cast<const MpiComm<int> > (comm);
        if (! mpiComm.is_null ()) {
          MPI_Datatype sendtype = getMpiDatatype<IntType> ();
          MPI_Datatype recvtype = sendtype;
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err = MPI_Gatherv (reinterpret_cast<void*> (const_cast<IntType*> (sendbuf)),
                                       sendcnt,
                                       sendtype,
                                       reinterpret_cast<void*> (recvbuf),
                                       const_cast<int*> (recvcnts),
                                       const_cast<int*> (displs),
                                       recvtype,
                                       root,
                                       rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(
            err != MPI_SUCCESS, std::runtime_error, "MPI_Gatherv failed");
          return;
        }
#endif // HAVE_MPI
        TEUCHOS_TEST_FOR_EXCEPTION(
          recvcnts[0] > sendcnt, 
	  std::invalid_argument, 
	  "gatherv: If the input communicator contains only one process, "
	  "then you cannot receive more entries than you send.  "
	  "You aim to receive " << recvcnts[0] << " entries, "
	  "but to send " << sendcnt << ".");
        // Serial communicator case: just copy.  recvcnts[0] is the
        // amount to receive, so it's the amount to copy.  Start
        // writing to recvbuf at the offset displs[0].
        std::copy (sendbuf, sendbuf + recvcnts[0], recvbuf + displs[0]);
      }

      // Given an arbitrary Map, compute a Map containing all the GIDs
      // in the same order as in (the one-to-one version of) map, but
      // all owned exclusively by Proc 0.
      template<class MapType>
      Teuchos::RCP<const MapType>
      computeGatherMap (Teuchos::RCP<const MapType> map)
      {
        using Tpetra::createOneToOne;
	using Tpetra::global_size_t;
        using Teuchos::Array;
        using Teuchos::ArrayRCP;
        using Teuchos::ArrayView;
        using Teuchos::as;
        using Teuchos::Comm;
        using Teuchos::RCP;
	using std::cerr;
	using std::endl;
        typedef typename MapType::local_ordinal_type LO;
        typedef typename MapType::global_ordinal_type GO;
        typedef typename MapType::node_type NT;

        RCP<const Comm<int> > comm = map->getComm ();
        const int numProcs = comm->getSize ();
        const int myRank = comm->getRank ();
	// mfh 12 Feb 2013: This should always be false, unless debugging.
	const bool debug = false;

	if (debug && myRank == 0) {
	  cerr << "computeGatherMap:" << endl
	       << "- Calling createOneToOne" << endl;
	}

        RCP<const MapType> oneToOneMap;
        if (map->isContiguous ()) {
          oneToOneMap = map; // contiguous Maps are always 1-to-1
        } else {
	  // It could be that Map is one-to-one, but the class doesn't
	  // give us a way to test this, other than to create the
	  // one-to-one Map.
          oneToOneMap = createOneToOne<LO, GO, NT> (map);
        }

        RCP<const MapType> gatherMap;
        if (numProcs == 1) {
          gatherMap = oneToOneMap;
        } else {
	  if (debug && myRank == 0) {
	    cerr << "- Gathering the Map elements" << endl;
	  }
          // Gather each process' count of Map elements to Proc 0,
          // into the recvCounts array.  This will tell Proc 0 how
          // many GIDs to expect from each process when calling
          // MPI_Gatherv.  Counts and offsets are all int, because
          // that's what MPI uses.  Teuchos::as will at least prevent
          // bad casts to int in a debug build.
          //
          // Yeah, it's not memory scalable.  Guess what: We're going
          // to bring ALL the data to Proc 0, not just the receive
          // counts.  The first MPI_Gather is only the beginning...
          // Seriously, if you want to make this memory scalable, the
          // right thing to do (after the Export to the one-to-one
          // Map) is to go round robin through the processes, having
          // each send a chunk of data (with its GIDs, to get the
          // order of rows right) at a time, and overlapping writing
          // to the file (resp. reading from it) with receiving (resp.
          // sending) the next chunk.
          const int myEltCount = as<int> (oneToOneMap->getNodeNumElements ());
          Array<int> recvCounts (numProcs);
          const int rootProc = 0;
          gather<int> (&myEltCount, 1, recvCounts.getRawPtr (), 1, rootProc, comm);

          ArrayView<const GO> myGlobalElts = oneToOneMap->getNodeElementList ();
          const int numMyGlobalElts = as<int> (myGlobalElts.size ());
          // Only Proc 0 needs to receive and store all the GIDs (from
          // all processes).
          ArrayRCP<GO> allGlobalElts;
          if (myRank == 0) {
            allGlobalElts = arcp<GO> (oneToOneMap->getGlobalNumElements ());
            std::fill (allGlobalElts.begin (), allGlobalElts.end (), 0);
          }
          //const int numAllGlobalElts = as<int> (allGlobalElts.size ());

	  if (debug) {
	    if (myRank == 0) {
	      cerr << "- Partial sum" << endl;
	    }
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	  }
          // Displacements for gatherv() in this case (where we are
          // gathering into a contiguous array) are an exclusive
          // partial sum (first entry is zero, second starts the
          // partial sum) of recvCounts.
          Array<int> displs (numProcs, 0);
          std::partial_sum (recvCounts.begin (), recvCounts.end () - 1, displs.begin () + 1);

	  if (debug) {
	    if (myRank == 0) {
	      cerr << "- Partial sum results:" << endl;
	    }
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	    std::ostringstream os;
	    os << "  - Proc " << myRank << ":" << endl
	       << "    - recvCounts (input): " << Teuchos::toString (recvCounts) << endl
	       << "    - displs (output): " << Teuchos::toString (displs) << endl;
	    for (int p = 0; p < numProcs; ++p) {
	      if (p == myRank) {
		cerr << os.str ();
	      }
	      comm->barrier ();
	      comm->barrier ();
	      comm->barrier ();
	    }
	  }

	  if (debug) {
	    if (myRank == 0) {
	      cerr << "- Gatherv" << endl;
	    }
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	    {
	      std::ostringstream os;
	      os << "  - Proc " << myRank << ":" << endl
		 << "    - myGlobalElts: " << Teuchos::toString (myGlobalElts) << endl
		 << "    - numMyGlobalElts: " << numMyGlobalElts << endl
		 << "    - myGlobalElts: " << Teuchos::toString (myGlobalElts) << endl
		 << "    - allGlobalElts: " << Teuchos::toString (allGlobalElts) << endl
		 << "    - recvCounts: " << Teuchos::toString (recvCounts) << endl
		 << "    - displs: " << Teuchos::toString (displs) << endl;
	      for (int p = 0; p < numProcs; ++p) {
		if (p == myRank) {
		  cerr << os.str ();
		}
		comm->barrier ();
		comm->barrier ();
		comm->barrier ();
	      }
	    }
	  }
          gatherv<GO> (myGlobalElts.getRawPtr (), numMyGlobalElts,
                       allGlobalElts.getRawPtr (), recvCounts.getRawPtr (), displs.getRawPtr (),
                       rootProc, comm);

          // Create a Map with all the GIDs, in the same order as in
          // the one-to-one Map, but owned by Proc 0.
          ArrayView<const GO> allElts (NULL, 0);
          if (myRank == 0) {
            allElts = allGlobalElts ();
          }
	  if (debug) {
	    if (myRank == 0) {
	      cerr << "- Creating Map:" << endl;
	    }
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	    std::ostringstream os;
	    os << "  - Proc " << myRank << endl
	       << "    - allElts: " << Teuchos::toString (allElts) << endl;
	    cerr << os.str ();
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	  } // if debug
          gatherMap = rcp (new MapType (Teuchos::OrdinalTraits<global_size_t>::invalid (),
					allElts, oneToOneMap->getIndexBase (),
					comm, map->getNode ()));
        }
	if (debug && myRank == 0) {
	  cerr << "- Done" << endl;
	}
        return gatherMap;
      }

    } // namespace (anonymous)

    /// \class Reader
    /// \brief Matrix Market file reader for CrsMatrix and MultiVector.
    /// \author Mark Hoemmen
    ///
    /// The Matrix Market (see its <a
    /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
    /// details) defines a human-readable text file format for storing
    /// sparse and dense matrices.  This class defines methods for
    /// reading sparse and dense matrices from a Matrix Market file or
    /// input stream.  It represents sparse matrices as CrsMatrix and
    /// dense vectors and matrices as MultiVector.  Reader can also
    /// read a Map (in the format produced by Writer) from a file or
    /// input stream.
    ///
    /// All methods of this class only open files or read from input
    /// streams on the MPI process with rank 0, with respect to the
    /// MPI communicator over which the given CrsMatrix or MultiVector
    /// is to be distributed.
    ///
    /// We define the MultiVector type returned by readDense() and
    /// readDenseFile() using the scalar_type, local_ordinal_type,
    /// global_ordinal_type, and node_type typedefs in
    /// SparseMatrixType.  This ensures that the multivectors returned
    /// by those methods have a type compatible with the CrsMatrix
    /// sparse matrices returned by readSparse() and readSparseFile().
    /// We do this because the typical use case of Matrix Market files
    /// in Trilinos is to test sparse matrix methods, which usually
    /// involves reading a sparse matrix A and perhaps also a dense
    /// right-hand side b.  Also, this lets you use CrsMatrix objects
    /// with non-default LocalMatOps template parameters.  (If we
    /// templated on Scalar, LocalOrdinal, GlobalOrdinal, and Node, we
    /// would also have to template on LocalMatOps in order to deal
    /// with CrsMatrix types with nondefault LocalMatOps.  That would
    /// tie Reader to CrsMatrix anyway, since MultiVector is not
    /// templated on LocalMatOps.  As a result, we might as well just
    /// template on the CrsMatrix type, in order to use arbitrary
    /// LocalMatOps types without additional code.)
    ///
    /// \tparam SparseMatrixType A specialization of CrsMatrix.
    ///
    /// Templating on the specialization of CrsMatrix means that the
    /// Reader expects matrix data of a type compatible with the
    /// CrsMatrix's scalar_type.  In general, Matrix Market files may
    /// contain data of integer, real, or complex type.  However, the
    /// reader methods have to return a CrsMatrix of a specific type,
    /// so we require that you declare a Reader with the CrsMatrix
    /// type that you want and that you expect the file(s) to contain.
    ///
    /// We didn't find any of the alternatives to this approach
    /// acceptable.  One possibility would have been to have the
    /// reader methods return a "container that can hold anything,"
    /// like a boost::any.  However, then you would have to know all
    /// five template arguments of the CrsMatrix in order to get the
    /// actual CrsMatrix object out.  C++ doesn't have algebraic data
    /// types (see the Wikipedia entry for a good definition) that are
    /// disjoint unions of different types.  Thus, we couldn't have
    /// had the readers return a CrsMatrix with scalar_type = "int or
    /// double or complex<double>."  While you can implement such a
    /// type in C++ (see e.g., boost::variant), it would not be
    /// interchangeable for its component types.  This is because it
    /// may not have the same memory layout (e.g., copying an array of
    /// boost::variant<int, double, complex<double> > bitwise into an
    /// array of int may not work).
    template<class SparseMatrixType>
    class Reader {
    public:
      //! This class' template parameter; a specialization of CrsMatrix.
      typedef SparseMatrixType sparse_matrix_type;
      typedef RCP<sparse_matrix_type> sparse_matrix_ptr;

      /// Type of the entries of the sparse matrix.
      /// The first template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// Type of the local indices of the sparse matrix.
      /// The second template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// Type of the global indices of the sparse matrix.
      ///
      /// The third template parameter of CrsMatrix and MultiVector.
      /// This is also the type of indices as read from the Matrix
      /// Market file.  Indices of the sparse matrix are read in as
      /// global ordinals, since Matrix Market files represent the
      /// whole matrix and don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type
        global_ordinal_type;
      /// The Kokkos Node type.
      /// The fourth template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::node_type node_type;

      //! The MultiVector specialization associated with SparseMatrixType.
      typedef MultiVector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> multivector_type;

      typedef RCP<node_type> node_ptr;
      typedef Comm<int> comm_type;
      typedef RCP<const comm_type> comm_ptr;
      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      typedef RCP<const map_type> map_ptr;

    private:
      /// \typedef size_type
      /// \brief Handy typedef for entries of arrays such as rowPtr.
      typedef typename ArrayRCP<global_ordinal_type>::size_type size_type;

      /// \brief Compute initial range map.
      ///
      /// Range maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pComm [in] Global communicator.
      /// \param pNode [in] Kokkos Node object.
      /// \param numRows [in] Global number of rows in the matrix.
      ///
      /// \return Range map to be used for constructing a CrsMatrix.
      static RCP<const map_type>
      makeRangeMap (const RCP<const comm_type>& pComm,
                    const RCP<node_type>& pNode,
                    const global_ordinal_type numRows)
      {
        // A conventional, uniformly partitioned, contiguous map.
        return rcp (new map_type (static_cast<global_size_t> (numRows),
                                  static_cast<global_ordinal_type> (0),
                                  pComm, GloballyDistributed, pNode));
      }

      /// \brief Compute initial row map, or verify an existing one.
      ///
      /// The typical case when reading a sparse matrix from a file is
      /// for the reader itself to create a new row map, in particular
      /// a standard uniform contiguous one-to-one row map.  However,
      /// we also give the option to use an existing row map, if you
      /// are already using a particular distribution for (say) vector
      /// data and don't want to stop using it.  In the latter case
      /// (pRowMap is not null), we validate the communicator and node
      /// of the existing row map that you pass in.  In either case,
      /// you need to know the (global) number of rows in the matrix.
      ///
      /// \param pRowMap [in] If non-null, test pRowMap for validity,
      ///   and return it if valid.  Otherwise, if pRowMap is null,
      ///   initialize and return a (uniform contiguous one-to-one)
      ///   row map.  "Validity" here means that the map's
      ///   communicator and node are the same objects (pointerwise)
      ///   as the corresponding arguments.  (Note that the global
      ///   number of elements may not be the same as the number of
      ///   rows; a row map is not required to be one-to-one.)  The
      ///   typical case is to pass in null here, which is why we call
      ///   this routine "makeRowMap".
      /// \param pComm [in] Global communicator.
      /// \param pNode [in] Kokkos Node object.
      /// \param numRows [in] Global number of rows in the matrix.  If
      ///   pRowMap is nonnull, used only for error checking.
      ///
      /// \return If pRowMap is null, a new row map, otherwise pRowMap.
      static RCP<const map_type>
      makeRowMap (const RCP<const map_type>& pRowMap,
                  const RCP<const comm_type>& pComm,
                  const RCP<node_type>& pNode,
                  const global_ordinal_type numRows)
      {
        // If the caller didn't provide a map, return a conventional,
        // uniformly partitioned, contiguous map.
        if (pRowMap.is_null()) {
          return rcp (new map_type (static_cast<global_size_t> (numRows),
                                    static_cast<global_ordinal_type> (0),
                                    pComm, GloballyDistributed, pNode));
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(! pRowMap->isDistributed() && pComm->getSize() > 1,
                             std::invalid_argument,
                             "The specified row map is not distributed, but "
                             "the given communicator includes more than one "
                             "rank (in fact, there are " << pComm->getSize()
                             << " ranks).");
          TEUCHOS_TEST_FOR_EXCEPTION(pRowMap->getComm() != pComm,
                             std::invalid_argument,
                             "The specified row map's communicator (pRowMap->"
                             "getComm()) is different than the given separately "
                             "supplied communicator pComm.");
          TEUCHOS_TEST_FOR_EXCEPTION(pRowMap->getNode() != pNode,
                             std::invalid_argument,
                             "The specified row map's node (pRowMap->getNode()) "
                             "is different than the given separately supplied "
                             "node pNode.");
          return pRowMap;
        }
      }

      /// \brief Compute domain map.
      ///
      /// Domain maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pRangeMap [in] Valid range map of the matrix,
      ///   as returned by \c makeRangeMap().
      /// \param numRows [in] Global number of rows in the matrix.
      /// \param numCols [in] Global number of columns in the matrix.
      ///
      /// \return The domain map.  If numRows == numCols, this is
      ///   identical to the range map, otherwise we make a new map
      ///   for the domain.
      static map_ptr
      makeDomainMap (const map_ptr& pRangeMap,
                     const global_ordinal_type numRows,
                     const global_ordinal_type numCols)
      {
        // Abbreviations so that the map creation call isn't too long.
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type Node;

        if (numRows == numCols) {
          return pRangeMap;
        } else {
          comm_ptr pComm = pRangeMap->getComm();
          node_ptr pNode = pRangeMap->getNode();
          return createUniformContigMapWithNode<LO,GO,Node> (numCols,
                                                             pComm,
                                                             pNode);
        }
      }

      /// \brief Distribute the sparse matrix entries.
      ///
      /// This is one of those routines that just has to be messy.  We try
      /// to comment everything so you can see what's going on.  All of the
      /// matrix data starts out on Rank 0, stored in CSR format (rowPtr,
      /// colInd, values).  numEntriesPerRow applies to this data; it can be
      /// computed from rowPtr, and has one less entry than rowPtr.  Rank 0
      /// then engages in a dialog with Ranks 1 .. numProcs-1 to determine
      /// which part of the matrix data is theirs, and to send them their
      /// data.
      ///
      /// \param myNumEntriesPerRow [out] For my row indices, the number of
      ///   entries per row.  This array has
      ///   pRowMap->getNodeElementList().size() entries, and is indexed in
      ///   the same way, so that myNumEntriesPerRow[k] is the number of
      ///   entries in row pRowMap->getNodeElementList()[k].
      ///
      /// \param myRowPtr [out] The row pointer array for the rows
      ///   that belong to me.  This array has one more entry than
      ///   myNumEntriesPerRow, and is the prefix sum (with
      ///   myRowPtr[0] = 0) of myNumEntriesPerRow.
      ///
      /// \param myColInd [out] My rows' column indices.  If myRows =
      ///   pRowMap->getNodeElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then
      ///   myColInd[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param myValues [out] My rows' stored matrix values.  If myRows =
      ///   pRowMap->getNodeElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then
      ///   myValues[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param pRowMap [in] Map describing the distribution of rows among
      ///   processors.
      ///
      /// \param numEntriesPerRow [in/out] For all row indices, the number
      ///   of entries per row.  You can construct this from the usual CSR
      ///   matrix "rowPtr" array by differences: numEntriesPerRow[k] =
      ///   rowPtr[k+1] - rowPtr[k].  As a result, numEntriesPerRow has one
      ///   less entry than rowPtr.  On output, the reference is invalidated
      ///   to save space.
      ///
      /// \param rowPtr [in/out] On input: the usual CSR matrix row pointer
      ///   array for the whole matrix.  On output, the reference is
      ///   invalidated to save space.
      ///
      /// \param colInd [in/out] On input: all rows' column indices.  If
      ///   start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then colInd[start .. end-1] are the column
      ///   indices for row k.  On output, the reference is invalidated to
      ///   save space.
      ///
      /// \param values [in/out] On input: all rows' stored matrix values.
      ///   If start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then values[start .. end-1] are the values
      ///   for row k.  On output, the reference is invalidated to save
      ///   space.
      ///
      /// \param debug [in] If true, print copious debugging output to
      ///   stderr on Rank 0.  This option is unlikely to be useful to
      ///   anyone but a Tpetra developer debugging this code.
      ///
      /// \note It was only after I wrote this routine that I learned
      ///   it was completely unnecessary.  All the data
      ///   redistribution can be done in about 10 lines of code by
      ///   using Tpetra::Map objects, and either Import or Export.
      ///   (For example, you could read the file into the sparse
      ///   matrix entirely owned by Proc 0, then export it using a
      ///   distributed Map.)  However, this routine works and I
      ///   haven't had time to rewrite it yet.  Just expect that the
      ///   implementation of \c readSparse() may become a lot shorter
      ///   in the future.
      static void
      distribute (ArrayRCP<size_t>& myNumEntriesPerRow,
                  ArrayRCP<size_t>& myRowPtr,
                  ArrayRCP<global_ordinal_type>& myColInd,
                  ArrayRCP<scalar_type>& myValues,
                  const RCP<const map_type>& pRowMap,
                  ArrayRCP<size_t>& numEntriesPerRow,
                  ArrayRCP<size_t>& rowPtr,
                  ArrayRCP<global_ordinal_type>& colInd,
                  ArrayRCP<scalar_type>& values,
                  const bool debug=false)
      {
         using Teuchos::as;
         using Teuchos::CommRequest;
         using Teuchos::receive;
         using Teuchos::send;
         using std::cerr;
         using std::endl;

         const bool extraDebug = false;
         comm_ptr pComm = pRowMap->getComm();
         const int numProcs = Teuchos::size (*pComm);
         const int myRank = Teuchos::rank (*pComm);
         const int rootRank = 0;

         // Type abbreviations to make the code more concise.
         typedef global_ordinal_type GO;
         typedef local_ordinal_type LO;

         // List of the global indices of my rows.  They may or may
         // not be contiguous, and the row map need not be one-to-one.
         ArrayView<const GO> myRows = pRowMap->getNodeElementList();
         const size_type myNumRows = myRows.size();
         TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(myNumRows) !=
                            pRowMap->getNodeNumElements(),
                            std::logic_error,
                            "pRowMap->getNodeElementList().size() = "
                            << myNumRows
                            << " != pRowMap->getNodeNumElements() = "
                            << pRowMap->getNodeNumElements() << ".  "
                            "Please report this bug to the Tpetra developers.");
         TEUCHOS_TEST_FOR_EXCEPTION(myRank == 0 && numEntriesPerRow.size() < myNumRows,
                            std::logic_error,
                            "On Proc 0: numEntriesPerRow.size() = "
                            << numEntriesPerRow.size()
                            << " != pRowMap->getNodeElementList().size() = "
                            << myNumRows << ".  Please report this bug to the "
                            "Tpetra developers.");

         // Space for my proc's number of entries per row.
         // Will be filled in below.
         myNumEntriesPerRow = arcp<size_t> (myNumRows);

         if (myRank != rootRank)
           {
             // Tell the root how many rows we have.  If we're sending
             // none, then we don't have anything else to send, nor
             // does the root have to receive anything else.
             send (*pComm, myNumRows, rootRank);
             if (myNumRows != 0)
               {
                 // Now send my rows' global indices.  Hopefully the
                 // cast to int doesn't overflow.  This is unlikely,
                 // since it should fit in a LO, even
                 // though it is a GO.
                 send (*pComm, static_cast<int> (myNumRows),
                       myRows.getRawPtr(), rootRank);

                 // I (this proc) don't care if my global row indices
                 // are contiguous, though the root proc does (since
                 // otherwise it needs to pack noncontiguous data into
                 // contiguous storage before sending).  That's why we
                 // don't check for contiguousness here.

                 // Ask the root processor for my part of the array of the
                 // number of entries per row.
                 receive (*pComm, rootRank,
                          static_cast<int> (myNumRows),
                          myNumEntriesPerRow.getRawPtr());

                 // Use the resulting array to figure out how many column
                 // indices and values for which I should ask from the root
                 // processor.
                 const local_ordinal_type myNumEntries =
                   std::accumulate (myNumEntriesPerRow.begin(),
                                    myNumEntriesPerRow.end(),
                                    0);

                 // Make space for my entries of the sparse matrix.
                 // Note that they don't have to be sorted by row
                 // index.  Iterating through all my rows requires
                 // computing a running sum over myNumEntriesPerRow.
                 myColInd = arcp<GO> (myNumEntries);
                 myValues = arcp<scalar_type> (myNumEntries);
                 if (myNumEntries > 0)
                   { // Ask for that many column indices and values,
                     // if there are any.
                     receive (*pComm, rootRank,
                              static_cast<int> (myNumEntries),
                              myColInd.getRawPtr());
                     receive (*pComm, rootRank,
                              static_cast<int> (myNumEntries),
                              myValues.getRawPtr());
                   }
               } // If I own at least one row
           } // If I am not the root processor
         else { // I _am_ the root processor
	   if (debug) {
	     cerr << "-- Proc 0: Copying my data from global arrays" << endl;
	   }
	   // Proc 0 still needs to (allocate, if not done already)
	   // and fill its part of the matrix (my*).
	   for (size_type k = 0; k < myNumRows; ++k) {
	     const GO myCurRow = myRows[k];
	     const local_ordinal_type numEntriesInThisRow = numEntriesPerRow[myCurRow];
	     //myNumEntriesPerRow[k] = numEntriesPerRow[myCurRow];
	     myNumEntriesPerRow[k] = numEntriesInThisRow;
	   }
	   if (extraDebug && debug) {
	     cerr << "Proc " << Teuchos::rank (*(pRowMap->getComm()))
		  << ": myNumEntriesPerRow[0.." << (myNumRows-1) << "] = [";
	     for (size_type k = 0; k < myNumRows; ++k) {
	       cerr << myNumEntriesPerRow[k];
	       if (k < myNumRows-1) {
		 cerr << " ";
	       }
	     }
	     cerr << "]" << endl;
	   }
	   // The total number of matrix entries that my proc owns.
	   const local_ordinal_type myNumEntries =
	     std::accumulate (myNumEntriesPerRow.begin(),
			      myNumEntriesPerRow.end(), 0);
	   if (debug) {
	     cerr << "-- Proc 0: I own " << myNumRows << " rows and "
		  << myNumEntries << " entries" << endl;
	   }
	   myColInd = arcp<GO> (myNumEntries);
	   myValues = arcp<scalar_type> (myNumEntries);

	   // Copy Proc 0's part of the matrix into the my* arrays.
	   // It's important that myCurPos be updated _before_ k,
	   // otherwise myCurPos will get the wrong number of entries
	   // per row (it should be for the row in the just-completed
	   // iteration, not for the next iteration's row).
	   local_ordinal_type myCurPos = 0;
	   for (size_type k = 0; k < myNumRows;
		myCurPos += myNumEntriesPerRow[k], ++k) {
	     const local_ordinal_type curNumEntries = myNumEntriesPerRow[k];
	     const GO myRow = myRows[k];
	     const size_t curPos = rowPtr[myRow];
	     if (extraDebug && debug) {
	       cerr << "k = " << k << ", myRow = " << myRow << ": colInd("
		    << curPos << "," << curNumEntries << "), myColInd("
		    << myCurPos << "," << curNumEntries << ")" << endl;
	     }
	     // Only copy if there are entries to copy, in order
	     // not to construct empty ranges for the ArrayRCP
	     // views.
	     if (curNumEntries > 0) {
	       ArrayView<GO> colIndView = colInd(curPos, curNumEntries);
	       ArrayView<GO> myColIndView =
		 myColInd(myCurPos, curNumEntries);
	       std::copy (colIndView.begin(), colIndView.end(),
			  myColIndView.begin());

	       ArrayView<scalar_type> valuesView =
		 values(curPos, curNumEntries);
	       ArrayView<scalar_type> myValuesView =
		 myValues(myCurPos, curNumEntries);
	       std::copy (valuesView.begin(), valuesView.end(),
			  myValuesView.begin());
	     }
	   }

	   // Proc 0 processes each other proc p in turn.
	   for (int p = 1; p < numProcs; ++p) {
	     if (debug) {
	       cerr << "-- Proc 0: Processing proc " << p << endl;
	     }

	     size_type theirNumRows = 0;
	     // Ask Proc p how many rows it has.  If it doesn't
	     // have any, we can move on to the next proc.  This
	     // has to be a standard receive so that we can avoid
	     // the degenerate case of sending zero data.
	     receive (*pComm, p, &theirNumRows);
	     if (debug) {
	       cerr << "-- Proc 0: Proc " << p << " owns "
		    << theirNumRows << " rows" << endl;
	     }
	     if (theirNumRows != 0) {
	       // Ask Proc p which rows it owns.  The resulting global
	       // row indices are not guaranteed to be contiguous or
	       // sorted.  Global row indices are themselves indices
	       // into the numEntriesPerRow array.
	       ArrayRCP<GO> theirRows = arcp<GO> (theirNumRows);
	       receive (*pComm, p, as<int> (theirNumRows),
			theirRows.getRawPtr ());
	       // Extra test to make sure that the rows we received
	       // are all sensible.  This is a good idea since we are
	       // going to use the global row indices we've received
	       // to index into the numEntriesPerRow array.  Better to
	       // catch any bugs here and print a sensible error
	       // message, rather than segfault and print a cryptic
	       // error message.
	       {
		 const global_size_t numRows = pRowMap->getGlobalNumElements();
		 bool theirRowsValid = true;
		 for (size_type k = 0; k < theirNumRows; ++k) {
		   // global_ordinal_type is generally signed.
		   if (theirRows[k] < 0) {
		     theirRowsValid = false;
		   }
		   // Same-size signed->unsigned cast never overflows.
		   else if (as<global_size_t> (theirRows[k]) >= numRows) {
		     theirRowsValid = false;
		   }
		 }
		 if (! theirRowsValid) {
		   TEUCHOS_TEST_FOR_EXCEPTION(
	             ! theirRowsValid, std::logic_error,
		     "Proc " << p << " has at least one invalid row index.  "
		     "Here are all of them: " << Teuchos::toString (theirRows));
		 }
	       }

	       // Perhaps we could save a little work if we check
	       // whether Proc p's row indices are contiguous.  That
	       // would make lookups in the global data arrays
	       // faster.  For now, we just implement the general
	       // case and don't prematurely optimize.  (Remember
	       // that you're making Proc 0 read the whole file, so
	       // you've already lost scalability.)

	       // Compute the number of entries in each of Proc p's
	       // rows.  (Proc p will compute its row pointer array
	       // on its own, after it gets the data from Proc 0.)
	       ArrayRCP<size_t> theirNumEntriesPerRow;
	       theirNumEntriesPerRow = arcp<size_t> (theirNumRows);
	       for (size_type k = 0; k < theirNumRows; ++k) {
		 theirNumEntriesPerRow[k] = numEntriesPerRow[theirRows[k]];
	       }

	       // Tell Proc p the number of entries in each of its
	       // rows.  Hopefully the cast to int doesn't overflow.
	       // This is unlikely, since it should fit in a LO,
	       // even though it is a GO.
	       send (*pComm, static_cast<int> (theirNumRows),
		     theirNumEntriesPerRow.getRawPtr(), p);

	       // Figure out how many entries Proc p owns.
	       const local_ordinal_type theirNumEntries =
		 std::accumulate (theirNumEntriesPerRow.begin(),
				  theirNumEntriesPerRow.end(), 0);

	       if (debug) {
		 cerr << "-- Proc 0: Proc " << p << " owns "
		      << theirNumEntries << " entries" << endl;
	       }

	       // If there are no entries to send, then we're done
	       // with Proc p.
	       if (theirNumEntries == 0) {
		 continue;
	       }

	       // Construct (views of) proc p's column indices
	       // and values.  Later, we might like to optimize
	       // for the (common) contiguous case, for which we
	       // don't need to copy data into separate "their*"
	       // arrays (we can just use contiguous views of
	       // the global arrays).
	       ArrayRCP<GO> theirColInd (theirNumEntries);
	       ArrayRCP<scalar_type> theirValues (theirNumEntries);
	       // Copy Proc p's part of the matrix into the their*
	       // arrays.  It's important that theirCurPos be
	       // updated _before_ k, otherwise theirCurPos will get
	       // the wrong number of entries per row (it should be
	       // for the row in the just-completed iteration, not
	       // for the next iteration's row).
	       local_ordinal_type theirCurPos = 0;
	       for (size_type k = 0; k < theirNumRows;
		    theirCurPos += theirNumEntriesPerRow[k], k++) {
		 const local_ordinal_type curNumEntries = theirNumEntriesPerRow[k];
		 const GO theirRow = theirRows[k];
		 const local_ordinal_type curPos = rowPtr[theirRow];

		 // Only copy if there are entries to copy, in order
		 // not to construct empty ranges for the ArrayRCP
		 // views.
		 if (curNumEntries > 0) {
		   ArrayView<GO> colIndView =
		     colInd (curPos, curNumEntries);
		   ArrayView<GO> theirColIndView =
		     theirColInd (theirCurPos, curNumEntries);
		   std::copy (colIndView.begin(), colIndView.end(),
			      theirColIndView.begin());

		   ArrayView<scalar_type> valuesView =
		     values (curPos, curNumEntries);
		   ArrayView<scalar_type> theirValuesView =
		     theirValues (theirCurPos, curNumEntries);
		   std::copy (valuesView.begin(), valuesView.end(),
			      theirValuesView.begin());
		 }
	       }
	       // Send Proc p its column indices and values.
	       // Hopefully the cast to int doesn't overflow.  This
	       // is unlikely, since it should fit in a LO, even
	       // though it is a GO.
	       send (*pComm, static_cast<int> (theirNumEntries),
		     theirColInd.getRawPtr(), p);
	       send (*pComm, static_cast<int> (theirNumEntries),
		     theirValues.getRawPtr(), p);

	       if (debug) {
		 cerr << "-- Proc 0: Finished with proc " << p << endl;
	       }
	     } // If proc p owns at least one row
	   } // For each proc p not the root proc 0
	 } // If I'm (not) the root proc 0

         // Invalidate the input data to save space, since we don't
         // need it anymore.
         numEntriesPerRow = null;
         rowPtr = null;
         colInd = null;
         values = null;

         if (debug && myRank == 0) {
           cerr << "-- Proc 0: About to fill in myRowPtr" << endl;
	 }

         // Allocate and fill in myRowPtr (the row pointer array for
         // my rank's rows).  We delay this until the end because we
         // don't need it to compute anything else in distribute().
         // Each proc can do this work for itself, since it only needs
         // myNumEntriesPerRow to do so.
         myRowPtr = arcp<size_t> (myNumRows+1);
         myRowPtr[0] = 0;
         for (size_type k = 1; k < myNumRows+1; ++k) {
           myRowPtr[k] = myRowPtr[k-1] + myNumEntriesPerRow[k-1];
         }
         if (extraDebug && debug) {
	   cerr << "Proc " << Teuchos::rank (*(pRowMap->getComm()))
		<< ": myRowPtr[0.." << myNumRows << "] = [";
	   for (size_type k = 0; k < myNumRows+1; ++k) {
	     cerr << myRowPtr[k];
	     if (k < myNumRows) {
	       cerr << " ";
	     }
	   }
	   cerr << "]" << endl << endl;
	 }

         if (debug && myRank == 0) {
           cerr << "-- Proc 0: Done with distribute" << endl;
	 }
      }

      /// \brief Given my proc's data, return the completed sparse matrix.
      ///
      /// Each proc inserts its data into the sparse matrix, and then,
      /// if callFillComplete is true, all procs call fillComplete().
      /// (For whatever reason, you might not be done with the matrix
      /// yet, so you might want to call fillComplete() yourself.
      /// CrsMatrix::fillResume() doesn't currently work as you might
      /// expect when storage optimization is enabled; it fixes the
      /// graph of the matrix, so that you can't add new entries.)
      ///
      static sparse_matrix_ptr
      makeMatrix (ArrayRCP<size_t>& myNumEntriesPerRow,
                  ArrayRCP<size_t>& myRowPtr,
                  ArrayRCP<global_ordinal_type>& myColInd,
                  ArrayRCP<scalar_type>& myValues,
                  const map_ptr& pRowMap,
                  const map_ptr& pRangeMap,
                  const map_ptr& pDomainMap,
                  const bool callFillComplete = true)
      {
        using std::cerr;
        using std::endl;
        // Typedef to make certain type declarations shorter.
        typedef global_ordinal_type GO;

        const bool extraDebug = false;
        const bool debug = false;

        // The row pointer array always has at least one entry, even
        // if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
        // and myValues would all be empty arrays in that degenerate
        // case, but the row and domain maps would still be nonnull
        // (though they would be trivial maps).
        TEUCHOS_TEST_FOR_EXCEPTION(myRowPtr.is_null(), std::logic_error,
                           "makeMatrix: myRowPtr array is null.  "
                           "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pDomainMap.is_null(), std::logic_error,
                           "makeMatrix: domain map is null.  "
                           "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRangeMap.is_null(), std::logic_error,
                           "makeMatrix: range map is null.  "
                           "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRowMap.is_null(), std::logic_error,
                           "makeMatrix: row map is null.  "
                           "Please report this bug to the Tpetra developers.");

        // Handy for debugging output; not needed otherwise.
        const int myRank = Teuchos::rank (*(pRangeMap->getComm()));

        if (extraDebug && debug)
          {
            cerr << "Proc " << myRank << ":" << endl
                 << "-- myRowPtr = [ ";
            std::copy (myRowPtr.begin(), myRowPtr.end(),
                       std::ostream_iterator<size_type>(cerr, " "));
            cerr << "]" << endl << "-- myColInd = [ ";
            std::copy (myColInd.begin(), myColInd.end(),
                       std::ostream_iterator<size_type>(cerr, " "));
            cerr << "]" << endl << endl;
          }

        // Go through all of my columns, and see if any are not in the
        // domain map.  This is possible if numProcs > 1, otherwise
        // not.
        if (extraDebug && debug)
          {
            size_type numRemote = 0;
            std::vector<GO> remoteGIDs;

            typedef typename ArrayRCP<GO>::const_iterator iter_type;
            for (iter_type it = myColInd.begin(); it != myColInd.end(); ++it)
              {
                if (! pDomainMap->isNodeGlobalElement (*it))
                  {
                    numRemote++;
                    remoteGIDs.push_back (*it);
                  }
              }
            if (numRemote > 0)
              {
                cerr << "Proc " << myRank << ": " << numRemote
                     << " remote GIDs = [ " << endl;
                std::copy (remoteGIDs.begin(), remoteGIDs.end(),
                           std::ostream_iterator<GO>(cerr, " "));
                cerr << "]" << endl;
              }
          }

        // Construct the CrsMatrix, using the row map, with the
        // constructor specifying the number of nonzeros for each row.
        // Create with DynamicProfile, so that the fillComplete() can
        // do first-touch reallocation (a NUMA (Non-Uniform Memory
        // Access) optimization on multicore CPUs).
        sparse_matrix_ptr A =
          rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow,
                                       DynamicProfile));
        TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::logic_error,
                           "makeMatrix: Initial allocation of CrsMatrix failed"
                           ".  Please report this bug to the Tpetra developers"
                           ".");
        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = pRowMap->getNodeElementList();
        const size_type myNumRows = myRows.size();

        // Add this processor's matrix entries to the CrsMatrix.
        for (size_type k = 0; k < myNumRows; ++k)
          {
            const size_type myCurPos = myRowPtr[k];
            const local_ordinal_type curNumEntries = myNumEntriesPerRow[k];

            if (extraDebug && debug)
              {
                cerr << "Proc " << myRank << ": k = " << k
                     << ", myCurPos = " << myCurPos
                     << ", curNumEntries = " << curNumEntries
                     << endl;
              }
            // Avoid constructing empty views of ArrayRCP objects.
            if (curNumEntries > 0)
              A->insertGlobalValues (myRows[k],
                                     myColInd(myCurPos, curNumEntries),
                                     myValues(myCurPos, curNumEntries));
          }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        if (callFillComplete)
          A->fillComplete (pDomainMap, pRangeMap);
        return A;
      }

    private:

      /// \brief Read in the Banner line from the given input stream.
      ///
      /// Only call this method on one (MPI communicator) process.
      ///
      /// \param in [in/out] Input stream from which to read the
      ///   Banner line.
      /// \param lineNumber [in/out] On input: Current line number of
      ///   the input stream.  On output: if any line(s) were
      ///   successfully read from the input stream, this is
      ///   incremented by the number of line(s) read.  (This includes
      ///   comment lines.)
      /// \param tolerant [in] Whether to parse tolerantly.
      /// \param debug [in] Whether to write debugging output to
      ///   stderr.
      ///
      /// \return Banner [non-null]
      static RCP<const Teuchos::MatrixMarket::Banner>
      readBanner (std::istream& in,
                  size_t& lineNumber,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        RCP<Banner> pBanner; // On output, if successful: the read-in Banner.
        std::string line; // If read from stream successful: the Banner line

        // Try to read a line from the input stream.
        const bool readFailed = ! getline(in, line);
        TEUCHOS_TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
          "Failed to get Matrix Market banner line from input.");

        // We read a line from the input stream.
        lineNumber++;

        // Assume that the line we found is the Banner line.
        try {
          pBanner = rcp (new Banner (line, tolerant));
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            "Matrix Market banner line contains syntax error(s): "
            << e.what());
        }
        TEUCHOS_TEST_FOR_EXCEPTION(pBanner->objectType() != "matrix",
          std::invalid_argument, "The Matrix Market file does not contain "
          "matrix data.  Its Banner (first) line says that its object type is \""
          << pBanner->matrixType() << "\", rather than the required \"matrix\".");

        // Validate the data type of the matrix, with respect to the
        // Scalar type of the CrsMatrix entries.
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! STS::isComplex && pBanner->dataType() == "complex",
          std::invalid_argument,
          "The Matrix Market file contains complex-valued data, but you are "
          "trying to read it into a matrix containing entries of the real-"
          "valued Scalar type \""
          << Teuchos::TypeNameTraits<scalar_type>::name() << "\".");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pBanner->dataType() != "real" &&
          pBanner->dataType() != "complex" &&
          pBanner->dataType() != "integer",
          std::invalid_argument,
          "When reading Matrix Market data into a Tpetra::CrsMatrix, the "
          "Matrix Market file may not contain a \"pattern\" matrix.  A "
          "pattern matrix is really just a graph with no weights.  It "
          "should be stored in a CrsGraph, not a CrsMatrix.");

        return pBanner;
      }

      /// \brief Read sparse matrix dimensions on Rank 0, and broadcast.
      ///
      /// Call on all MPI ranks.  MPI Rank 0 attempts to read in the
      /// coordinate dimensions from the input stream.  If it
      /// succeeds, it broadcasts them to all the other MPI ranks.
      /// (All ranks need to know the matrix dimensions in order to
      /// create domain, range, and column Maps.)
      ///
      /// \param in [in/out, valid only on Rank 0] Input stream from
      ///   which to read the sparse matrix dimensions.
      /// \param lineNumber [in/out, valid only on Rank 0] On input:
      ///   Current line number of the input stream.  On output: if
      ///   any line(s) were successfully read from the input stream,
      ///   this is incremented by the number of line(s) read.  (This
      ///   includes comment lines.)
      /// \param pComm [in, global] Communicator over which the matrix
      ///   will (eventually -- not here) be distributed.
      /// \param tolerant [in] Whether to parse tolerantly.
      /// \param debug [in] Whether to write debugging output to
      ///   stderr on MPI Proc 0.
      ///
      /// \return (numRows, numCols, numNonzeros)
      static Tuple<global_ordinal_type, 3>
      readCoordDims (std::istream& in,
                     size_t& lineNumber,
                     const RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                     const comm_ptr& pComm,
                     const bool tolerant = false,
                     const bool debug = false)
      {
        using Teuchos::MatrixMarket::readCoordinateDimensions;

        // Packed coordinate matrix dimensions (numRows, numCols,
        // numNonzeros); computed on Rank 0 and broadcasted to all MPI
        // ranks.
        Tuple<global_ordinal_type, 3> dims;

        // Read in the coordinate matrix dimensions from the input
        // stream.  "success" tells us whether reading in the
        // coordinate matrix dimensions succeeded ("Guilty unless
        // proven innocent").
        bool success = false;
        if (pComm->getRank() == 0) {
          TEUCHOS_TEST_FOR_EXCEPTION(pBanner->matrixType() != "coordinate",
            std::invalid_argument, "The Tpetra::CrsMatrix Matrix Market reader "
            "only accepts \"coordinate\" (sparse) matrix data.");
          // Unpacked coordinate matrix dimensions
          global_ordinal_type numRows, numCols, numNonzeros;
          // Only MPI Rank 0 reads from the input stream
          success = readCoordinateDimensions (in, numRows, numCols,
                                              numNonzeros, lineNumber,
                                              tolerant);
          // Pack up the data into a Tuple so we can send them with
          // one broadcast instead of three.
          dims[0] = numRows;
          dims[1] = numCols;
          dims[2] = numNonzeros;
        }
        // Only Rank 0 did the reading, so it decides success.
        //
        // FIXME (mfh 02 Feb 2011) Teuchos::broadcast doesn't know how
        // to send bools.  For now, we convert to/from int instead,
        // using the usual "true is 1, false is 0" encoding.
        {
          int the_success = success ? 1 : 0; // only matters on MPI Rank 0
          Teuchos::broadcast (*pComm, 0, &the_success);
          success = (the_success == 1);
        }
        if (success) {
          // Broadcast (numRows, numCols, numNonzeros) from Rank 0
          // to all the other MPI ranks.
          Teuchos::broadcast (*pComm, 0, dims);
        }
        else {
          // Perhaps in tolerant mode, we could set all the
          // dimensions to zero for now, and deduce correct
          // dimensions by reading all of the file's entries and
          // computing the max(row index) and max(column index).
          // However, for now we just error out in that case.
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            "Error reading Matrix Market sparse matrix: failed to read "
            "coordinate matrix dimensions.");
        }
        return dims;
      }

      /// \typedef adder_type
      /// \brief Type of object that "adds" entries to the sparse matrix.
      ///
      /// "Adds" here means that it collects and makes note of matrix
      /// entries read in from the input stream.  This object doesn't
      /// call insertGlobalEntries() or fillComplete() on the
      /// CrsMatrix.  Depending on the Matrix Market banner line
      /// information, it may "symmetrize" the matrix by adding entry
      /// A(j,i) (with the appropriate value depending on the symmetry
      /// type) if entry A(i,j) is seen.
      typedef Teuchos::MatrixMarket::SymmetrizingAdder<Teuchos::MatrixMarket::Raw::Adder<scalar_type, global_ordinal_type> > adder_type;

      /// \brief Make an "adder" object for processing matrix data.
      ///
      /// \param pComm [in] Communicator (across whose MPI ranks
      ///   the sparse matrix will be distributed)
      ///
      /// \param banner [in, nonnull and valid on Rank 0 only]
      ///   Object describing the type and symmetry of matrix data.
      ///
      /// \param dims [in] (numRows, numCols, numEntries).  These are
      ///   the "expected" values as read from the top of the Matrix
      ///   Market input stream.  Whether they are the final values
      ///   depends on the "tolerant" parameter and the actual matrix
      ///   data read from the input stream.
      ///
      /// \param tolerant [in] Whether the adder should be "tolerant"
      ///   of syntax errors and missing/incorrect metadata.  (In
      ///   particular, this refers to the number of rows, columns,
      ///   and entries in the matrix.)
      /// \param debug [in] Whether to print verbose debug output
      ///   to stderr on Rank 0.
      ///
      /// \return An adder_type object [nonnull and valid on Rank 0
      ///   only] that optionally symmetrizes the entries of the
      ///   sparse matrix.
      ///
      static RCP<adder_type>
      makeAdder (const RCP<const Comm<int> >& pComm,
                 RCP<const Teuchos::MatrixMarket::Banner>& pBanner,
                 const Tuple<global_ordinal_type, 3>& dims,
                 const bool tolerant=false,
                 const bool debug=false)
      {
        if (pComm->getRank() == 0) {
          typedef Teuchos::MatrixMarket::Raw::Adder<scalar_type, global_ordinal_type> raw_adder_type;
          RCP<raw_adder_type> pRaw =
            rcp (new raw_adder_type (dims[0], dims[1], dims[2], tolerant, debug));
          return rcp (new adder_type (pRaw, pBanner->symmType()));
        }
        else {
          return null;
        }
      }

    public:

      /// \brief Read sparse matrix from the given Matrix Market file.
      ///
      /// Open the given file on MPI Rank 0 (with respect to the given
      /// communicator).  The file should contain Matrix Market
      /// "coordinate" format sparse matrix data.  Read that data on
      /// Rank 0, and distribute it to all processors.  Return the
      /// resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 opens the
      ///   file and reads data from it, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param pNode [in] Kokkos Node object.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static sparse_matrix_ptr
      readSparseFile (const std::string& filename,
                      const RCP<const Comm<int> >& pComm,
                      const RCP<node_type>& pNode,
                      const bool callFillComplete=true,
                      const bool tolerant=false,
                      const bool debug=false)
      {
        const int myRank = Teuchos::rank (*pComm);
        std::ifstream in;

        // Only open the file on Rank 0.
        if (myRank == 0)
          in.open (filename.c_str());
        return readSparse (in, pComm, pNode, callFillComplete, tolerant, debug);
        // We can rely on the destructor of the input stream to close
        // the file on scope exit, even if readSparse() throws an
        // exception.
      }

      /// \brief Read sparse matrix from the given Matrix Market input stream.
      ///
      /// The given input stream need only be readable by MPI Rank 0
      /// (with respect to the given communicator).  The input stream
      /// should contain Matrix Market "coordinate" format sparse
      /// matrix data.  Read that data on Rank 0, and distribute it to
      /// all processors.  Return the resulting distributed CrsMatrix.
      ///
      /// \note This is a collective operation.  Only Rank 0 reads
      ///   data from the input stream, but all ranks participate and
      ///   wait for the final result.
      ///
      /// \param in [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param pNode [in] Kokkos Node object.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries
      ///   read in from the input stream.  (Not calling
      ///   fillComplete() may be useful if you want to change the
      ///   matrix after reading it from a file.)
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static sparse_matrix_ptr
      readSparse (std::istream& in,
                  const RCP<const Comm<int> >& pComm,
                  const RCP<node_type>& pNode,
                  const bool callFillComplete=true,
                  const bool tolerant=false,
                  const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::broadcast;
        using Teuchos::ptr;
        using Teuchos::reduceAll;
        using std::cerr;
        using std::endl;
        typedef Teuchos::ScalarTraits<scalar_type> STS;

        const bool extraDebug = false;
        const int myRank = pComm->getRank ();
        const int rootRank = 0;

        // Current line number in the input stream.  Various calls
        // will modify this depending on the number of lines that are
        // read from the input stream.  Only Rank 0 modifies this.
        size_t lineNumber = 1;

        if (debug && myRank == rootRank) {
          cerr << "Matrix Market reader: readSparse:" << endl
               << "-- Reading banner line" << endl;
        }

        // The "Banner" tells you whether the input stream represents
        // a sparse matrix, the symmetry type of the matrix, and the
        // type of the data it contains.
        //
        // pBanner will only be nonnull on MPI Rank 0.  It will be
        // null on all other MPI processes.
        RCP<const Banner> pBanner;
        {
          // We read and validate the Banner on Proc 0, but broadcast
          // the validation result to all processes.
          // Teuchos::broadcast doesn't currently work with bool, so
          // we use int (true -> 1, false -> 0).
          int bannerIsCorrect = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            // Read the Banner line from the input stream.
            try {
              pBanner = readBanner (in, lineNumber, tolerant, debug);
            }
            catch (std::exception& e) {
              errMsg << "Attempt to read the Matrix Market file's Banner line "
                "threw an exception: " << e.what();
              bannerIsCorrect = 0;
            }

            if (bannerIsCorrect) {
              // Validate the Banner for the case of a sparse matrix.
              // We validate on Proc 0, since it reads the Banner.

              // In intolerant mode, the matrix type must be "coordinate".
              if (! tolerant && pBanner->matrixType() != "coordinate") {
                bannerIsCorrect = 0;
                errMsg << "The Matrix Market input file must contain a "
                  "\"coordinate\"-format sparse matrix in order to create a "
                  "Tpetra::CrsMatrix object from it, but the file's matrix "
                  "type is \"" << pBanner->matrixType() << "\" instead.";
              }
              // In tolerant mode, we allow the matrix type to be
              // anything other than "array" (which would mean that
              // the file contains a dense matrix).
              if (tolerant && pBanner->matrixType() == "array") {
                bannerIsCorrect = 0;
                errMsg << "Matrix Market file must contain a \"coordinate\"-"
                  "format sparse matrix in order to create a Tpetra::CrsMatrix "
                  "object from it, but the file's matrix type is \"array\" "
                  "instead.  That probably means the file contains dense matrix "
                  "data.";
              }
            }
          } // Proc 0: Done reading the Banner, hopefully successfully.

          // Broadcast from Proc 0 whether the Banner was read correctly.
          broadcast (*pComm, rootRank, ptr (&bannerIsCorrect));

          // If the Banner is invalid, all processes throw an
          // exception.  Only Proc 0 gets the exception message, but
          // that's OK, since the main point is to "stop the world"
          // (rather than throw an exception on one process and leave
          // the others hanging).
          TEUCHOS_TEST_FOR_EXCEPTION(bannerIsCorrect == 0,
            std::invalid_argument, errMsg.str ());
        } // Done reading the Banner line and broadcasting success.
        if (debug && myRank == rootRank) {
          cerr << "-- Reading dimensions line" << endl;
        }

        // Read the matrix dimensions from the Matrix Market metadata.
        // dims = (numRows, numCols, numEntries).  Proc 0 does the
        // reading, but it broadcasts the results to all MPI
        // processes.  Thus, readCoordDims() is a collective
        // operation.  It does a collective check for correctness too.
        Tuple<global_ordinal_type, 3> dims =
          readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Making Adder for collecting matrix data" << endl;
        }

        // "Adder" object for collecting all the sparse matrix entries
        // from the input stream.  This is only nonnull on Proc 0.
        RCP<adder_type> pAdder =
          makeAdder (pComm, pBanner, dims, tolerant, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Reading matrix data" << endl;
        }
        //
        // Read the matrix entries from the input stream on Proc 0.
        //
        {
          // We use readSuccess to broadcast the results of the read
          // (succeeded or not) to all MPI processes.  Since
          // Teuchos::broadcast doesn't currently know how to send
          // bools, we convert to int (true -> 1, false -> 0).
          int readSuccess = 1;
          std::ostringstream errMsg; // Exception message (only valid on Proc 0)
          if (myRank == rootRank) {
            try {
              // Reader for "coordinate" format sparse matrix data.
              typedef Teuchos::MatrixMarket::CoordDataReader<adder_type,
                global_ordinal_type, scalar_type, STS::isComplex> reader_type;
              reader_type reader (pAdder);

              // Read the sparse matrix entries.
              std::pair<bool, std::vector<size_t> > results =
                reader.read (in, lineNumber, tolerant, debug);
              readSuccess = results.first ? 1 : 0;
            }
            catch (std::exception& e) {
              readSuccess = 0;
              errMsg << e.what();
            }
          }
          broadcast (*pComm, rootRank, ptr (&readSuccess));

          // It would be nice to add a "verbose" flag, so that in
          // tolerant mode, we could log any bad line number(s) on
          // Proc 0.  For now, we just throw if the read fails to
          // succeed.
          //
          // Question: If we're in tolerant mode, and if the read did
          // not succeed, should we attempt to call fillComplete()?
          TEUCHOS_TEST_FOR_EXCEPTION(readSuccess == 0, std::runtime_error,
            "Failed to read the Matrix Market sparse matrix file: "
            << errMsg.str());
        } // Done reading the matrix entries (stored on Proc 0 for now)

        if (debug && myRank == rootRank) {
          cerr << "-- Successfully read the Matrix Market data" << endl;
        }

        // In tolerant mode, we need to rebroadcast the matrix
        // dimensions, since they may be different after reading the
        // actual matrix data.  We only need to broadcast the number
        // of rows and columns.  Only Rank 0 needs to know the actual
        // global number of entries, since (a) we need to merge
        // duplicates on Rank 0 first anyway, and (b) when we
        // distribute the entries, each rank other than Rank 0 will
        // only need to know how many entries it owns, not the total
        // number of entries.
        if (tolerant) {
          if (debug && myRank == rootRank) {
            cerr << "-- Tolerant mode: rebroadcasting matrix dimensions"
                 << endl
                 << "----- Dimensions before: "
                 << dims[0] << " x " << dims[1]
                 << endl;
          }
          // Packed coordinate matrix dimensions (numRows, numCols).
          Tuple<global_ordinal_type, 2> updatedDims;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number
            // of rows is less than that specified in the
            // metadata.  We allow this case, and favor the
            // metadata so that the zero row(s) will be included.
            updatedDims[0] =
              std::max (dims[0], pAdder->getAdder()->numRows());
            updatedDims[1] = pAdder->getAdder()->numCols();
          }
          broadcast (*pComm, rootRank, updatedDims);
          dims[0] = updatedDims[0];
          dims[1] = updatedDims[1];
          if (debug && myRank == rootRank) {
            cerr << "----- Dimensions after: " << dims[0] << " x "
                 << dims[1] << endl;
          }
        }
        else {
          // In strict mode, we require that the matrix's metadata and
          // its actual data agree, at least somewhat.  In particular,
          // the number of rows must agree, since otherwise we cannot
          // distribute the matrix correctly.

          // Teuchos::broadcast() doesn't know how to broadcast bools,
          // so we use an int with the standard 1 == true, 0 == false
          // encoding.
          int dimsMatch = 1;
          if (myRank == rootRank) {
            // If one or more bottom rows of the matrix contain no
            // entries, then the Adder will report that the number of
            // rows is less than that specified in the metadata.  We
            // allow this case, and favor the metadata, but do not
            // allow the Adder to think there are more rows in the
            // matrix than the metadata says.
            if (dims[0] < pAdder->getAdder ()->numRows ()) {
              dimsMatch = 0;
            }
          }
          broadcast (*pComm, 0, ptr (&dimsMatch));
          if (dimsMatch == 0) {
            // We're in an error state anyway, so we might as well
            // work a little harder to print an informative error
            // message.
            //
            // Broadcast the Adder's idea of the matrix dimensions
            // from Proc 0 to all processes.
            Tuple<global_ordinal_type, 2> addersDims;
            if (myRank == rootRank) {
              addersDims[0] = pAdder->getAdder()->numRows();
              addersDims[1] = pAdder->getAdder()->numCols();
            }
            broadcast (*pComm, 0, addersDims);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dimsMatch == 0, std::runtime_error,
              "The matrix metadata says that the matrix is " << dims[0] << " x "
              << dims[1] << ", but the actual data says that the matrix is "
              << addersDims[0] << " x " << addersDims[1] << ".  That means the "
              "data includes more rows than reported in the metadata.  This "
              "is not allowed when parsing in strict mode.  Parse the matrix in "
              "tolerant mode to ignore the metadata when it disagrees with the "
              "data.");
          }
        } // Matrix dimensions (# rows, # cols, # entries) agree.

        if (debug && myRank == rootRank) {
          cerr << "-- Converting matrix data into CSR format on Proc 0" << endl;
        }

        // Now that we've read in all the matrix entries from the
        // input stream into the adder on Proc 0, post-process them
        // into CSR format (still on Proc 0).  This will facilitate
        // distributing them to all the processors.
        //
        // These arrays represent the global matrix data as a CSR
        // matrix (with numEntriesPerRow as redundant but convenient
        // metadata, since it's computable from rowPtr and vice
        // versa).  They are valid only on Proc 0.
        ArrayRCP<size_t> numEntriesPerRow;
        ArrayRCP<size_t> rowPtr;
        ArrayRCP<global_ordinal_type> colInd;
        ArrayRCP<scalar_type> values;

        // Proc 0 first merges duplicate entries, and then converts
        // the coordinate-format matrix data to CSR.
        {
          int mergeAndConvertSucceeded = 1;
          std::ostringstream errMsg;

          if (myRank == rootRank) {
            try {
              typedef Teuchos::MatrixMarket::Raw::Element<scalar_type,
                global_ordinal_type> element_type;
              typedef typename std::vector<element_type>::const_iterator iter_type;

              // Number of rows in the matrix.  If we are in tolerant
              // mode, we've already synchronized dims with the actual
              // matrix data.  If in strict mode, we should use dims
              // (as read from the file's metadata) rather than the
              // matrix data to determine the dimensions.  (The matrix
              // data will claim fewer rows than the metadata, if one
              // or more rows have no entries stored in the file.)
              const size_type numRows = dims[0];

              // Additively merge duplicate matrix entries.
              pAdder->getAdder()->merge ();

              // Get a temporary const view of the merged matrix entries.
              const std::vector<element_type>& entries =
                pAdder->getAdder()->getEntries();

              // Number of matrix entries (after merging).
              const size_t numEntries = (size_t)entries.size();

              if (debug) {
                cerr << "----- Proc 0: Matrix has numRows=" << numRows
                     << " rows and numEntries=" << numEntries
                     << " entries." << endl;
              }

              // Make space for the CSR matrix data.  Converting to
              // CSR is easier if we fill numEntriesPerRow with zeros
              // at first.
              numEntriesPerRow = arcp<size_t> (numRows);
              std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 0);
              rowPtr = arcp<size_t> (numRows+1);
              std::fill (rowPtr.begin(), rowPtr.end(), 0);
              colInd = arcp<global_ordinal_type> (numEntries);
              values = arcp<scalar_type> (numEntries);

              // Convert from array-of-structs coordinate format to CSR
              // (compressed sparse row) format.
              global_ordinal_type prvRow = 0;
              size_t curPos = 0;
              rowPtr[0] = 0;
              for (curPos = 0; curPos < numEntries; ++curPos) {
                const element_type& curEntry = entries[curPos];
                const global_ordinal_type curRow = curEntry.rowIndex();
                TEUCHOS_TEST_FOR_EXCEPTION(
                  curRow < prvRow, std::logic_error,
                  "Row indices are out of order, even though they are supposed "
                  "to be sorted.  curRow = " << curRow << ", prvRow = "
                  << prvRow << ", at curPos = " << curPos << ".  Please report "
                  "this bug to the Tpetra developers.");
                if (curRow > prvRow) {
                  for (global_ordinal_type r = prvRow+1; r <= curRow; ++r) {
                    rowPtr[r] = curPos;
                  }
                  prvRow = curRow;
                }
                numEntriesPerRow[curRow]++;
                colInd[curPos] = curEntry.colIndex();
                values[curPos] = curEntry.value();
              }
              // rowPtr has one more entry than numEntriesPerRow.  The
              // last entry of rowPtr is the number of entries in
              // colInd and values.
              rowPtr[numRows] = numEntries;
            } // Finished conversion to CSR format
            catch (std::exception& e) {
              mergeAndConvertSucceeded = 0;
              errMsg << "Failed to merge sparse matrix entries and convert to "
                "CSR format: " << e.what();
            }

            if (debug && mergeAndConvertSucceeded) {
              // Number of rows in the matrix.
              const size_type numRows = dims[0];
              const size_type maxToDisplay = 100;

              cerr << "----- Proc 0: numEntriesPerRow[0.."
                   << (numEntriesPerRow.size()-1) << "] ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > 0) {
                if (numRows > maxToDisplay) {
                  for (size_type k = 0; k < 2; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                  cerr << "... ";
                  for (size_type k = numRows-2; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                else {
                  for (size_type k = 0; k < numRows-1; ++k) {
                    cerr << numEntriesPerRow[k] << " ";
                  }
                }
                cerr << numEntriesPerRow[numRows-1];
              } // numRows > 0
              cerr << "]" << endl;

              cerr << "----- Proc 0: rowPtr ";
              if (numRows > maxToDisplay) {
                cerr << "(only showing first and last few entries) ";
              }
              cerr << "= [";
              if (numRows > maxToDisplay) {
                for (size_type k = 0; k < 2; ++k) {
                  cerr << rowPtr[k] << " ";
                }
                cerr << "... ";
                for (size_type k = numRows-2; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              else {
                for (size_type k = 0; k < numRows; ++k) {
                  cerr << rowPtr[k] << " ";
                }
              }
              cerr << rowPtr[numRows] << "]" << endl;
            }
          } // if myRank == rootRank
        } // Done converting sparse matrix data to CSR format

        // Now we're done with the Adder, so we can release the
        // reference ("free" it) to save space.  This only actually
        // does anything on Rank 0, since pAdder is null on all the
        // other MPI processes.
        pAdder = null;

        if (debug && myRank == rootRank) {
          cerr << "-- Making range, domain, and row maps" << endl;
        }

        // Make the maps that describe the matrix's range and domain,
        // and the distribution of its rows.  Creating a Map is a
        // collective operation, so we don't have to do a broadcast of
        // a success Boolean.
        map_ptr pRangeMap = makeRangeMap (pComm, pNode, dims[0]);
        map_ptr pDomainMap = makeDomainMap (pRangeMap, dims[0], dims[1]);
        map_ptr pRowMap = makeRowMap (null, pComm, pNode, dims[0]);

        if (debug && myRank == rootRank) {
          cerr << "-- Distributing the matrix data" << endl;
        }

        // Distribute the matrix data.  Each processor has to add the
        // rows that it owns.  If you try to make Proc 0 call
        // insertGlobalValues() for _all_ the rows, not just those it
        // owns, then fillComplete() will compute the number of
        // columns incorrectly.  That's why Proc 0 has to distribute
        // the matrix data and why we make all the processors (not
        // just Proc 0) call insertGlobalValues() on their own data.
        //
        // These arrays represent each processor's part of the matrix
        // data, in "CSR" format (sort of, since the row indices might
        // not be contiguous).
        ArrayRCP<size_t> myNumEntriesPerRow;
        ArrayRCP<size_t> myRowPtr;
        ArrayRCP<global_ordinal_type> myColInd;
        ArrayRCP<scalar_type> myValues;
        // Distribute the matrix data.  This is a collective operation.
        distribute (myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap,
                    numEntriesPerRow, rowPtr, colInd, values, debug);

        if (debug && myRank == rootRank) {
          cerr << "-- Inserting matrix entries on each processor";
          if (callFillComplete) {
            cerr << " and calling fillComplete()";
          }
          cerr << endl;
        }
        // Each processor inserts its part of the matrix data, and
        // then they all call fillComplete().  This method invalidates
        // the my* distributed matrix data before calling
        // fillComplete(), in order to save space.  In general, we
        // never store more than two copies of the matrix's entries in
        // memory at once, which is no worse than what Tpetra
        // promises.
        sparse_matrix_ptr pMatrix =
          makeMatrix (myNumEntriesPerRow, myRowPtr, myColInd, myValues,
                      pRowMap, pRangeMap, pDomainMap, callFillComplete);
        // Only use a reduce-all in debug mode to check if pMatrix is
        // null.  Otherwise, just throw an exception.  We never expect
        // a null pointer here, so we can save a communication.
        if (debug) {
          int localIsNull = pMatrix.is_null () ? 1 : 0;
          int globalIsNull = 0;
          reduceAll (*pComm, Teuchos::REDUCE_MAX, localIsNull, ptr (&globalIsNull));
          TEUCHOS_TEST_FOR_EXCEPTION(globalIsNull != 0, std::logic_error,
            "Reader::makeMatrix() returned a null pointer on at least one "
            "process.  Please report this bug to the Tpetra developers.");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(pMatrix.is_null(), std::logic_error,
            "Reader::makeMatrix() returned a null pointer.  "
            "Please report this bug to the Tpetra developers.");
        }

        // We can't get the dimensions of the matrix until after
        // fillComplete() is called.  Thus, we can't do the sanity
        // check (dimensions read from the Matrix Market data,
        // vs. dimensions reported by the CrsMatrix) unless the user
        // asked makeMatrix() to call fillComplete().
        //
        // Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
        // what one might think it does, so you have to ask the range
        // resp. domain map for the number of rows resp. columns.
        if (callFillComplete) {
          const int numProcs = pComm->getSize ();

          if (extraDebug && debug) {
            const global_size_t globalNumRows =
              pRangeMap->getGlobalNumElements();
            const global_size_t globalNumCols =
              pDomainMap->getGlobalNumElements();
            if (myRank == rootRank) {
              cerr << "-- Matrix is "
                   << globalNumRows << " x " << globalNumCols
                   << " with " << pMatrix->getGlobalNumEntries()
                   << " entries, and index base "
                   << pMatrix->getIndexBase() << "." << endl;
            }
            pComm->barrier ();
            for (int p = 0; p < numProcs; ++p) {
              if (myRank == p) {
                cerr << "-- Proc " << p << " owns "
                     << pMatrix->getNodeNumCols() << " columns, and "
                     << pMatrix->getNodeNumEntries() << " entries." << endl;
              }
              pComm->barrier ();
            }
          } // if (extraDebug && debug)
        } // if (callFillComplete)

        if (debug && myRank == rootRank) {
          cerr << "-- Done creating the CrsMatrix from the Matrix Market data"
               << endl;
        }
        return pMatrix;
      }

      /// \brief Read dense matrix (as a MultiVector) from the given
      ///   Matrix Market file.
      ///
      /// Open the given file on MPI Process 0 (with respect to the
      /// given communicator).  The file should contain Matrix Market
      /// "array" format dense matrix data.  Read that data on Process
      /// 0, and distribute it to all processors.  Return the
      /// resulting distributed MultiVector.
      ///
      /// See documentation of readDense() for details.
      ///
      /// \param filename [in] Name of the Matrix Market file from
      ///   which to read.  Both the filename and the file itself are
      ///   only accessed on Rank 0 of the given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the dense matrix will be distributed.
      /// \param node [in] Kokkos Node object.
      /// \param map [in/out] On input: if nonnull, the map describing
      ///   how to distribute the vector (not modified).  In this
      ///   case, the map's communicator and node must equal \c comm
      ///   resp. \c node.  If null on input, then on output, a
      ///   sensible (contiguous and uniformly distributed over the
      ///   given communicator) map describing the distribution of the
      ///   output multivector.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static RCP<multivector_type>
      readDenseFile (const std::string& filename,
                     const RCP<const comm_type>& comm,
                     const RCP<node_type>& node,
                     RCP<const map_type>& map,
                     const bool tolerant=false,
                     const bool debug=false)
      {
        std::ifstream in;
        if (comm->getRank () == 0) { // Only open the file on Proc 0.
          in.open (filename.c_str ()); // Destructor closes safely
	}
        return readDense (in, comm, node, map, tolerant, debug);
      }

      /// \brief Read dense matrix (as a MultiVector) from the given
      ///   Matrix Market input stream.
      ///
      /// The given input stream need only be readable by MPI Rank 0
      /// (with respect to the given communicator).  The input stream
      /// should contain Matrix Market "array" format dense matrix
      /// data.  Read that data on Process 0, and distribute it to all
      /// processes in the communicator.  Return the resulting
      /// distributed MultiVector.
      ///
      /// Unlike readSparse(), this method allows callers to supply a
      /// Map over which to distribute the resulting MultiVector.  The
      /// Map argument is optional; if null, we construct our own
      /// reasonable Map.  We let users supply their own Map, because
      /// a common case in Tpetra is to read in or construct a sparse
      /// matrix first, and then create dense (multi)vectors
      /// distributed with the sparse matrix's domain or range Map.
      ///
      /// \note This is a collective operation.  Only Process 0 in the
      ///   communicator opens the file and reads data from it, but
      ///   all processes in the communicator participate and wait for
      ///   the final result.
      ///
      /// \note "Tolerant" parsing mode means something different for
      ///   dense matrices than it does for sparse matrices.  Since
      ///   Matrix Market dense matrix files don't store indices with
      ///   each value, unlike sparse matrices, we can't determine the
      ///   matrix dimensions from the data alone.  Thus, we require
      ///   the metadata to include a valid number of rows and
      ///   columns.  "Tolerant" in the dense case refers to the data;
      ///   in tolerant mode, the number of stored matrix entries may
      ///   be more or less than the number reported by the metadata
      ///   (number of rows times number of columns).  If more, the
      ///   extra data are ignored; if less, the remainder is filled
      ///   in with zeros.
      ///
      /// \note On Map compatibility: Suppose that you write a
      ///   multivector X to a file.  Then, you read it back in as a
      ///   different multivector Y distributed over the same
      ///   communicator, but with a Map constructed by the input
      ///   routine (i.e., a null Map on input to readDenseFile() or
      ///   readDense()).  In that case, the only properties shared by
      ///   the maps of X and Y are that they have the same
      ///   communicator and the same number of GIDs.  The Maps need
      ///   not necessarily be compatible (in the sense of
      ///   Map::isCompatible()), and they certainly need not
      ///   necessarily be the same Map (in the sense of
      ///   Map::isSameAs()).
      ///
      /// \param in [in] The input stream from which to read.  The
      ///   stream is only accessed on Process 0 of the given
      ///   communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the dense matrix will be distributed.
      /// \param node [in] Kokkos Node object.
      /// \param map [in/out] On input: if nonnull, the map describing
      ///   how to distribute the vector (not modified).  In this
      ///   case, the map's communicator and node must equal comm
      ///   resp. node.  If null on input, then on output, a sensible
      ///   (contiguous and uniformly distributed over the given
      ///   communicator) map describing the distribution of the
      ///   output multivector.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the stream.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static RCP<multivector_type>
      readDense (std::istream& in,
                 const RCP<const comm_type>& comm,
                 const RCP<node_type>& node,
                 RCP<const map_type>& map,
                 const bool tolerant=false,
                 const bool debug=false)
      {
	return readDenseImpl<scalar_type> (in, comm, node, map, tolerant, debug);
      }

      /// \brief Read Map (as a MultiVector) from the given
      ///   Matrix Market file.
      ///
      /// Open the given file on MPI Process 0 (with respect to the
      /// given communicator).  The file should contain Matrix Market
      /// "array" format dense matrix data with two columns, as
      /// generated by Writer::writeMap() or Writer::writeMapFile().
      /// Read that data on Process 0, and distribute it to all
      /// processes.  Return the resulting Map.
      ///
      /// \param filename [in] Name of the Matrix Market file from
      ///   which to read.  Both the filename and the file itself are
      ///   only accessed on Process 0 of the given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the Map will be distributed.
      /// \param node [in] Kokkos Node object.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static RCP<const map_type> 
      readMapFile (const std::string& filename,
		   const RCP<const comm_type>& comm,
		   const RCP<node_type>& node,
		   const bool tolerant=false,
		   const bool debug=false)
      {
	std::ifstream in;
        if (comm->getRank () == 0) { // Only open the file on Proc 0.
          in.open (filename.c_str ()); // Destructor closes safely
	}
        return readMap (in, comm, node, tolerant, debug);
      }

    private:
      template<class MultiVectorScalarType>
      static RCP<Tpetra::MultiVector<MultiVectorScalarType, 
				     local_ordinal_type, 
				     global_ordinal_type, 
				     node_type> >
      readDenseImpl (std::istream& in,
		     const RCP<const comm_type>& comm,
		     const RCP<node_type>& node,
		     RCP<const map_type>& map,
		     const bool tolerant=false,
		     const bool debug=false)
      {
        using Teuchos::MatrixMarket::Banner;
        using Teuchos::MatrixMarket::checkCommentLine;
	using Teuchos::as;
	using Teuchos::broadcast;
	using Teuchos::outArg;
        using std::cerr;
        using std::endl;
        typedef MultiVectorScalarType ST;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;
        typedef Teuchos::ScalarTraits<ST> STS;
        typedef typename STS::magnitudeType MT;
        typedef Teuchos::ScalarTraits<MT> STM;
	typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;

        // Rank 0 is the only (MPI) process allowed to read from the
        // input stream.
        const int myRank = comm->getRank ();

        if (debug && myRank == 0) {
          cerr << "Matrix Market reader: readDenseImpl:" << endl;
        }

	// mfh 17 Feb 2013: It's not strictly necessary that the Comm
	// instances be identical and that the Node instances be
	// identical.  The essential condition is more complicated to
	// test and isn't the same for all Node types.  Thus, we just
	// leave it up to the user.

        // // If map is nonnull, check the precondition that its
        // // communicator resp. node equal comm resp. node.  Checking
        // // now avoids doing a lot of file reading before we detect the
        // // violated precondition.
        // TEUCHOS_TEST_FOR_EXCEPTION(
        //   ! map.is_null() && (map->getComm() != comm || map->getNode () != node,
        //   std::invalid_argument, "If you supply a nonnull Map, the Map's "
        //   "communicator and node must equal the supplied communicator resp. "
        //   "node.");

        // Process 0 will read in the matrix dimensions from the file,
        // and broadcast them to all ranks in the given communicator.
        // There are only 2 dimensions in the matrix, but we use the
        // third element of the Tuple to encode the banner's reported
        // data type: "real" == 0, "complex" == 1, and "integer" == 0
        // (same as "real").  We don't allow pattern matrices (i.e.,
        // graphs) since they only make sense for sparse data.
        Tuple<GO, 3> dims;
        dims[0] = 0;
        dims[1] = 0;

        // Current line number in the input stream.  Only valid on
        // Proc 0.  Various calls will modify this depending on the
        // number of lines that are read from the input stream.
        size_t lineNumber = 1;

        // Capture errors and their messages on Proc 0.
        std::ostringstream err;
	int localBannerReadSuccess = 1;
	int localDimsReadSuccess = 1;

        // Only Proc 0 gets to read matrix data from the input stream.
        if (myRank == 0) {
          if (debug && myRank == 0) {
            cerr << "-- Reading banner line (dense)" << endl;
          }

          // The "Banner" tells you whether the input stream
          // represents a dense matrix, the symmetry type of the
          // matrix, and the type of the data it contains.
	  RCP<const Banner> pBanner;
	  try {
	    pBanner = readBanner (in, lineNumber, tolerant, debug);
	  } catch (std::exception& e) {
	    err << e.what ();
	    localBannerReadSuccess = 0;
	  }
	  // Make sure the input stream is the right kind of data.
	  if (localBannerReadSuccess) {
	    if (pBanner->matrixType () != "array") {
	      err << "The Matrix Market file does not contain dense matrix "
		"data.  Its banner (first) line says that its matrix type is \""
		<< pBanner->matrixType () << "\", rather that the required "
		"\"array\".";
	      localBannerReadSuccess = 0;
	    } else if (pBanner->dataType() == "pattern") {
	      err << "The Matrix Market file's banner (first) "
		"line claims that the matrix's data type is \"pattern\".  This does "
		"not make sense for a dense matrix, yet the file reports the matrix "
		"as dense.  The only valid data types for a dense matrix are "
		"\"real\", \"complex\", and \"integer\".";
	      localBannerReadSuccess = 0;
	    } else {
	      // Encode the data type reported by the Banner as the
	      // third element of the dimensions Tuple.
	      dims[2] = encodeDataType (pBanner->dataType ());
	    }
	  } // if we successfully read the banner line

	  // At this point, we've successfully read the banner line.
	  // Now read the dimensions line.
	  if (localBannerReadSuccess) {
	    if (debug && myRank == 0) {
	      cerr << "-- Reading dimensions line (dense)" << endl;
	    }
	    // Keep reading lines from the input stream until we find
	    // a non-comment line, or until we run out of lines.  The
	    // latter is an error, since every "array" format Matrix
	    // Market file must have a dimensions line after the
	    // banner (even if the matrix has zero rows or columns, or
	    // zero entries).
	    std::string line;
	    bool commentLine = true;

	    while (commentLine) {
	      // Test whether it is even valid to read from the input
	      // stream wrapping the line.
	      if (in.eof () || in.fail ()) {
                err << "Unable to get array dimensions line (at all) from line "
		    << lineNumber << " of input stream.  The input stream "
		    << "claims that it is " 
		    << (in.eof() ? "at end-of-file." : "in a failed state.");
		localDimsReadSuccess = 0;
	      } else {
		// Try to get the next line from the input stream.
		if (getline (in, line)) {
		  ++lineNumber; // We did actually read a line.
		}
		// Is the current line a comment line?  Ignore start
		// and size; they are only useful for reading the
		// actual matrix entries.  (We could use them here as
		// an optimization, but we've chosen not to.)
		size_t start = 0, size = 0;
		commentLine = checkCommentLine (line, start, size, lineNumber, tolerant);
	      } // whether we failed to read the line at all
            } // while the line we just read is a comment line

	    //
	    // Get <numRows> <numCols> from the line we just read.
	    //
	    std::istringstream istr (line);

	    // Test whether it is even valid to read from the input
	    // stream wrapping the line.
            if (istr.eof () || istr.fail ()) {
	      err << "Unable to read any data from line " << lineNumber 
		  << " of input; the line should contain the matrix dimensions "
		  << "\"<numRows> <numCols>\".";
	      localDimsReadSuccess = 0;
	    } else { // It's valid to read from the line.
	      GO theNumRows = 0;
	      istr >> theNumRows; // Read in the number of rows.
	      if (istr.fail ()) {
		err << "Failed to get number of rows from line "
		    << lineNumber << " of input; the line should contains the "
		    << "matrix dimensions \"<numRows> <numCols>\".";
		localDimsReadSuccess = 0;
	      } else { // We successfully read the number of rows
		dims[0] = theNumRows; // Save the number of rows
		if (istr.eof ()) { // Do we still have data to read?
		  err << "No more data after number of rows on line "
		      << lineNumber << " of input; the line should contain the "
		      << "matrix dimensions \"<numRows> <numCols>\".";
		  localDimsReadSuccess = 0;
		} else { // Still data left to read; read in number of columns.
		  GO theNumCols = 0;
		  istr >> theNumCols; // Read in the number of columns
		  if (istr.fail ()) {
		    err << "Failed to get number of columns from line "
			<< lineNumber << " of input; the line should contain "
			<< "the matrix dimensions \"<numRows> <numCols>\".";
		    localDimsReadSuccess = 0;
		  } else { // We successfully read the number of columns
		    dims[1] = theNumCols; // Save the number of columns
		  } // if istr.fail ()
		} // if istr.eof ()
	      } // if we read the number of rows
	    } // if the input stream wrapping the dims line was (in)valid
	  } // if we successfully read the banner line
        } // if (myRank == 0)

	// Broadcast the matrix dimensions, the encoded data type, and
	// whether or not Proc 0 succeeded in reading the banner and
	// dimensions.
	Tuple<GO, 5> bannerDimsReadResult;
	if (myRank == 0) {
	  bannerDimsReadResult[0] = dims[0]; // numRows
	  bannerDimsReadResult[1] = dims[1]; // numCols
	  bannerDimsReadResult[2] = dims[2]; // encoded data type
	  bannerDimsReadResult[3] = localBannerReadSuccess;
	  bannerDimsReadResult[4] = localDimsReadSuccess;
	}
        // Broadcast matrix dimensions and the encoded data type from
        // Proc 0 to all the MPI processes.
        broadcast (*comm, 0, bannerDimsReadResult);

	TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[3] == 0, std::runtime_error,
	  "Failed to read banner line: " << err.str ());
	TEUCHOS_TEST_FOR_EXCEPTION(
          bannerDimsReadResult[4] == 0, std::runtime_error,
	  "Failed to read matrix dimensions line: " << err.str ());
	if (myRank != 0) {
	  dims[0] = bannerDimsReadResult[0];
	  dims[1] = bannerDimsReadResult[1];
	  dims[2] = bannerDimsReadResult[2];
	}

        // Tpetra objects want the matrix dimensions in these types.
        const global_size_t numRows = static_cast<global_size_t> (dims[0]);
        const size_t numCols = static_cast<size_t> (dims[1]);

        // Make a "Proc 0 owns everything" Map that we will use to
        // read in the multivector entries in the correct order on
        // Proc 0.  This must be a collective 
        RCP<const map_type> proc0Map; // "Proc 0 owns everything" Map
        if (map.is_null ()) {
          // The user didn't supply a Map.  Make a contiguous
          // distributed Map for them, using the read-in multivector
          // dimensions.
          map = createUniformContigMapWithNode<LO, GO, NT> (numRows, comm, node);
	  const size_t localNumRows = (myRank == 0) ? numRows : 0;
          proc0Map = createContigMapWithNode<LO, GO, NT> (numRows, localNumRows, comm, node);
        }
        else { // The user supplied a Map.
          proc0Map = computeGatherMap<map_type> (map);
        }

        // Make a multivector X owned entirely by Proc 0.
        RCP<MV> X = createMultiVector<ST, LO, GO, NT> (proc0Map, numCols);

        //
        // On Proc 0, read the Matrix Market data from the input
        // stream into the multivector X.
        //
	int localReadDataSuccess = 1;
        if (myRank == 0) {
	  try {
	    if (debug && myRank == 0) {
	      cerr << "-- Reading matrix data (dense)" << endl;
	    }

	    // Make sure that we can get a 1-D view of X.
            TEUCHOS_TEST_FOR_EXCEPTION(
              ! X->isConstantStride (), std::logic_error,
  	      "Can't get a 1-D view of the entries of the MultiVector X on "
	      "Process 0, because the stride between the columns of X is not "
	      "constant.  This shouldn't happen because we just created X and "
	      "haven't filled it in yet.  Please report this bug to the Tpetra "
	      "developers.");

	    // Get a writeable 1-D view of the entries of X.  Rank 0
	    // owns all of them.  The view will expire at the end of
	    // scope, so (if necessary) it will be written back to X
	    // at this time.
	    ArrayRCP<ST> X_view = X->get1dViewNonConst ();
            TEUCHOS_TEST_FOR_EXCEPTION(
              as<global_size_t> (X_view.size ()) < numRows * numCols,
              std::logic_error,
              "The view of X has size " << X_view << " which is not enough to "
              "accommodate the expected number of entries numRows*numCols = "
              << numRows << "*" << numCols << " = " << numRows*numCols << ".  "
              "Please report this bug to the Tpetra developers.");
	    const size_t stride = X->getStride ();

	    // The third element of the dimensions Tuple encodes the data
	    // type reported by the Banner: "real" == 0, "complex" == 1,
	    // "integer" == 0 (same as "real"), "pattern" == 2.  We do not
	    // allow dense matrices to be pattern matrices, so dims[2] ==
	    // 0 or 1.  We've already checked for this above.
	    const bool isComplex = (dims[2] == 1);
	    size_type count = 0, curRow = 0, curCol = 0;

	    std::string line;
	    while (getline (in, line)) {
	      ++lineNumber;
	      // Is the current line a comment line?  If it's not,
	      // line.substr(start,size) contains the data.
	      size_t start = 0, size = 0;
	      const bool commentLine =
		checkCommentLine (line, start, size, lineNumber, tolerant);
	      if (! commentLine) {
		// Make sure we have room in which to put the new matrix
		// entry.  We check this only after checking for a
		// comment line, because there may be one or more
		// comment lines at the end of the file.  In tolerant
		// mode, we simply ignore any extra data.
		if (count >= X_view.size()) {
		  if (tolerant) {
		    break;
		  }
		  else {
		    TEUCHOS_TEST_FOR_EXCEPTION(
                       count >= X_view.size(),
                       std::runtime_error,
                       "The Matrix Market input stream has more data in it than "
                       "its metadata reported.  Current line number is "
                       << lineNumber << ".");
		  }
		}

		// mfh 19 Dec 2012: Ignore everything up to the initial
		// colon.  writeDense() has the option to print out the
		// global row index in front of each entry, followed by
		// a colon and space.
		{
		  const size_t pos = line.substr (start, size).find (':');
		  if (pos != std::string::npos) {
		    start = pos+1;
		  }
		}
		std::istringstream istr (line.substr (start, size));
		// Does the line contain anything at all?  Can we
		// safely read from the input stream wrapping the
		// line?
		if (istr.eof() || istr.fail()) {
		  // In tolerant mode, simply ignore the line.
		  if (tolerant) {
		    break;
		  }
		  // We repeat the full test here so the exception
		  // message is more informative.
		  TEUCHOS_TEST_FOR_EXCEPTION(
                    ! tolerant && (istr.eof() || istr.fail()),
                    std::runtime_error, 
  		    "Line " << lineNumber << " of the Matrix Market file is "
                    "empty, or we cannot read from it for some other reason.");
		}
		// Current matrix entry to read in.
		ST val = STS::zero();
		// Real and imaginary parts of the current matrix entry.
		// The imaginary part is zero if the matrix is real-valued.
		MT real = STM::zero(), imag = STM::zero();

		// isComplex refers to the input stream's data, not to
		// the scalar type S.  It's OK to read real-valued
		// data into a matrix storing complex-valued data; in
		// that case, all entries' imaginary parts are zero.
		if (isComplex) {
		  // STS::real() and STS::imag() return a copy of
		  // their respective components, not a writeable
		  // reference.  Otherwise we could just assign to
		  // them using the istream extraction operator (>>).
		  // That's why we have separate magnitude type "real"
		  // and "imag" variables.

		  // Attempt to read the real part of the current entry.
		  istr >> real;
		  if (istr.fail()) {
		    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.eof(), std::runtime_error,
  		      "Failed to get the real part of a complex-valued matrix "
		      "entry from line " << lineNumber << " of the Matrix Market "
		      "file.");
		    // In tolerant mode, just skip bad lines.
		    if (tolerant) {
		      break;
		    }
		  } else if (istr.eof()) {
		    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.eof(), std::runtime_error,
		      "Missing imaginary part of a complex-valued matrix entry "
		      "on line " << lineNumber << " of the Matrix Market file.");
		    // In tolerant mode, let any missing imaginary part be 0.
		  } else {
		    // Attempt to read the imaginary part of the current
		    // matrix entry.
		    istr >> imag;
		    TEUCHOS_TEST_FOR_EXCEPTION(
                      ! tolerant && istr.fail(), std::runtime_error,
		      "Failed to get the imaginary part of a complex-valued "
		      "matrix entry from line " << lineNumber << " of the "
		      "Matrix Market file.");
		    // In tolerant mode, let any missing or corrupted
		    // imaginary part be 0.
		  }
		} else { // Matrix Market file contains real-valued data.
		  // Attempt to read the current matrix entry.
		  istr >> real;
		  TEUCHOS_TEST_FOR_EXCEPTION(
                    ! tolerant && istr.fail(), std::runtime_error,
		    "Failed to get a real-valued matrix entry from line "
		    << lineNumber << " of the Matrix Market file.");
		  // In tolerant mode, simply ignore the line if
		  // we failed to read a matrix entry.
		  if (istr.fail() && tolerant) {
		    break;
		  }
		}
		// In tolerant mode, we simply let pass through whatever
		// data we got.
		TEUCHOS_TEST_FOR_EXCEPTION(
                  ! tolerant && istr.fail(), std::runtime_error,
		  "Failed to read matrix data from line " << lineNumber 
		  << " of the Matrix Market file.");

		// Assign val = ST(real, imag).
		Teuchos::MatrixMarket::details::assignScalar<ST> (val, real, imag);

		curRow = count % numRows;
		curCol = count / numRows;
		X_view[curRow + curCol*stride] = val;
		++count;
	      } // if not a comment line
	    } // while there are still lines in the file, get the next one

	    TEUCHOS_TEST_FOR_EXCEPTION(
              ! tolerant && static_cast<global_size_t> (count) < numRows * numCols,
	      std::runtime_error,
	      "The Matrix Market metadata reports that the dense matrix is " 
	      << numRows <<  " x " << numCols << ", and thus has " 
	      << numRows*numCols << " total entries, but we only found " << count
	      << " entr" << (count == 1 ? "y" : "ies") << " in the file.");
	  } catch (std::exception& e) {
	    err << e.what ();
	    localReadDataSuccess = 0;
	  }
        } // if (myRank == 0)

	// Synchronize on whether Proc 0 successfully read the data.
	int globalReadDataSuccess = localReadDataSuccess;
        broadcast (*comm, 0, outArg (globalReadDataSuccess));
	TEUCHOS_TEST_FOR_EXCEPTION(
          globalReadDataSuccess == 0, std::runtime_error,
	  "Failed to read the multivector's data: " << err.str ());

        // If there's only one MPI process and the user didn't supply
        // a Map (i.e., pMap is null), we're done.  Set pMap to the
        // Map used to distribute X, and return X.
        if (comm->getSize () == 1 && map.is_null ()) {
          map = proc0Map;
          if (debug && myRank == 0) {
            cerr << "-- Done reading multivector" << endl;
          }
          return X;
        }

        if (debug && myRank == 0) {
          cerr << "-- Creating distributed Map and target MultiVector" << endl;
        }

        // Make a multivector Y with the distributed map pMap.
        RCP<MV> Y = createMultiVector<ST, LO, GO, NT> (map, numCols);

        // Make an Export object that will export X to Y.  First
        // argument is the source map, second argument is the target
        // map.
        Export<LO, GO, NT> exporter (proc0Map, map);

        if (debug && myRank == 0) {
          cerr << "-- Exporting from X (owned by Rank 0) to globally owned Y"
               << endl;
        }
        // Export X into Y.
        Y->doExport (*X, exporter, INSERT);

        if (debug && myRank == 0) {
          cerr << "-- Done reading multivector" << endl;
        }

        // Y is distributed over all process(es) in the communicator.
        return Y;
      }

    public:
      /// \brief Read Map (as a MultiVector) from the given input stream.
      ///
      /// Read the given input stream on MPI Process 0 (with respect
      /// to the given communicator).  The stream should contain
      /// Matrix Market "array" format dense matrix data with two
      /// columns, as generated by Writer::writeMap() or
      /// Writer::writeMapFile().  Distribute the data from Process 0
      /// to all processes.  Return the resulting Map.
      ///
      /// \param in [in/out] Input stream of Matrix Market data from
      ///   which to read.  This is only accessed on Process 0 of the
      ///   given communicator.
      /// \param comm [in] Communicator containing all process(es)
      ///   over which the Map will be distributed.
      /// \param node [in] Kokkos Node object.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static RCP<const map_type> 
      readMap (std::istream& in,
	       const RCP<const comm_type>& comm,
	       const RCP<node_type>& node,
	       const bool tolerant=false,
	       const bool debug=false)
      {
	using Tpetra::global_size_t;
	using Teuchos::Array;
	using Teuchos::ArrayRCP;
	using Teuchos::as;
	using Teuchos::Comm;
	using Teuchos::CommRequest;
	using Teuchos::outArg;
	using Teuchos::reduceAll;
	using Teuchos::REDUCE_MIN;
	using Teuchos::SerialComm;
	using Teuchos::toString;
	using std::cerr;
	using std::endl;
	typedef ptrdiff_t int_type; // Can hold int and GO
	typedef local_ordinal_type LO;
	typedef global_ordinal_type GO;
	typedef node_type NT;
	typedef Tpetra::MultiVector<int_type, LO, GO, NT> MV;

	const int numProcs = comm->getSize ();
	const int myRank = comm->getRank ();

	if (debug && myRank == 0) {
	  cerr << "readMap:" << endl
	       << "- Reading the Map as a dense multivector" << endl;
	}
	// This is currently the only place where we use 'tolerant'
	// and 'debug'.  Later, if we want to be clever, we could have
	// tolerant mode allow PIDs out of order.
	int localReadSuccess = 1;
	std::string readExMsg;

	RCP<MV> data; // Will only be valid on Proc 0
	if (myRank == 0) {
	  // If we want to reuse readDenseImpl, we have to make a
	  // communicator that only contains Proc 0.  Otherwise,
	  // readDenseImpl will redistribute the data to all
	  // processes.  While we eventually want that, neither we nor
	  // readDenseImpl know the correct Map to use at the moment.
	  // That depends on the second column of the multivector.
	  RCP<const Comm<int> > proc0Comm (new SerialComm<int> ());
	  try {
	    RCP<const map_type> dataMap;
	    data = readDenseImpl<int_type> (in, proc0Comm, node, dataMap, 
					    tolerant, debug);
	    (void) dataMap; // Silence "unused" warnings
	  } catch (std::exception& e) {
	    localReadSuccess = 0;
	    readExMsg = e.what ();
	  }
	}
	int globalReadSuccess = 1;
	reduceAll (*comm, REDUCE_MIN, localReadSuccess, 
		   outArg (globalReadSuccess));
	TEUCHOS_TEST_FOR_EXCEPTION(
          globalReadSuccess == 0, std::runtime_error,
	  "Tpetra::MatrixMarket::readMap: "
	  "Reading the Map failed with the following exception message: "
	  << readExMsg);

	if (debug && myRank == 0) {
	  cerr << "- Successfully read the dense multivector" << endl
	       << "- Extracting the Map's data on Proc 0" << endl;
	}

	ArrayRCP<const GO> myGids;
	GO indexBase = 0; // must be global min GID
	if (myRank == 0) {
	  // Assume that the Map's data are ordered by PID (2nd column).
	  ArrayRCP<const int_type> GIDs = data->getData (0);
	  ArrayRCP<const int_type> PIDs = data->getData (1);
	  TEUCHOS_TEST_FOR_EXCEPTION(
            GIDs.size () != PIDs.size (), std::logic_error, 
	    "GIDs.size() = " << GIDs.size() << " != PIDs.size() = " 
	    << PIDs.size() << ".  This should never happen.  "
	    "Please report this bug to the Tpetra developers.");				     
	  // Count of data in each process.
	  const size_type globalNumGIDs = GIDs.size ();
	  ArrayRCP<size_type> gidsPerProcess (numProcs, 0);
	  // Don't throw in the loop; just accumulate an error string.
	  // That way, we won't stall the other processes, which have
	  // posted receives.
	  std::ostringstream err;
	  // Error conditions.  If any are nonzero, there was an
	  // error.  We use int, not bool, because bool doesn't have
	  // an MPI_Datatype.
	  int numNegPids = 0;
	  int numTooBigPids = 0;
	  int pidsOutOfOrder = 0;
	  Array<size_type> badRows;
	  Array<size_type> startIndices (numProcs+1, 0);
	  startIndices[0] = 0;
	  int lastPid = 0;
	  for (size_type k = 0; k < globalNumGIDs; ++k) {
	    const GO gid = as<GO> (GIDs[k]);
	    const int pid = as<int> (PIDs[k]);
	    if (debug) {
	      cerr << "k: " << k << ", gid: " << gid << ", pid: " << pid 
		   << ", lastPid: " << lastPid << endl;
	    }
	    if (pid < 0) {
	      ++numNegPids;
	      badRows.push_back (k);
	      // err << "Row k=" << k << " (zero-based) of the Map data file, "
	      //   "corresponding to global index GID=" << gid << ", contains an "
	      //   "invalid negative process rank PID=" << pid << ".  This probably "
	      //   "means that the data file is corrupt or does not encode a Map.";
	    } else if (pid >= numProcs) {
	      ++numTooBigPids;
	      badRows.push_back (k);
	      // err << "Row k=" << k << " (zero-based) of the Map data file, "
	      //   "corresponding to global index GID=" << gid << ", contains a "
	      //   "process rank PID=" << pid << " >= the number of processes " 
	      //     << numProcs << " in the communicator.  This may mean that the "
	      //   "data file is corrupt or does not encode a Map.  It may also "
	      //   "mean that the communicator over which the Map was distributed "
	      //   "had a different process count than the given communicator.";
	    } else if (pid < lastPid) {
	      // Did the PID occur out of order? 
	      // We allow PIDs to have zero GIDs.
	      ++pidsOutOfOrder;
	      badRows.push_back (k);
	    } else { // We know now that pid is valid.
	      ++gidsPerProcess[pid];
	      if (k == 0 || GIDs[k] < indexBase) {
		indexBase = GIDs[k]; // indexBase must be the global min GID
	      }
	      // It could be that Proc 0 owns no GIDs.  In that case,
	      // PID[0] will be > 0 and the code below will fill in
	      // startIndices[p] = 0 for p = 1, 2, ..., PIDs[0].
	      if (pid > lastPid) {
		// startIndices is analogous to the 'ptr' array in CSR.
		// Fill in offsets for processes with no GIDs.
		for (int p = lastPid+1; p <= pid; ++p) {
		  startIndices[p] = k;
		}
	      }
	      lastPid = PIDs[k];
	    }
	  } // for each GID
	  startIndices[numProcs] = globalNumGIDs;

	  const bool readSucceeded = 	
    numNegPids == 0 && numTooBigPids == 0 && pidsOutOfOrder == 0;

	  if (debug && myRank == 0) {
	    if (readSucceeded) {
	      cerr << "- Proc 0: The Map's data are valid" << endl;
	    } else {
	      cerr << "- Proc 0: The Map's data are invalid" << endl;
	    }
	    cerr << "- startIndices: " << toString (startIndices) << endl
		 << "- gidsPerProcess: " << toString (gidsPerProcess ()) << endl;
	  }
	  if (readSucceeded) {
	    // We have to be tricky because we read in GIDs as
	    // int_type, which may differ from GO.  (For example, on
	    // 64-bit Mac or Linux, int_type is long; GO might be int
	    // or long.)  If GO == int_type, we can just make myGids a
	    // view of the relevant section of GIDs.  However, the
	    // code still has to compile when GO != int_type.  Since
	    // we know that GIDs will be valid until the end of this
	    // method, we can use a raw pointer in that case to avoid
	    // the copy.
	    if (typeid (GO) == typeid (int_type)) {
	      ArrayRCP<const int_type> myGidsAsIntType = 
		GIDs.persistingView (startIndices[0], gidsPerProcess[0]);
	      const int_type* myGidsRaw = myGidsAsIntType.getRawPtr ();
	      myGids = arcp<const GO> (reinterpret_cast<const GO*> (myGidsRaw),
				       0, gidsPerProcess[0], false);
	    } else {
	      ArrayRCP<const int_type> myGidsAsIntType = 
		GIDs.persistingView (startIndices[0], gidsPerProcess[0]);
	      ArrayRCP<GO> myGidsAsGO (myGidsAsIntType.size ());
	      for (size_type k = 0; k < myGidsAsIntType.size (); ++k) {
		myGidsAsGO[k] = as<GO> (myGidsAsIntType[k]);
	      }
	      myGids = myGidsAsGO.getConst ();
	    }
	  } else {
	    myGids = Teuchos::null;
	  }
	  // MPI guarantees ordering of messages, so we can use
	  // nonblocking sends for both (count, indexBase) and the
	  // content (the GIDs).  We don't send the GIDs if that
	  // process will own zero of them.  We use a count of -1 as a
	  // flag that something went wrong while reading.
	  Array<RCP<CommRequest<int> > > countRequests (numProcs-1);
	  for (int p = 1; p < numProcs; ++p) {
	    ArrayRCP<size_type> sendBuf (2);
	    if (readSucceeded) {
	      sendBuf[0] = gidsPerProcess[p];
	    } else { // read did not succeed
	      // Tell all the other processes that the read failed,
	      // before throwing an exception.
	      sendBuf[0] = as<size_type> (-1);
	    }
	    sendBuf[1] = as<size_type> (indexBase);
	    countRequests[p-1] = isend (*comm, sendBuf.getConst (), p);
          }
	  // This matches the "gidCount < 0" throw test on other processes.
	  TEUCHOS_TEST_FOR_EXCEPTION(! readSucceeded, std::runtime_error, 
            "Reading the Map from the given file on Process 0 failed.");

	  if (debug) {
	    cerr << "- Proc 0: indexBase = " << indexBase << endl;
	  }
	  Array<RCP<CommRequest<int> > > dataRequests;
	  ArrayRCP<GO> gidsToSend;
	  for (int p = 1; p < numProcs; ++p) {
	    if (gidsPerProcess[p] > 0) {
	      gidsToSend.resize (gidsPerProcess[p]);
	      for (size_type k = 0; k < gidsPerProcess[p]; ++k) {
		gidsToSend[k] = as<GO> (GIDs[startIndices[p] + k]);
	      }
	      dataRequests.push_back (isend (*comm, gidsToSend.getConst (), p));
	    }
	  }

	  waitAll (*comm, countRequests ());
	  waitAll (*comm, dataRequests ());

	} else { // if (myRank != 0)

	  const int rootRank = 0;

	  // Receive the count of GIDs to receive.
	  // If -1, there was an error.  If 0, don't do a second receive.
	  Tuple<size_type, 2> gidCountAndMinAllGid;
	  gidCountAndMinAllGid[0] = 0; // gidCount
	  receive (*comm, rootRank, 2, gidCountAndMinAllGid.getRawPtr ());
	  indexBase = as<GO> (gidCountAndMinAllGid[1]);
	  const size_type gidCount = gidCountAndMinAllGid[0];
	  indexBase = as<GO> (gidCountAndMinAllGid[1]);
	  if (gidCount < 0) {
	    // This matches the "! readSucceeded" throw test on Proc 0.
  	    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
              "Reading the Map from the given file on Process 0 failed.");
	  } else if (gidCount > 0) {
	    if (debug) {
	      std::ostringstream err;
	      err << "- Proc " << myRank << ": indexBase = " << indexBase << endl;
	      cerr << err.str ();
	    }
	    ArrayRCP<GO> myGidsCopy (gidCount);
	    receive (*comm, rootRank, as<int> (gidCount), myGidsCopy.getRawPtr ());
	    myGids = myGidsCopy.getConst ();
	  }
	} // whether myRank is 0

	if (debug) {
	  comm->barrier ();
	  if (myRank == 0) {
	    cerr << "- Successfully completed GID list distribution" << endl;
	  }
	  comm->barrier ();
	  for (int p = 0; p < numProcs; ++p) {
	    if (myRank == p) {
	      cerr << "-- Proc " << myRank << " GID list: " 
		   << toString (myGids ()) << endl;
	    }
	    comm->barrier ();
	    comm->barrier ();
	    comm->barrier ();
	  }
	}
	const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();
	return rcp (new map_type (INVALID, myGids (), indexBase, comm, node));
      }

    private:

      /// \brief Encode the Matrix Market data type as an int.
      ///
      /// We assume that the Banner has already validated the data
      /// type string, so the string is one of four canonical values.
      /// We encode these values as follows: "real" as 0, "complex" as
      /// 1, "integer" as 0 (same as "real", since the way we
      /// implement reading integer values is the same as the way we
      /// implement reading real values), and "pattern" as 2.
      ///
      /// We use this encoding for communicating the data type.
      static int
      encodeDataType (const std::string& dataType)
      {
        if (dataType == "real" || dataType == "integer") {
          return 0;
        } else if (dataType == "complex") {
          return 1;
        } else if (dataType == "pattern") {
          return 2;
        } else {
          // We should never get here, since Banner validates the
          // reported data type and ensures it is one of the accepted
          // values.
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Unrecognized Matrix Market data type \"" << dataType 
            << "\".  We should never get here.  "
            "Please report this bug to the Tpetra developers.");
        }
      }
    };

    /// \class Writer
    /// \brief Matrix Market file writer for CrsMatrix and MultiVector.
    /// \author Mark Hoemmen
    ///
    /// The Matrix Market (see their <a
    /// href="http://math.nist.gov/MatrixMarket"> web site </a> for
    /// details) defines a human-readable ASCII text file format for
    /// interchange of sparse and dense matrices.  This class defines
    /// methods for writing sparse and dense matrices to a Matrix
    /// Market file or input stream.
    ///
    /// All methods of this class assume that the file is only
    /// openable resp. the input stream is only writeable, on the MPI
    /// process with Rank 0 (with respect to the MPI communicator over
    /// which the given CrsMatrix or MultiVector is to be
    /// distributed).
    ///
    /// We define the MultiVector type accepted by \c writeDense() and
    /// \c writeDenseFile() using the scalar_type, local_ordinal_type,
    /// global_ordinal_type, and node_type typedefs in
    /// SparseMatrixType.  This ensures that the multivectors returned
    /// by those methods have a type compatible with the CrsMatrix
    /// sparse matrices accepted by \c writeSparse() and \c
    /// writeSparseFile().  We do this because the typical use case of
    /// Matrix Market files in Trilinos is to test sparse matrix
    /// methods, which usually involves reading a sparse matrix A and
    /// perhaps also a dense right-hand side b.  Also, this lets you
    /// use CrsMatrix objects with non-default LocalMatOps template
    /// parameters.  (If we templated on Scalar, LocalOrdinal,
    /// GlobalOrdinal, and Node, we would also have to template on
    /// LocalMatOps in order to deal with CrsMatrix types with
    /// nondefault LocalMatOps.  That would tie Writer to CrsMatrix
    /// anyway, since MultiVector is not templated on LocalMatOps.  As
    /// a result, we might as well just template on the CrsMatrix
    /// type, in order to use arbitrary LocalMatOps types without
    /// additional code.)
    ///
    /// \tparam SparseMatrixType A specialization of \c Tpetra::CrsMatrix.
    ///
    template<class SparseMatrixType>
    class Writer {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef RCP<sparse_matrix_type> sparse_matrix_ptr;

      /// \typedef scalar_type
      /// \brief Type of the entries of the sparse matrix.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// \typedef local_ordinal_type
      /// \brief Only used to define map_type.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// \typedef global_ordinal_type
      /// \brief Type of indices as read from the Matrix Market file.
      ///
      /// Indices of the sparse matrix are stored as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      /// \typedef node_type
      /// \brief The Kokkos Node type.
      typedef typename SparseMatrixType::node_type node_type;

      /// \typedef multivector_type
      /// \brief The MultiVector type associated with SparseMatrixType.
      typedef MultiVector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> multivector_type;
      /// \typedef map_type
      /// \brief Tpetra::Map specialization associated with SparseMatrixType.
      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

      /// \brief Print the sparse matrix in Matrix Market format, with comments.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to the given
      /// file, using the Matrix Market "coordinate" format.  MPI Proc
      /// 0 is the only MPI process that opens or writes to the file.
      /// Include the matrix name and description in the comments
      /// section of the file (after the initial banner line, but
      /// before the matrix metadata and data).
      ///
      /// \param filename [in] Name of the file to which to write the
      ///   given sparse matrix.  The matrix is distributed, but only
      ///   Proc 0 opens the file and writes to it.
      ///
      /// \param pMatrix [in] The sparse matrix to write to the file.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output file.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output file.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeSparseFile (const std::string& filename,
                       const RCP<const sparse_matrix_type>& pMatrix,
                       const std::string& matrixName,
                       const std::string& matrixDescription,
                       const bool debug=false)
      {
        const int myRank = Teuchos::rank (*(pMatrix->getComm()));
        std::ofstream out;

        // Only open the file on Rank 0.
        if (myRank == 0) out.open (filename.c_str());
        writeSparse (out, pMatrix, matrixName, matrixDescription, debug);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeSparse() throws an
        // exception.
      }

      /// \brief Print the sparse matrix in Matrix Market format.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to the given
      /// file, using the Matrix Market "coordinate" format.  MPI Proc
      /// 0 is the only MPI process that opens or writes to the file.
      ///
      /// \param filename [in] Name of the file to which to write the
      ///   given sparse matrix.  The matrix is distributed, but only
      ///   Proc 0 opens the file and writes to it.
      ///
      /// \param pMatrix [in] The sparse matrix to write to the file.
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeSparseFile (const std::string& filename,
                       const RCP<const sparse_matrix_type>& pMatrix,
                       const bool debug=false)
      {
        writeSparseFile (filename, pMatrix, "", "", debug);
      }

    public:

      /// \brief Print the sparse matrix in Matrix Market format, with comments.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to an output
      /// stream, using the Matrix Market "coordinate" format.  MPI
      /// Proc 0 is the only MPI process that writes to the output
      /// stream.
      ///
      /// \param out [out] Name of the output stream to which to write
      ///   the given sparse matrix.  The matrix is distributed, but
      ///   only Proc 0 writes to the output stream.
      ///
      /// \param pMatrix [in] The sparse matrix to write to the given
      ///   output stream.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeSparse (std::ostream& out,
                   const RCP<const sparse_matrix_type>& pMatrix,
                   const std::string& matrixName,
                   const std::string& matrixDescription,
                   const bool debug=false)
      {
        using std::cerr;
        using std::endl;
	typedef scalar_type ST;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef typename Teuchos::ScalarTraits<ST> STS;
        typedef typename STS::magnitudeType MT;
        typedef typename Teuchos::ScalarTraits<MT> STM;
        typedef typename ArrayView<ST>::size_type size_type;
	typedef typename ArrayView<const LO>::const_iterator lo_iter;
	typedef typename ArrayView<const GO>::const_iterator go_iter;
	typedef typename ArrayView<const ST>::const_iterator st_iter;

        // Make the output stream write floating-point numbers in
        // scientific notation.  It will politely put the output
        // stream back to its state on input, when this scope
        // terminates.
        Teuchos::MatrixMarket::details::SetScientific<ST> sci (out);

        RCP<const Comm<int> > comm = pMatrix->getComm();
        const int myRank = comm->getRank ();
        // Whether to print debugging output to stderr.
        const bool debugPrint = debug && myRank == 0;

	RCP<const map_type> rowMap = pMatrix->getRowMap ();
	RCP<const map_type> colMap = pMatrix->getColMap ();
	RCP<const map_type> domainMap = pMatrix->getDomainMap ();
	RCP<const map_type> rangeMap = pMatrix->getRangeMap ();

        const global_size_t numRows = rangeMap->getGlobalNumElements ();
        const global_size_t numCols = domainMap->getGlobalNumElements ();
        if (debugPrint) {
          cerr << "writeSparse:" << endl
               << "-- Input sparse matrix is:"
               << "---- " << numRows << " x " << numCols << " with "
               << pMatrix->getGlobalNumEntries() << " entries;" << endl
               << "---- "
               << (pMatrix->isGloballyIndexed() ? "Globally" : "Locally")
               << " indexed." << endl
	       << "---- Its row map has " << rowMap->getGlobalNumElements ()
               << " elements." << endl
               << "---- Its col map has " << colMap->getGlobalNumElements ()
               << " elements." << endl;
        }
        // Make the "gather" row map, where Proc 0 owns all rows and
        // the other procs own no rows.
        const size_t localNumRows = (myRank == 0) ? numRows : 0;
        RCP<node_type> node = rowMap->getNode();
        RCP<const map_type> gatherRowMap = 
	  rcp (new map_type (numRows, localNumRows, 
			     rowMap->getIndexBase (), 
			     comm, node));
        // Since the matrix may in general be non-square, we need to
        // make a column map as well.  In this case, the column map
        // contains all the columns of the original matrix, because we
        // are gathering the whole matrix onto Proc 0.
        const size_t localNumCols = (myRank == 0) ? numCols : 0;
        RCP<const map_type> gatherColMap =
	  rcp (new map_type (numCols, localNumCols, 
			     colMap->getIndexBase (), 
			     comm, node));
        // Current map is the source map, gather map is the target map.
        typedef Import<LO, GO, node_type> import_type;
        import_type importer (rowMap, gatherRowMap);

        // Create a new CrsMatrix to hold the result of the import.
        // The constructor needs a column map as well as a row map,
        // for the case that the matrix is not square.
        RCP<sparse_matrix_type> newMatrix =
          rcp (new sparse_matrix_type (gatherRowMap, gatherColMap,
                                       static_cast<size_t> (0)));
        // Import the sparse matrix onto Proc 0.
        newMatrix->doImport (*pMatrix, importer, INSERT);

        // fillComplete() needs the domain and range maps for the case
        // that the matrix is not square.
        {
          RCP<const map_type> gatherDomainMap =
	    rcp (new map_type (numCols, localNumCols, 
			       domainMap->getIndexBase (),
			       comm, node));
          RCP<const map_type> gatherRangeMap =
	    rcp (new map_type (numRows, localNumRows, 
			       rangeMap->getIndexBase (),
			       comm, node));
          newMatrix->fillComplete (gatherDomainMap, gatherRangeMap);
        }

        if (debugPrint) {
          cerr << "-- Output sparse matrix is:"
               << "---- " << newMatrix->getRangeMap ()->getGlobalNumElements ()
               << " x " 
	       << newMatrix->getDomainMap ()->getGlobalNumElements ()
               << " with "
               << newMatrix->getGlobalNumEntries () << " entries;" << endl
               << "---- "
               << (newMatrix->isGloballyIndexed () ? "Globally" : "Locally")
               << " indexed." << endl
               << "---- Its row map has "
               << newMatrix->getRowMap ()->getGlobalNumElements ()
               << " elements, with index base " 
	       << newMatrix->getRowMap ()->getIndexBase () << "." << endl
               << "---- Its col map has "
               << newMatrix->getColMap ()->getGlobalNumElements ()
               << " elements, with index base "
	       << newMatrix->getColMap ()->getIndexBase () << "." << endl
	       << "---- Element count of output matrix's column Map may differ "
	       << "from that of the input matrix's column Map, if some columns "
	       << "of the matrix contain no entries." << endl;
        }

        //
        // Print the metadata and the matrix entries on Rank 0.
        //
        if (myRank == 0) {
          // Print the Matrix Market banner line.  CrsMatrix stores
          // data nonsymmetrically ("general").  This implies that
          // readSparse() on a symmetrically stored input file,
          // followed by writeSparse() on the resulting sparse matrix,
          // will result in an output file with a different banner
          // line than the original input file.
          out << "%%MatrixMarket matrix coordinate "
              << (STS::isComplex ? "complex" : "real")
              << " general" << endl;

          // Print comments (the matrix name and / or description).
          if (matrixName != "") {
            printAsComment (out, matrixName);
	  }
          if (matrixDescription != "") {
            printAsComment (out, matrixDescription);
	  }

          // Print the Matrix Market header (# rows, # columns, #
          // nonzeros).  Use the range resp. domain map for the number
          // of rows resp. columns, since Tpetra::CrsMatrix uses the
          // column map for the number of columns.  That only
          // corresponds to the "linear-algebraic" number of columns
          // when the column map is uniquely owned (a.k.a. one-to-one),
          // which only happens if the matrix is (block) diagonal.
          out << newMatrix->getRangeMap ()->getGlobalNumElements () << " "
              << newMatrix->getDomainMap ()->getGlobalNumElements () << " "
              << newMatrix->getGlobalNumEntries () << endl;

	  RCP<const map_type> gatherColMap = newMatrix->getColMap ();
	  // The Matrix Market format expects one-based row and column
	  // indices.  We'll convert the indices on output from
	  // whatever index base they use to one-based indices.
          const GO rowIndexBase = gatherRowMap->getIndexBase ();
          const GO colIndexBase = gatherColMap->getIndexBase ();
          //
          // Print the entries of the matrix.
          //
          // newMatrix can never be globally indexed, since we called
          // fillComplete() on it.  We include code for both cases
          // (globally or locally indexed) just in case that ever
          // changes.
          if (newMatrix->isGloballyIndexed()) {
	    // We know that the "gather" row Map is contiguous, so we
	    // don't need to get the list of GIDs.
	    const GO minAllGlobalIndex = gatherRowMap->getMinAllGlobalIndex ();
	    const GO maxAllGlobalIndex = gatherRowMap->getMaxAllGlobalIndex ();
            for (GO globalRowIndex = minAllGlobalIndex; 
		 globalRowIndex <= maxAllGlobalIndex; // inclusive range
		 ++globalRowIndex) {
	      ArrayView<const GO> ind;
	      ArrayView<const ST> val;
	      newMatrix->getGlobalRowView (globalRowIndex, ind, val);
	      go_iter indIter = ind.begin ();
	      st_iter valIter = val.begin ();
	      for (; indIter != ind.end() && valIter != val.end();
		   ++indIter, ++valIter) {
		const GO globalColIndex = *indIter;
		// Convert row and column indices to 1-based.
		// This works because the global index type is signed.
		out << (globalRowIndex + 1 - rowIndexBase) << " "
		    << (globalColIndex + 1 - colIndexBase) << " ";
		if (STS::isComplex) {
		  out << STS::real (*valIter) << " " << STS::imag (*valIter);
		} else {
		  out << *valIter;
		}
		out << endl;
	      } // For each entry in the current row
	    } // For each row of the "gather" matrix
          } else { // newMatrix is locally indexed
	    typedef OrdinalTraits<GO> OTG;
	    for (LO localRowIndex = gatherRowMap->getMinLocalIndex();
		 localRowIndex <= gatherRowMap->getMaxLocalIndex();
		 ++localRowIndex) {
	      // Convert from local to global row index.
	      const GO globalRowIndex =
		gatherRowMap->getGlobalElement (localRowIndex);
	      TEUCHOS_TEST_FOR_EXCEPTION(
	        globalRowIndex == OTG::invalid(), std::logic_error,
		"Failed to convert the supposed local row index " 
		<< localRowIndex << " into a global row index.  "
		"Please report this bug to the Tpetra developers.");
	      ArrayView<const LO> ind;
	      ArrayView<const ST> val;
	      newMatrix->getLocalRowView (localRowIndex, ind, val);
	      lo_iter indIter = ind.begin ();
	      st_iter valIter = val.begin ();
	      for (; indIter != ind.end() && valIter != val.end();
		   ++indIter, ++valIter) {
		// Convert the column index from local to global.
		const GO globalColIndex = 
		  gatherColMap->getGlobalElement (*indIter);
		TEUCHOS_TEST_FOR_EXCEPTION(
                  globalColIndex == OTG::invalid(), std::logic_error,
		  "On local row " << localRowIndex << " of the sparse matrix: "
		  "Failed to convert the supposed local column index " 
		  << *indIter << " into a global column index.  Please report "
		  "this bug to the Tpetra developers.");
		// Convert row and column indices to 1-based.
		// This works because the global index type is signed.
		out << (globalRowIndex + 1 - rowIndexBase) << " "
		    << (globalColIndex + 1 - colIndexBase) << " ";
		if (STS::isComplex) {
		  out << STS::real (*valIter) << " " << STS::imag (*valIter);
		} else {
		  out << *valIter;
		}
		out << endl;
	      } // For each entry in the current row
	    } // For each row of the "gather" matrix
	  } // Whether the "gather" matrix is locally or globally indexed
        } // If my process' rank is 0
      }

      /// \brief Print the sparse matrix in Matrix Market format.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to an output
      /// stream, using the Matrix Market "coordinate" format.  MPI
      /// Proc 0 is the only MPI process that writes to the output
      /// stream.
      ///
      /// \param out [out] Name of the output stream to which to write
      ///   the given sparse matrix.  The matrix is distributed, but
      ///   only Proc 0 writes to the output stream.
      ///
      /// \param pMatrix [in] The sparse matrix to write to the given
      ///   output stream.
      ///
      /// \param debug [in] Whether to print possibly copious
      ///   debugging output to stderr on Proc 0.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeSparse (std::ostream& out,
                   const RCP<const sparse_matrix_type>& pMatrix,
                   const bool debug=false)
      {
        writeSparse (out, pMatrix, "", "", debug);
      }

      /// \brief Print the multivector in Matrix Market format, with comments.
      ///
      /// Write the given Tpetra::MultiVector matrix to the given
      /// file, using the Matrix Market "array" format for dense
      /// matrices.  MPI Proc 0 is the only MPI process that opens or
      /// writes to the file.
      ///
      /// \param filename [in] Name of the output file to create (on
      ///   MPI Proc 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output file.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output file.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output file.  If empty, we
      ///   don't print anything (not even an empty line).
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeDenseFile (const std::string& filename,
                      const RCP<const multivector_type>& X,
                      const std::string& matrixName,
                      const std::string& matrixDescription)
      {
        const int myRank = Teuchos::rank (*(X->getMap()->getComm()));
        std::ofstream out;

        if (myRank == 0) // Only open the file on Rank 0.
          out.open (filename.c_str());

        writeDense (out, X, matrixName, matrixDescription);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeDense() throws an
        // exception.
      }

      /// \brief Print the multivector in Matrix Market format.
      ///
      /// Write the given Tpetra::MultiVector matrix to the given
      /// file, using the Matrix Market "array" format for dense
      /// matrices.  MPI Proc 0 is the only MPI process that opens or
      /// writes to the file.
      ///
      /// \param filename [in] Name of the output file to create (on
      ///   MPI Proc 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output file.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeDenseFile (const std::string& filename,
                      const RCP<const multivector_type>& X)
      {
        writeDenseFile (filename, X, "", "");
      }

      /// \brief Print the multivector in Matrix Market format.
      ///
      /// Write the given Tpetra::MultiVector matrix to an output
      /// stream, using the Matrix Market "array" format for dense
      /// matrices.  MPI Proc 0 is the only MPI process that writes to
      /// the output stream.
      ///
      /// \param out [out] The output stream to which to write (on MPI
      ///   Proc 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output stream.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeDense (std::ostream& out,
                  const RCP<const multivector_type>& X,
                  const std::string& matrixName,
                  const std::string& matrixDescription)
      {
	writeDenseImpl<scalar_type> (out, *X, matrixName, matrixDescription);
      }

      /// \brief Print the multivector in Matrix Market format.
      ///
      /// Write the given Tpetra::MultiVector matrix to an output
      /// stream, using the Matrix Market "array" format for dense
      /// matrices.  MPI Proc 0 is the only MPI process that writes to
      /// the output stream.
      ///
      /// \param out [out] The output stream to which to write (on MPI
      ///   Proc 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output stream.
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto MPI Proc 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      ///
      static void
      writeDense (std::ostream& out,
                  const RCP<const multivector_type>& X)
      {
        writeDense (out, X, "", "");
      }

      //! Print the Map to the given output stream.
      static void
      writeMap (std::ostream& out, const map_type& map)
      {
	using Teuchos::ArrayRCP;
	using Teuchos::ArrayView;
	using Teuchos::as;
	using Teuchos::rcp_const_cast;
	using Teuchos::rcpFromRef;
	using Teuchos::TypeNameTraits;
	using std::cerr;
	using std::endl;
	typedef local_ordinal_type LO;
	typedef global_ordinal_type GO;
	typedef node_type NT;
	typedef int pid_type;

	// ptrdiff_t should be the biggest signed built-in integer
	// type that can hold any GO or pid_type (= int) quantity
	// without overflow.  Nevertheless, we'll test this.
	typedef ptrdiff_t int_type;
	TEUCHOS_TEST_FOR_EXCEPTION(
	  sizeof(GO) > sizeof(int_type), std::logic_error,
	  "The global ordinal type GO=" << TypeNameTraits<GO>::name () 
	  << " is too big for ptrdiff_t.  sizeof(GO) = " << sizeof (GO) 
	  << " > sizeof(ptrdiff_t) = " << sizeof (ptrdiff_t) << ".");
	TEUCHOS_TEST_FOR_EXCEPTION(
	  sizeof(pid_type) > sizeof(int_type), std::logic_error,
	  "The process rank type pid_type=" << TypeNameTraits<pid_type>::name ()
	  << " is too big for ptrdiff_t.  sizeof(pid_type) = " << sizeof (pid_type)
	  << " > sizeof(ptrdiff_t) = " << sizeof (ptrdiff_t) << ".");

	const int myRank = map.getComm ()->getRank ();
	const bool debug = false;

	// We pack the Map into a 2-column MultiVector.  Column 0
	// holds the GIDs, and Column 1 holds the PIDs.
	typedef Tpetra::MultiVector<int_type, LO, GO, NT> MV;
	MV data (rcpFromRef (map), 2);

	if (debug) {
	  if (myRank == 0) {
	    cerr << endl << "writeMap:" << endl << endl;
	  }
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	}

	// Column 0: GIDs for each process in sequence (by their
	// process ranks), in their original order on that process.
	// We put the GIDs first because the ordering of GIDs defines
	// the permutation that the Map represents, no matter how the
	// GIDs are distributed over processes.
	ArrayView<const GO> gidList = map.getNodeElementList ();
	{
	  ArrayRCP<int_type> V0 = data.getDataNonConst (0);
	  for (int_type k = 0; k < gidList.size (); ++k) {
	    V0[k] = as<int_type> (gidList[k]);
	  }
	}

	if (debug) {
	  if (myRank == 0) {
	    cerr << endl << "  - Done with Column 0" << endl << endl;
	  }
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	}

	// Column 1: PID (process rank) that owns the GID in the same
	// row.  The same GID might be owned by multiple processes; in
	// that case, it will appear multiple times in the first
	// column, with a different PID each time.
	{
	  const pid_type myRank = map.getComm ()->getRank ();
	  ArrayRCP<int_type> V1 = data.getDataNonConst (1);
	  for (int_type k = 0; k < gidList.size (); ++k) {
	    V1[k] = as<int_type> (myRank);
	  }
	}

	if (debug) {
	  if (myRank == 0) {
	    cerr << endl << "  - Done with Column 1" << endl << endl;
	  }
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	  map.getComm ()->barrier ();
	}

	// Map description explains each column.
	std::ostringstream desc; 
	desc << "This file encodes a Tpetra::Map." << endl
	     << "It is stored as a dense matrix with 2 columns:" << endl
	     << "  - Column 1: GID (global index)" << endl
	     << "  - Column 2: PID (rank) of process owning that GID" << endl;
	writeDenseImpl<int_type> (out, data, map.getObjectLabel (), desc.str ());
      }

      //! Write the Map to the given file.
      static void
      writeMapFile (const std::string& filename,
		    const map_type& map)
      {
        const int myRank = map.getComm ()->getRank ();
        std::ofstream out;
        if (myRank == 0) { // Only open the file on Proc 0.
          out.open (filename.c_str());
	}
        writeMap (out, map);
        // We can rely on the destructor of the output stream to close
        // the file on scope exit, even if writeDense() throws an
        // exception.
      }

    private:
      /// \brief Print the multivector in Matrix Market format.
      ///
      /// Write the given Tpetra::MultiVector matrix to an output
      /// stream, using the Matrix Market "array" format for dense
      /// matrices.  Process 0 is the only MPI process that writes to
      /// the output stream.
      ///
      /// \param out [out] The output stream to which to write (on
      ///   Process 0 only).
      ///
      /// \param X [in] The dense matrix (stored as a multivector) to
      ///   write to the output stream.
      ///
      /// \param matrixName [in] Name of the matrix, to print in the
      ///   comments section of the output stream.  If empty, we don't
      ///   print anything (not even an empty line).
      ///
      /// \param matrixDescription [in] Matrix description, to print
      ///   in the comments section of the output stream.  If empty,
      ///   we don't print anything (not even an empty line).
      ///
      /// \warning The current implementation gathers the whole matrix
      ///   onto Process 0.  This will cause out-of-memory errors if
      ///   the matrix is too big to fit on one process.  This will be
      ///   fixed in the future.
      template<class MultiVectorScalarType>
      static void
      writeDenseImpl (std::ostream& out,
		      const Tpetra::MultiVector<MultiVectorScalarType, local_ordinal_type, global_ordinal_type, node_type>& X,
		      const std::string& matrixName,
		      const std::string& matrixDescription)
      {
        using Tpetra::createOneToOne;
        using Tpetra::createNonContigMapWithNode;
        using Teuchos::as;
        using Teuchos::rcp;
        using Teuchos::rcp_const_cast;
        using Teuchos::rcp_dynamic_cast;
        using std::endl;

	typedef MultiVectorScalarType ST;
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;

        typedef typename Teuchos::ScalarTraits<ST> STS;
        typedef typename STS::magnitudeType MT;
        typedef typename Teuchos::ScalarTraits<MT> STM;
        typedef typename ArrayView<ST>::size_type size_type;
	typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;

        // Make the output stream write floating-point numbers in
        // scientific notation.  It will politely put the output
        // stream back to its state on input, when this scope
        // terminates.
        Teuchos::MatrixMarket::details::SetScientific<ST> sci (out);

        RCP<const Comm<int> > comm = X.getMap ()->getComm ();
        const int myRank = comm->getRank ();
        RCP<const map_type> map = X.getMap ();
        const global_size_t numRows = map->getGlobalNumElements ();
        // Promote to global_size_t.
        const global_size_t numCols = X.getNumVectors ();

	const bool debug = false;
	using std::cerr;
	if (debug && myRank == 0) {
	  cerr << "  writeDenseImpl:" << endl
	       << "  - Computing gather Map" << endl;
	}

        // Make the "gather" map, where Proc 0 owns all rows of X, and
        // the other procs own no rows.
        RCP<NT> node = map->getNode();
        RCP<const map_type> gatherMap = computeGatherMap<map_type> (map);

	if (debug && myRank == 0) {
	  cerr << "  - Computing Import" << endl;
	}

        // Create an Import object to import X's data into a
        // multivector Y owned entirely by Proc 0.  In the Import
        // constructor, X's map is the source map, and the "gather
        // map" is the target map.
        Import<LO, GO, NT> importer (map, gatherMap);

	if (debug && myRank == 0) {
	  cerr << "  - Creating target MultiVector" << endl;
	}

        // Create a new multivector Y to hold the result of the import.
        RCP<MV> Y = createMultiVector<ST, LO, GO, NT> (gatherMap, numCols);

	if (debug && myRank == 0) {
	  cerr << "  - Doing Import" << endl;
	}

        // Import the multivector onto Proc 0.
        Y->doImport (X, importer, INSERT);

        //
        // Print the matrix in Matrix Market format on Proc 0.
        //
	std::ostringstream err; // for exception messages on Proc 0
	int localWriteSuccess = 1; // whether the write succeeded on Proc 0
        if (myRank == 0) {
	  try {
	    std::string dataType;
	    if (STS::isComplex) {
	      dataType = "complex";
	    } else if (STS::isOrdinal) {
	      dataType = "integer";
	    } else {
	      dataType = "real";
	    }

	    // Print the Matrix Market banner line.  MultiVector
	    // stores data nonsymmetrically, hence "general".
	    out << "%%MatrixMarket matrix array " 
		<< dataType << " general" << endl;

	    // Print comments (the matrix name and / or description).
	    if (matrixName != "") {
	      printAsComment (out, matrixName);
	    }
	    if (matrixDescription != "") {
	      printAsComment (out, matrixDescription);
	    }
	    // Print the Matrix Market dimensions header for dense matrices.
	    out << numRows << " " << numCols << endl;

	    // Get a read-only view of the entries of Y.
	    // Rank 0 owns all of the entries of Y.
	    ArrayRCP<ArrayRCP<const ST> > Y_view = Y->get2dView ();
	    TEUCHOS_TEST_FOR_EXCEPTION(
              as<global_size_t> (Y_view.size ()) < numCols, std::logic_error,
	      "Y_view has size " << Y_view.size () << " < numCols = " 
	      << numCols << ".");

	    // Print the entries of the matrix, in column-major order.
	    for (global_size_t j = 0; j < numCols; ++j) {
	      ArrayRCP<const ST> Y_j = Y_view[j];
	      TEUCHOS_TEST_FOR_EXCEPTION(
                as<global_size_t> (Y_j.size ()) < numRows, std::logic_error,
	        "The current column's data has length " << Y_j.size () 
		<< " < numRows=" << numRows << ".");
	      for (global_size_t i = 0; i < numRows; ++i) {
		const ST Y_ij = Y_j[i];
		if (STS::isComplex) {
		  out << STS::real (Y_ij) << " " << STS::imag (Y_ij) << endl;
		} else {
		  out << Y_ij << endl;
		}
	      }
	    }
	  } catch (std::exception& e) {
	    err << e.what ();
	    localWriteSuccess = 0;
	  }
        } // if (myRank == 0)

	// Synchronize on error.  This is not necessary here, but it
	// will prevent later code from hanging in Proc 0 threw an
	// exception above.
	int globalWriteSuccess = localWriteSuccess;
	broadcast (*comm, 0, outArg (globalWriteSuccess));
	TEUCHOS_TEST_FOR_EXCEPTION(
          globalWriteSuccess == 0, std::runtime_error, 
	  "Failed to write data: " << err.str ());

	if (debug && myRank == 0) {
	  cerr << "  - Done" << endl;
	}
      }

      /// \brief Print the given possibly multiline string as a comment.
      ///
      /// If the string is empty, don't print anything (not even an
      /// empty line).  Otherwise, print each line of the string (one
      /// or more lines) as a comment in the comments section of a
      /// Matrix Market file (the part of the file after the initial
      /// banner line, but before the matrix's size metadata and
      /// data).
      ///
      /// \param out [out] The output string to which to print.  This
      ///   function is <i>not</i> a collective operation; whichever
      ///   MPI process calls it will print to the given output
      ///   stream.
      ///
      /// \param str [in] The string to print.  It consists of zero or
      ///   more lines.  If empty, nothing is printed, not even an
      ///   empty line.
      ///
      /// \note Printing comments is tricky only because the string
      ///   might contain newlines.  We have to ensure that all the
      ///   lines start with a comment character.  If they already do,
      ///   we print each line as is; otherwise, we append a comment
      ///   character and a space to each line.
      static void
      printAsComment (std::ostream& out, const std::string& str)
      {
        using std::endl;
        std::istringstream istream (str);
        std::string line;

        while (getline (istream, line)) {
          if (! line.empty()) {
            // Note that getline() doesn't store '\n', so we have to
            // append the endline ourselves.
            if (line[0] == '%') { // Line starts with a comment character.
              out << line << endl;
            }
            else { // Line doesn't start with a comment character.
              out << "%% " << line << endl;
            }
          }
        }
      }
    }; // class Writer

  } // namespace MatrixMarket
} // namespace Tpetra

#endif // __MatrixMarket_Tpetra_hpp
