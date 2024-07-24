// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_READTRIPLES_HPP
#define TPETRA_DETAILS_READTRIPLES_HPP

/// \file Tpetra_Details_ReadTriples.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::readAndDealOutTriples, which reads a Matrix
///   Market file or input stream on one process, and distributes the
///   resulting sparse matrix entries to the other processes.
///
/// \warning This is an implementation detail of Tpetra.
///   Users must not rely on this file or its contents.

#include "TpetraCore_config.h"
#include "Tpetra_Details_PackTriples.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_MatrixMarket_generic.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <iostream>
#include <typeinfo> // for debugging

namespace Tpetra {
namespace Details {

//
// Search for "SKIP DOWN TO HERE" (omit quotes) for the "public"
// interface.  I put "public" in quotes because it's public only for
// Tpetra developers, NOT for Tpetra users.
//

namespace Impl {

// mfh 01 Feb 2017: Unfortunately,
// Teuchos::MatrixMarket::readComplexData requires Teuchos to have
// complex arithmetic support enabled.  To avoid this issue, I
// reimplement the function here.  It's not very long.

/// \brief Read "<rowIndex> <colIndex> <realPart> <imagPart>" from a line.
///
/// Matrix Market files that store a sparse matrix with complex values
/// do so with one sparse matrix entry per line.  It is stored as
/// space-delimited ASCII text: the row index, the column index, the
/// real part, and the imaginary part, in that order.  Both the row
/// and column indices are 1-based.  This function attempts to read
/// one line from the given input stream istr, extract the row and
/// column indices and the real and imaginary parts, and write them to
/// the corresponding output variables.
///
/// \param istr [in/out] Input stream from which to attempt to read
///   one line.
/// \param rowIndex [out] On output: if successful, the row index read
///   from the line.
/// \param colIndex [out] On output: if successful, the column index
///   read from the line.
/// \param realPart [out] On output: if successful, the real part of
///   the matrix entry's value read from the line.
/// \param imagPart [out] On output: if successful, the imaginary part
///   of the matrix entry's value read from the line.
/// \param lineNumber [in] The current line number.  Used only for
///   diagnostic error messages.
/// \param tolerant [in] Whether to parse tolerantly.  In tolerant
///   mode, if this function fails in any way to read any of the data,
///   it will return false without throwing an exception.  Otherwise,
///   this function will either throw an exception or return true.
///
/// \return True if this function successfully read the line from istr
///   and extracted all the output data, false otherwise.  If
///   tolerant==false, this function never returns false; it either
///   returns true or throws an exception.
template<class OrdinalType, class RealType>
bool
readComplexData (std::istream& istr,
                 OrdinalType& rowIndex,
                 OrdinalType& colIndex,
                 RealType& realPart,
                 RealType& imagPart,
                 const std::size_t lineNumber,
                 const bool tolerant)
{
  using ::Teuchos::MatrixMarket::readRealData;

  RealType the_realPart, the_imagPart;
  if (! readRealData (istr, rowIndex, colIndex, the_realPart, lineNumber, tolerant)) {
    if (tolerant) {
      return false;
    }
    else {
      std::ostringstream os;
      os << "Failed to read pattern data and/or real value from line "
         << lineNumber << " of input";
      throw std::invalid_argument(os.str());
    }
  }
  if (istr.eof ()) {
    if (tolerant) {
      return false;
    }
    else {
      std::ostringstream os;
      os << "No more data after real value on line "
         << lineNumber << " of input";
      throw std::invalid_argument (os.str ());
    }
  }
  istr >> the_imagPart;
  if (istr.fail ()) {
    if (tolerant) {
      return false;
    }
    else {
      std::ostringstream os;
      os << "Failed to get imaginary value from line "
         << lineNumber << " of input";
      throw std::invalid_argument (os.str ());
    }
  }
  realPart = the_realPart;
  imagPart = the_imagPart;
  return true;
}


/// \brief Implementation of the readLine stand-alone function in this
///   namespace (see below).
///
/// Implementations are specialized on whether or not SC is a
/// complex-valued type.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
/// \tparam isComplex Whether SC is a complex-valued type.
template<class SC,
         class GO,
         const bool isComplex = ::Kokkos::ArithTraits<SC>::is_complex>
struct ReadLine {
  /// \brief Take a line from the Matrix Market file or input stream,
  ///   and process the sparse matrix entry in that line.
  ///
  /// \param processTriple [in] Closure, generally with side effects,
  ///   that takes in and stores off a sparse matrix entry.  First
  ///   argument is the (global) row index, second argument is the
  ///   (global) column index, and third argument is the value of the
  ///   entry.  The closure must NOT do MPI communication.  Return
  ///   value is an error code, that is zero if and only if the
  ///   closure succeeded.
  /// \param line [in] Current line of the Matrix Market file or input
  ///   stream to read.
  /// \param lineNumber [in] Current line number in the file or input
  ///   stream.
  /// \param tolerant [in] Whether to read tolerantly.
  /// \param errStrm [in] If not NULL, print any error messages to
  ///   this stream.
  /// \param debug [in] If true, print debug messages to \c *errStrm.
  ///
  /// \return Error code; 0 if and only if success.
  static int
  readLine (std::function<int (const GO, const GO, const SC&)> processTriple,
            const std::string& line,
            const std::size_t lineNumber,
            const bool tolerant = false,
            std::ostream* errStrm = NULL,
            const bool debug = false);
};

/// \brief Complex-arithmetic partial specialization of ReadLine.
///
/// This helps implement the readLine stand-alone function in this
/// namespace (see below).
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
template<class SC, class GO>
struct ReadLine<SC, GO, true> {
  /// \brief Take a line from the Matrix Market file or input stream,
  ///   and process the sparse matrix entry in that line.
  ///
  /// \param processTriple [in] Closure, generally with side effects,
  ///   that takes in and stores off a sparse matrix entry.  First
  ///   argument is the (global) row index, second argument is the
  ///   (global) column index, and third argument is the value of the
  ///   entry.  The closure must NOT do MPI communication.  Return
  ///   value is an error code, that is zero if and only if the
  ///   closure succeeded.
  /// \param line [in] Current line of the Matrix Market file or input
  ///   stream to read.
  /// \param lineNumber [in] Current line number in the file or input
  ///   stream.
  /// \param tolerant [in] Whether to read tolerantly.
  /// \param errStrm [in] If not NULL, print any error messages to
  ///   this stream.
  /// \param debug [in] If true, print debug messages to \c *errStrm.
  ///
  /// \return Error code; 0 if and only if success.
  static int
  readLine (std::function<int (const GO, const GO, const SC&)> processTriple,
            const std::string& line,
            const std::size_t lineNumber,
            const bool tolerant = false,
            std::ostream* errStrm = NULL,
            const bool debug = false)
  {
    using ::Teuchos::MatrixMarket::checkCommentLine;
    typedef typename ::Kokkos::ArithTraits<SC>::mag_type real_type;
    using std::endl;

    GO rowInd, colInd;
    real_type realPart, imagPart;
    std::istringstream istr (line);
    bool success = true;
    try {
      // Use the version of this function in this file, not the
      // version in Teuchos_MatrixMarket_generic.hpp, because the
      // latter only exists if HAVE_TEUCHOS_COMPLEX is defined.
      success = readComplexData (istr, rowInd, colInd, realPart, imagPart,
                                 lineNumber, tolerant);
    }
    catch (std::exception& e) {
      success = false;
      if (errStrm != NULL) {
        std::ostringstream os;
        os << "readLine: readComplexData threw an exception: " << e.what ()
           << endl;
        *errStrm << os.str ();
      }
    }

    if (success) {
      // if (debug && errStrm != NULL) {
      //   std::ostringstream os;
      //   os << "readLine: Got entry: row=" << rowInd << ", col=" << colInd
      //      << ", realPart=" << realPart << ", imagPart=" << imagPart
      //      << std::endl;
      //   *errStrm << os.str ();
      // }
      // This line may have side effects.
      const int errCode =
        processTriple (rowInd, colInd, SC (realPart, imagPart));
      if (errCode != 0 && errStrm != NULL) {
        std::ostringstream os;
        os << "readLine: processTriple returned " << errCode << " != 0."
           << endl;
        *errStrm << os.str ();
      }
      return errCode;
    }
    else {
      return -1;
    }
  }
};

/// \brief Real-arithmetic partial specialization of ReadLine.
///
/// This helps implement the readLine stand-alone function in this
/// namespace (see below).
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
template<class SC, class GO>
struct ReadLine<SC, GO, false> {
  /// \brief Take a line from the Matrix Market file or input stream,
  ///   and process the sparse matrix entry in that line.
  ///
  /// \param processTriple [in] Closure, generally with side effects,
  ///   that takes in and stores off a sparse matrix entry.  First
  ///   argument is the (global) row index, second argument is the
  ///   (global) column index, and third argument is the value of the
  ///   entry.  The closure must NOT do MPI communication.  Return
  ///   value is an error code, that is zero if and only if the
  ///   closure succeeded.
  /// \param line [in] Current line of the Matrix Market file or input
  ///   stream to read.
  /// \param lineNumber [in] Current line number in the file or input
  ///   stream.
  /// \param tolerant [in] Whether to read tolerantly.
  /// \param errStrm [in] If not NULL, print any error messages to
  ///   this stream.
  /// \param debug [in] If true, print debug messages to \c *errStrm.
  ///
  /// \return Error code; 0 if and only if success.
  static int
  readLine (std::function<int (const GO, const GO, const SC&)> processTriple,
            const std::string& line,
            const std::size_t lineNumber,
            const bool tolerant = false,
            std::ostream* errStrm = NULL,
            const bool debug = false)
  {
    using ::Teuchos::MatrixMarket::checkCommentLine;
    using ::Teuchos::MatrixMarket::readRealData;
    using std::endl;

    GO rowInd, colInd;
    SC val;
    std::istringstream istr (line);
    bool success = true;
    try {
      success = readRealData (istr, rowInd, colInd, val,
                              lineNumber, tolerant);
    }
    catch (std::exception& e) {
      success = false;
      if (errStrm != NULL) {
        std::ostringstream os;
        os << "readLine: readRealData threw an exception: " << e.what ()
           << endl;
        *errStrm << os.str ();
      }
    }

    if (success) {
      if (debug && errStrm != NULL) {
        std::ostringstream os;
        os << "readLine: Got entry: row=" << rowInd << ", col=" << colInd
           << ", val=" << val << std::endl;
        *errStrm << os.str ();
      }
      // This line may have side effects.
      const int errCode = processTriple (rowInd, colInd, val);
      if (errCode != 0 && errStrm != NULL) {
        std::ostringstream os;
        os << "readLine: processTriple returned " << errCode << " != 0."
           << endl;
        *errStrm << os.str ();
      }
      return errCode;
    }
    else {
      return -1;
    }
  }
};

/// \brief Take a line from the Matrix Market file or input stream,
///   and process the sparse matrix entry in that line.
///
/// The line must be a valid Matrix Market line, not a comment.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
///
/// \param processTriple [in] Closure, generally with side effects,
///   that takes in and stores off a sparse matrix entry.  First
///   argument is the (global) row index, second argument is the
///   (global) column index, and third argument is the value of the
///   entry.  The closure must NOT do MPI communication.  Return value
///   is an error code, that is zero if and only if the closure
///   succeeded.
/// \param line [in] The line from the Matrix Market file or input
///   stream to read.
/// \param lineNumber [in] Current line number in the file or input
///   stream.
/// \param tolerant [in] Whether to read tolerantly.
/// \param errStrm [in] If not NULL, print any error messages to this
///   stream.
/// \param debug [in] If true, print debug messages to \c *errStrm.
///
/// \return Error code; 0 if and only if success.
template<class SC, class GO>
int
readLine (std::function<int (const GO, const GO, const SC&)> processTriple,
          const std::string& line,
          const std::size_t lineNumber,
          const bool tolerant = false,
          std::ostream* errStrm = NULL,
          const bool debug = false)
{
  return ReadLine<SC, GO>::readLine (processTriple, line, lineNumber,
                                     tolerant, errStrm, debug);
}

/// \brief Read at most numTriplesToRead triples from the given Matrix
///   Market input stream, and pass along any resulting matrix entries
///   to the given closure.
///
/// The line must be a valid Matrix Market line, not a comment.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
///
/// \param inputStream [in/out] Input stream from which to read.
/// \param curLineNum [in/out] Current line number in the input stream.
/// \param numTriplesRead [out] On output: Number of matrix triples
///   (row index, column index, value) successfully read from the
///   input stream on this call to the function.x
/// \param processTriple [in] Closure, generally with side effects,
///   that takes in and stores off a sparse matrix entry.  First
///   argument is the (global) row index, second argument is the
///   (global) column index, and third argument is the value of the
///   entry.  The closure must NOT do MPI communication.  Return value
///   is an error code, that is zero if and only if the closure
///   succeeded.
/// \param maxNumTriplesToRead [in] Maximum number of triples to read
///   from the input stream on this call of the function.  This is a
///   strict upper bound for numTriplesRead (see above).
/// \param tolerant [in] Whether to read tolerantly.
/// \param errStrm [in] If not NULL, print any error messages to this
///   stream.
/// \param debug [in] If true, print debug messages to \c *errStrm.
///
/// \return Error code; 0 if and only if success.
template<class SC, class GO>
int
readTriples (std::istream& inputStream,
             std::size_t& curLineNum,
             std::size_t& numTriplesRead,
             std::function<int (const GO, const GO, const SC&)> processTriple,
             const std::size_t maxNumTriplesToRead,
             const bool tolerant = false,
             std::ostream* errStrm = NULL,
             const bool debug = false)
{
  using Teuchos::MatrixMarket::checkCommentLine;
  using std::endl;
  using std::size_t;

  numTriplesRead = 0; // output argument only
  if (inputStream.eof ()) {
    return 0; // no error, just nothing left to read
  }
  else if (inputStream.fail ()) {
    if (errStrm != NULL) {
      *errStrm << "Input stream reports a failure (not the same as "
        "end-of-file)." << endl;
    }
    return -1;
  }

  std::string line;
  std::vector<size_t> badLineNumbers;
  int errCode = 0; // 0 means success

  bool inputStreamCanStillBeRead = std::getline (inputStream, line).good ();
  ++curLineNum; // we read the line; we can't put it back
  while (inputStreamCanStillBeRead && numTriplesRead < maxNumTriplesToRead) {
    // if (debug && errStrm != NULL) {
    //   std::ostringstream os;
    //   os << "readTriples: Got line: \"" << line << "\"" << std::endl;
    //   *errStrm << os.str ();
    // }
    size_t start, size;

    const bool isCommentLine =
      checkCommentLine (line, start, size, curLineNum, tolerant);
    if (isCommentLine) {
      // Move on to the next line, if there is a next line.
      inputStreamCanStillBeRead = std::getline (inputStream, line).good ();
      ++curLineNum; // we read another line; we can't put it back
      continue; // move on to the next line
    }
    else { // not a comment line; should have a sparse matrix entry
      const std::string theLine = line.substr (start, size);
      // If the line has a valid sparse matrix entry, extract it and
      // hand it off to the processTriple closure.
      const int curErrCode =
        readLine (processTriple, theLine, curLineNum, tolerant, errStrm, debug);
      if (curErrCode != 0) {
        errCode = curErrCode;
        badLineNumbers.push_back (curLineNum);
      }
      else {
        ++numTriplesRead;
      }
      if (numTriplesRead < maxNumTriplesToRead) {
        inputStreamCanStillBeRead = std::getline (inputStream, line).good ();
      }
    }
  } // while there are lines to read and we need more triples

  if (errCode != 0 && errStrm != NULL) {
    const size_t numBadLines = badLineNumbers.size ();
    *errStrm << "Encountered " << numBadLines << " bad line"
             << (numBadLines != size_t (1) ? "s" : "")
             << ": [";
    for (size_t k = 0; k < numBadLines; ++k) {
      *errStrm << badLineNumbers[k];
      if (k + 1 < numBadLines) {
        *errStrm << ", ";
      }
    }
    *errStrm << "]" << endl;
  }
  if (! inputStream.eof () && inputStream.fail ()) {
    if (errCode == 0) {
      errCode = -1;
    }
    if (errStrm != NULL) {
      *errStrm << "The input stream is not at end-of-file, "
        "but is in a bad state." << endl;
    }
  }
  return errCode;
}

/// \brief Read at most maxNumEntPerMsg sparse matrix entries from the
///   input stream, and send them to the process with rank destRank.
///
/// To be called only by the sending process.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
///
/// \param inputStream [in/out] Input stream from which to read.
/// \param curLineNum [in/out] Current line number in the input stream.
/// \param numEntRead [out] Number of matrix entries successfully read
///   from the input stream on this call of the function.
/// \param sizeBuf [in/out] Array of length 1, for sending message size.
/// \param msgBuf [in/out] Message buffer; to be resized as needed.
/// \param rowInds [out] Row indices read from the file.
/// \param colInds [out] Column indices read from the file.
/// \param vals [out] Matrix values read from the file.
/// \param maxNumEntPerMsg [in] Maximum number of sparse matrix
///   entries to read from the input stream on this call to the
///   function.
/// \param destRank [in] Rank of the process where to send the triples.
/// \param comm [in] Communicator to use for sending the triples.
/// \param tolerant [in] Whether to read tolerantly.
/// \param errStrm [in] If not NULL, print any error messages to this
///   stream.
template<class SC, class GO>
int
readAndSendOneBatchOfTriples (std::istream& inputStream,
                              std::size_t& curLineNum,
                              std::size_t& numEntRead,
                              ::Teuchos::ArrayRCP<int>& sizeBuf,
                              ::Teuchos::ArrayRCP<char>& msgBuf,
                              std::vector<GO>& rowInds,
                              std::vector<GO>& colInds,
                              std::vector<SC>& vals,
                              const std::size_t maxNumEntPerMsg,
                              const int destRank,
                              const ::Teuchos::Comm<int>& comm,
                              const bool tolerant = false,
                              std::ostream* errStrm = NULL,
                              const bool debug = false)
{
  using ::Tpetra::Details::countPackTriplesCount;
  using ::Tpetra::Details::countPackTriples;
  using ::Tpetra::Details::packTriplesCount;
  using ::Tpetra::Details::packTriples;
  using ::Teuchos::isend;
  using std::endl;

  using ::Kokkos::ArithTraits;
  // constexpr int sizeTag = 42 + (ArithTraits<SC>::is_complex ? 100 : 0);
  // constexpr int msgTag = 43 + (ArithTraits<SC>::is_complex ? 100 : 0);
  constexpr int sizeTag = 42;
  constexpr int msgTag = 43;
  //constexpr int srcRank = 0;
  int errCode = 0;

  // This doesn't actually deallocate memory; it just changes the size
  // back to zero, so that push_back starts over from the beginning.
  rowInds.resize (0);
  colInds.resize (0);
  vals.resize (0);
  // Closure that adds the new matrix entry to the above temp arrays.
  auto processTriple = [&rowInds, &colInds, &vals]
    (const GO rowInd, const GO colInd, const SC& val) {
      try {
        rowInds.push_back (rowInd);
        colInds.push_back (colInd);
        vals.push_back (val);
      }
      catch (...) {
        return -1;
      }
      return 0;
    };
  numEntRead = 0; // output argument
  errCode = readTriples<SC, GO> (inputStream, curLineNum, numEntRead,
                                 processTriple, maxNumEntPerMsg, tolerant,
                                 errStrm, debug);
  if (debug && errStrm != NULL) {
    std::ostringstream os;
    os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
       << ", GO=" << typeid (GO).name () << ": "
       << "readAndSendOneBatchOfTriples: readTriples read "
       << numEntRead << " matrix entries, and returned errCode="
       << errCode << "." << std::endl;
    *errStrm << os.str ();
  }
  if (numEntRead != rowInds.size () ||
      numEntRead != colInds.size () ||
      numEntRead != vals.size ()) {
    if (errStrm != NULL) {
      *errStrm << "readTriples size results are not consistent.  "
               << "numEntRead = " << numEntRead
               << ", rowInds.size() = " << rowInds.size ()
               << ", colInds.size() = " << colInds.size ()
               << ", and vals.size() = " << vals.size () << "."
               << std::endl;
    }
    if (errCode == 0) {
      errCode = -1;
    }
  }

  // We don't consider reading having "failed" if we've reached
  // end-of-file before reading maxNumEntPerMsg entries.  It's OK if
  // we got fewer triples than that.  Furthermore, we have to send at
  // least one message to the destination process, even if the read
  // from the file failed.

  if (numEntRead == 0 || errCode != 0) {
    // Send a message size of zero to the receiving process, to tell
    // it that we have no triples to send, or that there was an error
    // reading.  The latter just means that we "go through the
    // motions," then broadcast the error code.
    sizeBuf[0] = 0;
    if (debug && errStrm != NULL) {
      std::ostringstream os;
      os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
         << ", GO=" << typeid (GO).name () << ": "
         << "Post send (size=0, errCode=" << errCode << ") "
         << "to " << destRank << " with tag " << sizeTag << endl;
      *errStrm << os.str ();
    }
    send (sizeBuf.getRawPtr (), 1, destRank, sizeTag, comm);
    return errCode;
  }
  else { // we read a nonzero # of triples, without error
    const int numEnt = static_cast<int> (numEntRead);
    int countSize = 0; // output argument
    int triplesSize = 0; // output argument

    errCode = countPackTriplesCount (comm, countSize, errStrm);
    // countSize should never be nonpositive, since we have to pack an
    // integer size.
    if (countSize <= 0 && errCode == 0) {
      errCode = -1;
    }

    if (errCode != 0) {
      // Send zero to the receiving process, to tell it about the error.
      sizeBuf[0] = 0;
      if (debug && errStrm != NULL) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
           << ", GO=" << typeid (GO).name () << ": "
           << "Post send (size=0, error case) to " << destRank
           << " with tag " << sizeTag << endl;
        *errStrm << os.str ();
      }
      send (sizeBuf.getRawPtr (), 1, destRank, sizeTag, comm);
      return errCode;
    }
    else { // countPackTriplesCount succeeded
      errCode = countPackTriples<SC, GO> (numEnt, comm, triplesSize, errStrm);
      if (errCode != 0) {
        // Send a message size of zero to the receiving process, to
        // tell it that there was an error counting.
        sizeBuf[0] = 0;
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Post send (size=0, error case) to " << destRank
             << " with tag " << sizeTag << endl;
          *errStrm << os.str ();
        }
        send (sizeBuf.getRawPtr (), 1, destRank, sizeTag, comm);
        return errCode;
      }
      else { // countPackTriples succeeded; message packed & ready to send
        // Send the message size (in bytes).  We can use a nonblocking
        // send here, and try to overlap with message packing.
        const int outBufSize = countSize + triplesSize;
        sizeBuf[0] = outBufSize;
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Post isend (size=" << sizeBuf[0] << ") to " << destRank
             << " with tag " << sizeTag << endl;
          *errStrm << os.str ();
        }
        auto sizeReq = isend<int, int> (sizeBuf, destRank, sizeTag, comm);

        msgBuf.resize (outBufSize);
        char* outBuf = msgBuf.getRawPtr ();

        // If anything goes wrong with packing, send the pack buffer
        // anyway, since the receiving process expects a message.
        int outBufCurPos = 0; // input/output argument
        errCode = packTriplesCount (numEnt, outBuf, outBufSize,
                                    outBufCurPos, comm, errStrm);
        if (errCode == 0) {
          errCode = packTriples<SC, GO> (rowInds.data (), colInds.data (),
                                         vals.data (), numEnt, outBuf,
                                         outBufSize, outBufCurPos, comm,
                                         errStrm);
        }
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Post isend (packed data) to " << destRank
             << " with tag " << msgTag << endl;
          *errStrm << os.str ();
        }
        auto msgReq = isend<int, char> (msgBuf, destRank, msgTag, comm);

        // Wait on the two messages.  It doesn't matter in what order
        // we send them, because they have different tags.  The
        // receiving process will wait on the first message first, in
        // order to get the size of the second message.
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Wait on isend (size)" << endl;
          *errStrm << os.str ();
        }
        sizeReq->wait ();
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Wait on isend (packed data)" << endl;
          *errStrm << os.str ();
        }
        msgReq->wait ();

        // This doesn't actually deallocate; it just resets sizes to zero.
        rowInds.clear ();
        colInds.clear ();
        vals.clear ();
      }
    }
  }
  return errCode;
}

/// \brief Read at most maxNumEntPerMsg sparse matrix entries from the
///   input stream, and send them to the process with rank destRank.
///
/// To be called only by the sending process.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.

/// \tparam CommRequestPtr The type of a (smart) pointer to a
///   "communication request" returned by, e.g., ::Teuchos::ireceive.
///   It must implement operator= and operator->, and the thing to
///   which it points must implement <tt>void wait()</tt>.  A model
///   for this is ::Teuchos::RCP< ::Teuchos::CommRequest<int> >.
///
/// \param rowInds [out] Row indices to receive.  Will be resized as
///   needed.
/// \param colInds [out] Column indices to receive.  Will be resized
///   as needed.
/// \param vals [out] Matrix values to receive.  Will be resized as
///   needed.
/// \param numEnt [out] Number of matrix entries (triples) received.
/// \param sizeBuf [in/out] Array of length 1, for receiving message
///   size (size in bytes of \c msgBuf).
/// \param msgBuf [in/out] Message buffer; to be resized as needed.
/// \param sizeReq [in/out] Preposed receive request for message size.
///   After waiting on this, you may read the contents of sizeBuf.
///   This is a nonconst (smart) pointer reference, so that we can
///   assign to it.  A model for this is
///   ::Teuchos::RCP< ::Teuchos::CommRequest<int> >&.
/// \param srcRank [in] Rank of the process from which to receive the
///   matrix entries (triples).
/// \param comm [in] Communicator to use for receiving the triples.
/// \param tolerant [in] Whether to read tolerantly.
/// \param errStrm [in] If not NULL, print any error messages to this
///   stream.
template<class SC, class GO, class CommRequestPtr>
int
recvOneBatchOfTriples (std::vector<GO>& rowInds,
                       std::vector<GO>& colInds,
                       std::vector<SC>& vals,
                       int& numEnt,
                       ::Teuchos::ArrayRCP<int>& sizeBuf,
                       ::Teuchos::ArrayRCP<char>& msgBuf,
                       CommRequestPtr& sizeReq,
                       const int srcRank,
                       const ::Teuchos::Comm<int>& comm,
                       const bool tolerant = false,
                       std::ostream* errStrm = NULL,
                       const bool debug = false)
{
  using ::Tpetra::Details::unpackTriplesCount;
  using ::Tpetra::Details::unpackTriples;
  using ::Kokkos::ArithTraits;

  ////constexpr int sizeTag = 42 + (ArithTraits<SC>::is_complex ? 100 : 0);
  //constexpr int msgTag = 43 + (ArithTraits<SC>::is_complex ? 100 : 0);
  //constexpr int sizeTag = 42;
  constexpr int msgTag = 43;
  int errCode = 0; // return value
  numEnt = 0; // output argument

  // Wait on the ireceive we preposted before calling this function.
  if (debug && errStrm != NULL) {
    std::ostringstream os;
    os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
       << ", GO=" << typeid (GO).name () << ": "
       << "Wait on irecv (size)" << std::endl;
    *errStrm << os.str ();
  }
  sizeReq->wait ();
  sizeReq = CommRequestPtr (NULL);
  const int inBufSize = sizeBuf[0];
  if (debug && errStrm != NULL) {
    std::ostringstream os;
    os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
       << ", GO=" << typeid (GO).name () << ": "
       << "Received size: sizeBuf[0]=" << sizeBuf[0] << std::endl;
    *errStrm << os.str ();
  }

  if (inBufSize == 0) {
    numEnt = 0;
    rowInds.resize (0);
    colInds.resize (0);
    vals.resize (0);
  }
  else {
    msgBuf.resize (inBufSize);
    char* inBuf = msgBuf.getRawPtr ();

    if (debug && errStrm != NULL) {
      std::ostringstream os;
      os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
         << ", GO=" << typeid (GO).name () << ": "
         << "Post irecv (packed data) " << "from " << srcRank
         << " with tag " << msgTag << std::endl;
      *errStrm << os.str ();
    }
    auto msgReq = ::Teuchos::ireceive (msgBuf, srcRank, msgTag, comm);
    if (debug && errStrm != NULL) {
      std::ostringstream os;
      os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
         << ", GO=" << typeid (GO).name () << ": "
         << "Wait on irecv (packed data)" << std::endl;
      *errStrm << os.str ();
    }
    msgReq->wait ();

    int inBufCurPos = 0; // output argument
    errCode = unpackTriplesCount (inBuf, inBufSize, inBufCurPos,
                                  numEnt, comm, errStrm);
    if (errCode == 0) {
      rowInds.resize (numEnt);
      colInds.resize (numEnt);
      vals.resize (numEnt);
      errCode = unpackTriples<SC, GO> (inBuf, inBufSize, inBufCurPos,
                                       rowInds.data (), colInds.data (),
                                       vals.data (), numEnt, comm, errStrm);
    }
  }
  return errCode;
}

} // namespace Impl

//
// SKIP DOWN TO HERE FOR "PUBLIC" INTERFACE
//

/// \brief On Process 0 in the given communicator, read sparse matrix
///   entries (in chunks of at most maxNumEntPerMsg entries at a time)
///   from the input stream, and "deal them out" to all other
///   processes in the communicator.
///
/// This is a collective over the communicator.
///
/// \tparam SC The type of the value of each matrix entry.
/// \tparam GO The type of each (global) index of each matrix entry.
///
/// \param inputStream [in/out] Input stream from which to read Matrix
///   Market - format matrix entries ("triples").  Only Process 0 in
///   the communicator needs to be able to access this.
/// \param curLineNum [in/out] On both input and output, the
///   current line number in the input stream.  (In the Matrix Market
///   format, sparse matrix entries cannot start until at least line 3
///   of the file.)  This is only valid on Process 0.
/// \param totalNumEntRead [out] Total number of matrix entries
///   (triples) read on Process 0.  This is only valid on Process 0.
/// \param processTriple [in] Closure, generally with side effects,
///   that takes in and stores off a sparse matrix entry.  First
///   argument is the (global) row index, second argument is the
///   (global) column index, and third argument is the value of the
///   entry.  The closure must NOT do MPI communication.  Return value
///   is an error code, that is zero if and only if the closure
///   succeeded.  We intend for you to use this to call
///   CooMatrix::insertEntry.
/// \param comm [in] Communicator to use for receiving the triples.
/// \param tolerant [in] Whether to read tolerantly.
/// \param errStrm [in] If not NULL, print any error messages to this
///   stream.
///
/// \return Error code; 0 if and only if success.
template<class SC, class GO>
int
readAndDealOutTriples (std::istream& inputStream, // only valid on Proc 0
                       std::size_t& curLineNum, // only valid on Proc 0
                       std::size_t& totalNumEntRead, // only valid on Proc 0
                       std::function<int (const GO, const GO, const SC&)> processTriple,
                       const std::size_t maxNumEntPerMsg,
                       const ::Teuchos::Comm<int>& comm,
                       const bool tolerant = false,
                       std::ostream* errStrm = NULL,
                       const bool debug = false)
{
  using Kokkos::ArithTraits;
  using std::endl;
  using std::size_t;

  constexpr int srcRank = 0;
  //constexpr int sizeTag = 42 + (ArithTraits<SC>::is_complex ? 100 : 0);
  ////constexpr int msgTag = 43 + (ArithTraits<SC>::is_complex ? 100 : 0);
  constexpr int sizeTag = 42;
  //constexpr int msgTag = 43;
  const int myRank = comm.getRank ();
  const int numProcs = comm.getSize ();
  int errCode = 0;

  ::Teuchos::ArrayRCP<int> sizeBuf (1);
  ::Teuchos::ArrayRCP<char> msgBuf; // to be resized as needed

  // Temporary storage for reading & packing (on Process srcRank) or
  // unpacking (every other process) triples.
  std::vector<GO> rowInds;
  std::vector<GO> colInds;
  std::vector<SC> vals;
  rowInds.reserve (maxNumEntPerMsg);
  colInds.reserve (maxNumEntPerMsg);
  vals.reserve (maxNumEntPerMsg);

  totalNumEntRead = 0;
  if (myRank == srcRank) {
    // Loop around through all the processes, including this one, over
    // and over until we reach the end of the file, or an error occurs.
    int destRank = 0;
    bool lastMessageWasLegitZero = false;
    for ( ;
          ! inputStream.eof () && errCode == 0;
         destRank = (destRank + 1) % numProcs) {

      size_t curNumEntRead = 0; // output argument of below
      if (destRank == srcRank) {
        // We can read and process the triples directly.  We don't
        // need to use intermediate storage, because we don't need to
        // pack and send the triples.
        const int readErrCode =
          Impl::readTriples<SC, GO> (inputStream, curLineNum, curNumEntRead,
                                     processTriple, maxNumEntPerMsg, tolerant,
                                     errStrm, debug);
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "(dest=src) readTriples returned curNumEntRead="
             << curNumEntRead << ", errCode=" << readErrCode << endl;
          *errStrm << os.str ();
        }
        errCode = (readErrCode != 0) ? readErrCode : errCode;
      }
      else {
        if (false && debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Calling readAndSend... with destRank=" << destRank << endl;
          *errStrm << os.str ();
        }
        // Read, pack, and send the triples to destRank.
        const int readAndSendErrCode =
          Impl::readAndSendOneBatchOfTriples<SC, GO> (inputStream, curLineNum,
                                                      curNumEntRead,
                                                      sizeBuf, msgBuf,
                                                      rowInds, colInds, vals,
                                                      maxNumEntPerMsg, destRank,
                                                      comm, tolerant, errStrm,
                                                      debug);
        totalNumEntRead += curNumEntRead;
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "readAndSend... with destRank=" << destRank
             << " returned curNumEntRead=" << curNumEntRead
             << ", errCode=" << readAndSendErrCode << endl;
          *errStrm << os.str ();
        }
        errCode = (readAndSendErrCode != 0) ? readAndSendErrCode : errCode;
        if (readAndSendErrCode == 0 && curNumEntRead == 0) {
          lastMessageWasLegitZero = true;
          if (debug && errStrm != NULL) {
            std::ostringstream os;
            os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
               << ", GO=" << typeid (GO).name () << ": "
               << "Last send to " << destRank << " with tag " << sizeTag
               << " was legit zero, counts as termination" << endl;
            *errStrm << os.str ();
          }
        }
      }
    } // loop around through processes until done reading file, or error

    // Loop around through the remaining processes, and tell them that
    // we're done, by sending zero.  If the last message we sent to
    // destRank was zero, then skip that process, since it only
    // expects one message of size zero.  Note that destRank got
    // incremented mod numProcs at end of loop, so we have to
    // decrement it mod numProcs.
    destRank = (destRank - 1) % numProcs;
    if (destRank < 0) { // C mod operator does not promise positivity
      destRank = destRank + numProcs;
    }

    const int startRank = lastMessageWasLegitZero ? (destRank+1) : destRank;
    for (int outRank = startRank; outRank < numProcs; ++outRank) {
      if (outRank != srcRank) {
        if (debug && errStrm != NULL) {
          std::ostringstream os;
          os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
             << ", GO=" << typeid (GO).name () << ": "
             << "Post send (size, termination msg) to " << outRank
             << " with tag " << sizeTag << "(was last message legit zero? "
             << (lastMessageWasLegitZero ? "true" : "false") << ")" << endl;
          *errStrm << os.str ();
        }
        sizeBuf[0] = 0;
        ::Teuchos::send (sizeBuf.getRawPtr (), 1, outRank, sizeTag, comm);
      }
    }
  }
  else {
    while (true) {
      // Prepost a message to receive the size (in bytes) of the
      // incoming packet.
      sizeBuf[0] = 0; // superfluous, but safe
      if (debug && errStrm != NULL) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
           << ", GO=" << typeid (GO).name () << ": "
           << "Post irecv (size) from " << srcRank
           << " with tag " << sizeTag << std::endl;
        *errStrm << os.str ();
      }
      auto sizeReq = ::Teuchos::ireceive (sizeBuf, srcRank, sizeTag, comm);

      int numEnt = 0; // output argument
      const int recvErrCode =
        Impl::recvOneBatchOfTriples (rowInds, colInds, vals, numEnt, sizeBuf,
                                     msgBuf, sizeReq, srcRank, comm, tolerant,
                                     errStrm, debug);
      if (debug && errStrm != NULL) {
        std::ostringstream os;
        os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
           << ", GO=" << typeid (GO).name () << ": "
           << "recvOneBatchOfTriples returned numEnt=" << numEnt
           << ", errCode=" << recvErrCode << endl;
        *errStrm << os.str ();
      }
      errCode = (recvErrCode != 0) ? recvErrCode : errCode;

      if (numEnt != static_cast<int> (rowInds.size ()) ||
          numEnt != static_cast<int> (colInds.size ()) ||
          numEnt != static_cast<int> (vals.size ())) {
        errCode = (errCode == 0) ? -1 : errCode;
        if (errStrm != NULL) {
          *errStrm << "recvOneBatchOfTriples produced inconsistent data sizes.  "
                   << "numEnt = " << numEnt
                   << ", rowInds.size() = " << rowInds.size ()
                   << ", colInds.size() = " << colInds.size ()
                   << ", vals.size() = " << vals.size () << "."
                   << endl;
        }
      } // if sizes inconsistent

      // Sending zero items is how Process srcRank tells this process
      // that it (Process srcRank) is done sending out data.
      if (numEnt == 0) {
        break;
      }

      for (int k = 0; k < numEnt && errCode == 0; ++k) {
        const int curErrCode = processTriple (rowInds[k], colInds[k], vals[k]);
        errCode = (curErrCode == 0) ? errCode : curErrCode;
      }
    } // while we still get messages from srcRank
  }

  if (debug && errStrm != NULL) {
    std::ostringstream os;
    os << "Proc " << comm.getRank () << ", SC=" << typeid (SC).name ()
       << ", GO=" << typeid (GO).name () << ": "
       << "Done with send/recv loop" << endl;
    *errStrm << os.str ();
  }
  // Do a bitwise OR to get an error code that is nonzero if and only
  // if any process' local error code is nonzero.
  using ::Teuchos::outArg;
  using ::Teuchos::REDUCE_BOR;
  using ::Teuchos::reduceAll;
  const int lclErrCode = errCode;
  reduceAll<int, int> (comm, REDUCE_BOR, lclErrCode, outArg (errCode));
  return errCode;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_READTRIPLES_HPP
