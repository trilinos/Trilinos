// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"

namespace Teuchos {

  namespace {

    // Pack the given array of strings into a single string with an
    // offsets array.  This is a helper routine of \c sendStrings().
    // For strings[k], offsets[k] gives its start position in
    // packedString, and offsets[k+1]-ofsets[k] gives its length.
    // Thus, packedString.substr (offsets[k], offsets[k+1]-offsets[k])
    // == strings[k].
    //
    // \param packedString [out] The packed string.  It will be
    //   resized on output as necessary.
    //
    // \param offsets [out] Array of offsets, of length one more than
    //   the number of strings in the \c strings array.  Thus, the
    //   offsets array always has positive length.  
    //
    // \param strings [in] The array of strings to pack.  It may have
    //   zero length.  In that case, on output, the offsets array will
    //   have length one, offsets[0] = 0, and packedString will have
    //   length zero.
    //
    // Although std::string could be considered an array, it does not
    // have a size_type typedef.  Instead, it uses size_t for that
    // purpose.
    void
    packStringsForSend (std::string& packedString, 
			Array<size_t>& offsets,
			const Array<std::string>& strings)
    {
      using std::cerr;
      using std::endl;
      using std::string;

      const bool debug = false;

      // Compute index offsets in the packed string.
      offsets.resize (strings.size() + 1);
      size_t totalLength = 0;
      Array<size_t>::size_type offsetsIndex = 0;
      for (Array<string>::const_iterator it = strings.begin();
	   it != strings.end(); ++it, ++offsetsIndex)
	{
	  offsets[offsetsIndex] = totalLength;
	  totalLength += it->size();
	}
      offsets[offsetsIndex] = totalLength;

      // Pack the array of strings into the packed string.
      packedString.resize (totalLength);
      string::iterator packedStringIter = packedString.begin();
      for (Array<string>::const_iterator it = strings.begin(); 
	   it != strings.end(); ++it)
	packedStringIter = std::copy (it->begin(), it->end(), packedStringIter);

      if (debug)
	{
	  std::ostringstream out;
	  RCP<const Comm<int> > pComm = DefaultComm<int>::getComm ();
	  out << "Proc " << pComm->getRank() << ": in pack: offsets = [";
	  for (Array<size_t>::const_iterator it = offsets.begin(); 
	       it != offsets.end(); ++it)
	    {
	      out << *it;
	      if (it + 1 != offsets.end())
		out << ", ";
	    }
	  out << "], packedString = " << packedString << endl;
	  cerr << out.str();
	}
    }

    // \brief Send an array of strings.
    //
    // Teuchos::send() (or rather, Teuchos::SerializationTraits)
    // doesn't know how to send an array of strings.  This function
    // packs an array of strings into a single string with an offsets
    // array, and sends the offsets array (and the packed string, if
    // it is not empty).
    void
    sendStrings (const Comm<int>& comm, // in
		 const Array<std::string>& strings, // in
		 const int destRank) // in
    {
      // Pack the string array into the packed string, and compute
      // offsets.
      std::string packedString;
      Array<size_t> offsets;
      packStringsForSend (packedString, offsets, strings);
      TEST_FOR_EXCEPTION(offsets.size() == 0, std::logic_error, 
			 "packStringsForSend() returned a zero-length offsets "
			 "array on MPI Proc " << comm.getRank() << ", to be "
			 "sent to Proc " << destRank << ".  The offsets array "
			 "should always have positive length.  Please report "
			 "this bug to the Teuchos developers.");

      // Send the count of offsets.
      send (comm, offsets.size(), destRank);

      // Send the array of offsets.  There is always at least one
      // element in the offsets array, so we can always take the
      // address of the first element.
      const int offsetsSendCount = static_cast<int> (offsets.size());
      send (comm, offsetsSendCount, &offsets[0], destRank);

      // Now send the packed string.  It may be empty if the strings
      // array has zero elements or if all the strings in the array
      // are empty.  If the packed string is empty, we don't send
      // anything, since the receiving process already knows (from the
      // offsets array) not to expect anything.
      const int stringSendCount = static_cast<int> (packedString.size());
      if (stringSendCount > 0)
	send (comm, stringSendCount, &packedString[0], destRank);
    }

    void
    unpackStringsAfterReceive (Array<std::string>& strings,
			       const std::string& packedString,
			       const Array<size_t> offsets)
    {
      const bool debug = false;
      if (debug)
	{
	  using std::cerr;
	  using std::endl;

	  std::ostringstream out;
	  RCP<const Comm<int> > pComm = DefaultComm<int>::getComm ();
	  out << "Proc " << pComm->getRank() << ": in unpack: offsets = [";
	  for (Array<size_t>::const_iterator it = offsets.begin(); 
	       it != offsets.end(); ++it)
	    {
	      out << *it;
	      if (it + 1 != offsets.end())
		out << ", ";
	    }
	  out << "], packedString = " << packedString << endl;
	  cerr << out.str();
	}
      TEST_FOR_EXCEPTION(offsets.size() == 0, std::logic_error, 
			 "The offsets array has length zero, which does not "
			 "make sense.  Even when sending / receiving zero "
			 "strings, the offsets array should have one entry "
			 "(namely, zero).");
      const Array<size_t>::size_type numStrings = offsets.size() - 1;
      strings.resize (numStrings);
      for (Array<size_t>::size_type k = 0; k < numStrings; ++k)
	{ // Exclusive index range in the packed string in which to
	  // find the current string.
	  const size_t start = offsets[k];
	  const size_t end = offsets[k+1];
	  strings[k] = packedString.substr (start, end - start);
	}
    }

    // Function corresponding to \c sendStrings() that receives an
    // array of strings (Array<std::string>) in packed form.
    void
    receiveStrings (const Comm<int>& comm,
		    const int sourceRank,
		    Array<std::string>& strings)
    {
      // Receive the number of offsets.  There should always be at
      // least 1 offset.
      Array<size_t>::size_type numOffsets = 0;
      receive (comm, sourceRank, &numOffsets);
      TEST_FOR_EXCEPTION(numOffsets == 0, std::logic_error, 
			 "Invalid number of offsets numOffsets=" << numOffsets 
			 << " received on MPI Rank " << comm.getRank() 
			 << " from Rank " << sourceRank << ".  Please report "
			 "this bug to the Teuchos developers.");

      // Receive the array of offsets.
      Array<size_t> offsets (numOffsets);
      const int offsetsRecvCount = static_cast<int> (numOffsets);
      receive (comm, sourceRank, offsetsRecvCount, &offsets[0]);

      // If the packed string is nonempty, receive the packed string,
      // and unpack it.  The last entry of offsets is the length of
      // the packed string.
      std::string packedString (offsets.back(), ' ');
      const int stringRecvCount = static_cast<int> (offsets.back());
      if (stringRecvCount > 0)
	{
	  receive (comm, sourceRank, stringRecvCount, &packedString[0]);
	  unpackStringsAfterReceive (strings, packedString, offsets);
	}
    }

    // \brief Helper function for \c mergeCounterNamesHelper().
    //
    // The \c mergeCounterNamesHelper() function implements (using a
    // parallel reduction) the set union resp. intersection (depending
    // on the \c setOp argument) of the MPI process' sets of counter
    // names.  This function implements the binary associative
    // operator which computes the set union resp. intersection of two
    // sets: the "left" process' intermediate reduction result (global
    // counter names), and the "mid" process' local counter names.
    // 
    // \param comm [in] Communicator for which \c
    //   mergeCounterNamesHelper() was called.
    //
    // \param myRank [in] Rank of the calling MPI process; must be
    //   either == left or == mid.
    //
    // \param left [in] The "left" input argument of
    //   \c mergeCounterNamesHelper().
    //
    // \param mid [in] The value of "mid" in the implementation of \c
    //   mergeCounterNamesHelper().
    //
    // \param localNames [in] List of counter names belonging to the
    //   calling MPI process.
    //
    // \param globalNames [in/out] Only accessed if myRank == left.
    //   If so, on input: the intermediate reduction result of the
    //   union resp. intersection (depending on \c setOp).  On output:
    //   the union resp. intersection of the input value of the "left"
    //   MPI process' globalNames with the "mid" MPI process'
    //   localNames.
    //
    // \param setOp [in] If Intersection, compute the set intersection
    //   of counter names, else if Union, compute the set union of
    //   counter names.
    void
    mergeCounterNamesPair (const Comm<int>& comm, 
			   const int myRank,
			   const int left,
			   const int mid,
			   const Array<std::string>& localNames,
			   Array<std::string>& globalNames,
			   const ECounterSetOp setOp)
    {
      using std::cerr;
      using std::endl;
      using std::string;

      const bool debug = false;

      if (myRank == left)
	{ // Receive names from the other process, and merge its names
	  // with the names on this process.
	  Array<string> otherNames;
	  receiveStrings (comm, mid, otherNames);

	  if (debug)
	    {
	      // Buffering locally in an ostringstream before writing to
	      // the shared stderr sometimes helps avoid interleaved
	      // debugging output.
	      std::ostringstream out;
	      out << "Proc " << myRank << ": in mergePair: otherNames = [";
	      for (Array<std::string>::const_iterator it = otherNames.begin(); 
		   it != otherNames.end(); ++it)
		{
		  out << "\"" << *it << "\"";
		  if (it + 1 != otherNames.end())
		    out << ", ";
		}
	      out << "]" << endl;
	      cerr << out.str();
	    }

	  // Assume that both globalNames and otherNames are sorted.
	  // Compute the set intersection / union as specified by the
	  // enum.
	  Array<string> newNames;
	  if (setOp == Intersection)
	    std::set_intersection (globalNames.begin(), globalNames.end(),
				   otherNames.begin(), otherNames.end(),
				   std::back_inserter (newNames));
	  else if (setOp == Union)
	    std::set_union (globalNames.begin(), globalNames.end(),
			    otherNames.begin(), otherNames.end(),
			    std::back_inserter (newNames));
	  else
	    TEST_FOR_EXCEPTION(setOp != Intersection && setOp != Union,
			       std::logic_error,
			       "Invalid set operation enum value.  Please "
			       "report this bug to the Teuchos developers.");
	  globalNames.swap (newNames);
	}
      else if (myRank == mid)
	sendStrings (comm, localNames, left);
      else
	TEST_FOR_EXCEPTION(myRank != left && myRank != mid, 
			   std::logic_error,
			   "myRank=" << myRank << " is neither left=" << left
			   << " nor mid=" << mid << ".  Please report this "
			   "bug to the Teuchos developers.");
    }

    // Recursive helper function for \c mergeCounterNames().
    //
    // This function implements the set union resp. intersection
    // (depending on the \c setOp argument) of the MPI process' sets
    // of counter names, using a parallel reduction. (Since the
    // Teuchos comm wrappers as of 11 July 2011 lack a wrapper for
    // MPI_Reduce(), we hand-roll the reduction using a binary tree
    // via recursion.  We don't need an all-reduce in this case.)
    void
    mergeCounterNamesHelper (const Comm<int>& comm, 
			 const int myRank,
			 const int left,
			 const int right, // inclusive range [left, right]
			 const Array<std::string>& localNames,
			 Array<std::string>& globalNames,
			 const ECounterSetOp setOp)
    {
      // Correctness proof:
      //
      // 1. Both set intersection and set union are associative (and
      //    indeed even commutative) operations.
      // 2. mergeCounterNamesHelper() is just a reduction by binary tree.
      // 3. Reductions may use any tree shape as long as the binary
      //    operation is associative.
      //
      // Recursive "reduction" algorithm:
      //
      // Let mid be the midpoint of the inclusive interval [left,
      // right].  If the (intersection, union) of [left, mid-1] and
      // the (intersection, union) of [mid, right] are both computed
      // correctly, then the (intersection, union) of these two sets
      // is the (intersection, union) of [left, right].
      //
      // The base case is left == right: the (intersection, union) of
      // one set is simply that set, so copy localNames into
      // globalNames.
      //
      // We include another base case for safety: left > right,
      // meaning that the set of processes is empty, so we do nothing
      // (the (intersection, union) of an empty set of sets is the
      // empty set).
      if (left > right)
	return;
      else if (left == right)
	{
	  Array<string> newNames;
	  newNames.reserve (localNames.size());
	  std::copy (localNames.begin(), localNames.end(), 
		     std::back_inserter (newNames));
	  globalNames.swap (newNames);
	}
      else
	{ // You're sending messages across the network, so don't
	  // bother to optimize away a few branches here.
	  //
	  // Recurse on [left, mid-1] or [mid, right], depending on myRank.
	  const int mid = left + (right - left + 1) / 2;
	  if (myRank >= left && myRank <= mid-1)
	    mergeCounterNamesHelper (comm, myRank, left, mid-1, 
				     localNames, globalNames, setOp);
	  else if (myRank >= mid && myRank <= right)
	    mergeCounterNamesHelper (comm, myRank, mid, right,
				     localNames, globalNames, setOp);
	  else
	    return; // Don't call mergeCounterNamesPair() if not participating.

	  // Combine the results of the recursive step.
	  if (myRank == left || myRank == mid)
	    mergeCounterNamesPair (comm, myRank, left, mid, 
				   localNames, globalNames, setOp);
	}
    }

  } // namespace (anonymous)

  void
  mergeCounterNames (const Comm<int>& comm, 
		     const Array<std::string>& localNames,
		     Array<std::string>& globalNames,
		     const ECounterSetOp setOp)
  {
    const int myRank = comm.getRank();
    const int left = 0;
    const int right = comm.getSize() - 1;
    Array<std::string> theGlobalNames;
    mergeCounterNamesHelper (comm, myRank, left, right, 
			     localNames, theGlobalNames, setOp);
    // "Transactional" semantics ensure strong exception safety for
    // output.
    globalNames.swap (theGlobalNames);
  }

} // namespace Teuchos
