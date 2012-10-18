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

#include "Teuchos_CommHelpers.hpp"

namespace Teuchos {
namespace { // (anonymous)

#ifdef HAVE_MPI
MPI_Op getMpiOpForEReductionType (const enum EReductionType reductionType) {
  switch (reductionType) {
  case REDUCE_SUM: return MPI_SUM;
  case REDUCE_MIN: return MPI_MIN;
  case REDUCE_MAX: return MPI_MAX;
  case REDUCE_AND: return MPI_LAND; // logical AND, not bitwise AND
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The given EReductionType value is invalid."); 
  }
}
  
std::string getMpiErrorString (const int errCode) {
  // Space for storing the error string returned by MPI.
  // Leave room for null termination, since I don't know if MPI does this.
  char errString [MPI_MAX_ERROR_STRING+1];
  int errStringLen = MPI_MAX_ERROR_STRING; // output argument
  (void) MPI_Error_string (errCode, errString, &errStringLen);
  // errStringLen on output is the number of characters written.
  // I'm not sure (the MPI 3.0 Standard doesn't say) if this
  // includes the '\0', so I'll make sure.  We reserved space for
  // the extra '\0' if needed.
  if (errString[errStringLen-1] != '\0') {
    errString[errStringLen] = '\0';
  }
  return std::string (errString); // This copies the original string.
}
#endif // HAVE_MPI

} // namespace (anonymous)

// mfh 18 Oct 2012: Note on full template specializations
//
// To make Windows builds happy, declarations of full template
// specializations (as found in Teuchos_CommHelpers.hpp) must use the
// TEUCHOS_LIB_DLL_EXPORT macro.  However, _definitions_ of the
// specializations (as found in this file) must _not_ use the macro.
// That's why we don't use that macro here.

#ifdef TEUCHOS_HAVE_COMPLEX
// Specialization for Ordinal=int and Packet=std::complex<double>.
template<>
void 
reduceAll<int, std::complex<double> > (const Comm<int>& comm, 
				       const EReductionType reductType,
				       const int count, 
				       const std::complex<double> sendBuffer[], 
				       std::complex<double> globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, std::complex<double> > (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int, std::complex<double> > > reductOp (createOp<int, std::complex<double> > (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<std::complex<double>* > (sendBuffer), 
      globalReducts, count, MPI_C_DOUBLE_COMPLEX, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}

// Specialization for Ordinal=int and Packet=std::complex<float>.
template<>
void 
reduceAll<int, std::complex<float> > (const Comm<int>& comm, 
				       const EReductionType reductType,
				       const int count, 
				       const std::complex<float> sendBuffer[], 
				       std::complex<float> globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, std::complex<float> > (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int, std::complex<float> > > reductOp (createOp<int, std::complex<float> > (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<std::complex<float>* > (sendBuffer), 
      globalReducts, count, MPI_C_FLOAT_COMPLEX, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}
#endif // TEUCHOS_HAVE_COMPLEX


// Specialization for Ordinal=int and Packet=double.
template<>
void 
reduceAll<int, double> (const Comm<int>& comm, 
			const EReductionType reductType,
			const int count, 
			const double sendBuffer[], 
			double globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, double> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,double> > reductOp (createOp<int,double> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<double*> (sendBuffer), globalReducts, count, MPI_DOUBLE, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


// Specialization for Ordinal=int and Packet=float.
template<>
void 
reduceAll<int, float> (const Comm<int>& comm, 
		       const EReductionType reductType,
		       const int count, 
		       const float sendBuffer[], 
		       float globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, float> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,float> > reductOp (createOp<int,float> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<float*> (sendBuffer), globalReducts, count, MPI_FLOAT, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


#ifdef TEUCHOS_HAVE_LONG_LONG_INT
// Specialization for Ordinal=int and Packet=long long.
template<>
void 
reduceAll<int, long long> (const Comm<int>& comm, 
			   const EReductionType reductType,
			   const int count, 
			   const long long sendBuffer[], 
			   long long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, long long> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,long long> > reductOp (createOp<int,long long> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // SerialComm means one process only.  
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // The Comm is an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<long long*> (sendBuffer), globalReducts, count, MPI_LONG_LONG, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}
#endif // TEUCHOS_HAVE_LONG_LONG_INT


// Specialization for Ordinal=int and Packet=long.
template<>
void 
reduceAll<int, long> (const Comm<int>& comm, 
		      const EReductionType reductType,
		      const int count, 
		      const long sendBuffer[], 
		      long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, long> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,long> > reductOp (createOp<int,long> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // SerialComm means one process only.  
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // The Comm is an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<long*> (sendBuffer), globalReducts, count, MPI_LONG, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


// Specialization for Ordinal=int and Packet=int.
template<>
void 
reduceAll<int, int> (const Comm<int>& comm, 
		     const EReductionType reductType,
		     const int count, 
		     const int sendBuffer[], 
		     int globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, int> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,int> > reductOp (createOp<int,int> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // SerialComm means one process only.  
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // The Comm is an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<int*> (sendBuffer), globalReducts, count, MPI_INT, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


// Specialization for Ordinal=int and Packet=short.
template<>
void 
reduceAll<int, short> (const Comm<int>& comm, 
		       const EReductionType reductType,
		       const int count, 
		       const short sendBuffer[], 
		       short globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, short> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,short> > reductOp (createOp<int,short> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // SerialComm means one process only.  
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // The Comm is an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<short*> (sendBuffer), globalReducts, count, MPI_SHORT, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


// Specialization for Ordinal=int and Packet=char.
template<>
void 
reduceAll<int, char> (const Comm<int>& comm, 
		      const EReductionType reductType,
		      const int count, 
		      const char sendBuffer[], 
		      char globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, char> (" << count << ", " << toString (reductType) << ")"
    );
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) { 
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      std::auto_ptr<ValueTypeReductionOp<int,char> > reductOp (createOp<int,char> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    } 
    else { // SerialComm means one process only.  
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // The Comm is an MpiComm.  Invoke MPI directly.
    MPI_Op op = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    const int err = MPI_Allreduce (const_cast<char*> (sendBuffer), globalReducts, count, MPI_CHAR, op, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "MPI_Allreduce failed with the following error: " 
      << getMpiErrorString (err));
  }
#else 
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


} // namespace Teuchos
