// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Details_MpiTypeTraits.hpp"
#include "mpi.h"

// mfh 10 Nov 2016: This file depends on MPI >= 2 functions and
// predefined constants.  Teuchos::Details::MpiTypeTraits should still
// work with MPI 1, so we simply disable the test in that case.  We
// don't normally test with very old MPI implementations, so I would
// rather not claim to support MPI 1.
#if MPI_VERSION >= 2

namespace { // (anonymous)

#define REPORT_MPI_ERR( errCode, mpiFunctionName ) \
do { \
  if (errCode != MPI_SUCCESS) { \
    char errString[MPI_MAX_ERROR_STRING]; \
    int len = 0; \
    (void) MPI_Error_string (errCode, errString, &len); \
    TEUCHOS_TEST_FOR_EXCEPTION \
      (true, std::logic_error, "MPI routine " << mpiFunctionName << " returned " \
       "err != MPI_SUCCESS.  Reported error: \"" << errString << "\"."); \
  } \
} while (false)


//! Is the given MPI_Datatype a custom (derived, not built-in) datatype?
bool
mpiDatatypeIsCustom (MPI_Datatype dt)
{
#if MPI_VERSION < 2
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "mpiDatatypeIsCustom: This function requires MPI "
     "2.0 at least, since it relies on MPI_Type_get_envelope.  MPI 2.0 came "
     "out in 1996, so I think it's time to upgrade your MPI implementation.");
#else // MPI_VERSION >= 2
  int numInts = 0;
  int numAddrs = 0;
  int numTypes = 0;
  int combiner = 0;
  // This does not exist in MPI 1.1, but does exist in MPI 2.0.
  const int err =
    MPI_Type_get_envelope (dt, &numInts, &numAddrs, &numTypes, &combiner);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_get_envelope "
     "returned err != MPI_SUCCESS.");

  switch (combiner) {
  case MPI_COMBINER_NAMED:
    return false;
  case MPI_COMBINER_DUP:
  case MPI_COMBINER_CONTIGUOUS:
  case MPI_COMBINER_VECTOR:
  case MPI_COMBINER_HVECTOR:
  case MPI_COMBINER_INDEXED:
  case MPI_COMBINER_HINDEXED:
  case MPI_COMBINER_INDEXED_BLOCK:
#if MPI_VERSION >= 3
  case MPI_COMBINER_HINDEXED_BLOCK: // in MPI 3.0, not in MPI 2.{0,1,2}
#endif // 0
  case MPI_COMBINER_STRUCT:
  case MPI_COMBINER_SUBARRAY:
  case MPI_COMBINER_DARRAY:
  case MPI_COMBINER_F90_REAL:
  case MPI_COMBINER_F90_COMPLEX:
  case MPI_COMBINER_F90_INTEGER:
  case MPI_COMBINER_RESIZED:
  default:
    return true;
  }
#endif // MPI_VERSION >= 2
}

// mfh 09 Nov 2016: MPI (at least as of 3.1) has no built-in function
// that lets you compare MPI_Datatype instances.  Direct equality
// comparison (==) certainly does not work, because it's possible to
// create two distinct MPI_Datatype instances that represent the same
// type.  The right way to do this is to decode the datatype, using
// MPI_Type_get_envelope (to get sizes of arrays) and
// MPI_Type_get_contents recursively, until one reaches primitive
// (built-in, not custom) MPI_Datatype instances.  This recursion
// takes stack space and could, in theory, overflow.  However, most
// MPI_Datatype instances are not deeply nested.

bool
mpiDatatypeIsSame (MPI_Datatype t1, MPI_Datatype t2)
{
#if MPI_VERSION < 2
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "mpiDatatypeIsSame: This function requires MPI "
     "2.0 at least, since it relies on MPI_Type_get_envelope and "
     "MPI_Type_get_contents.  MPI 2.0 came out in 1996, so I think it's time "
     "to upgrade your MPI implementation.");
#else // MPI_VERSION >= 2
  int err = MPI_SUCCESS;

  int numInts1, numAddrs1, numTypes1, combiner1;
  // This does not exist in MPI 1.1, but does exist in MPI 2.0.
  err = MPI_Type_get_envelope (t1, &numInts1, &numAddrs1, &numTypes1, &combiner1);
  REPORT_MPI_ERR(err, "MPI_Type_get_envelope");
  int numInts2, numAddrs2, numTypes2, combiner2;
  // This does not exist in MPI 1.1, but does exist in MPI 2.0.
  err = MPI_Type_get_envelope (t2, &numInts2, &numAddrs2, &numTypes2, &combiner2);
  REPORT_MPI_ERR(err, "MPI_Type_get_envelope");

  if (combiner1 != combiner2 ||
      numInts1 != numInts2 ||
      numAddrs1 != numAddrs2 ||
      numTypes1 != numTypes2) {
    return false;
  }
  else if (combiner1 == MPI_COMBINER_NAMED) {
    // For built-in MPI_Datatype, OpenMPI 1.10.1 reports an internal
    // error when one attempts to call MPI_Type_get_contents.  On
    // the other hand, for built-in MPI_Datatype, one may compare
    // the "combiner" values directly.
    return t1 == t2;
  }
  else if (numTypes1 == 0) {
    return true; // no types to check
  }
  else {
    // The most general recursive case takes some stack space, and
    // thus could, in theory, overflow.  However, most MPI_Datatype
    // instances in practice are not deeply nested.
    std::vector<MPI_Datatype> theTypes1 (numTypes1);
    std::vector<MPI_Datatype> theTypes2 (numTypes2);

    // We don't actually need anything but the arrays of addresses
    // for the recursion.  Restricting the other std::vector
    // instances to this inner scope that closes before the
    // recursion ensures that they don't stick around on the stack.
    {
      // Minimum of one entry, to please MPI_Type_get_contents.
      std::vector<int> theInts1 (numInts1);
      std::vector<int> theInts2 (numInts2);
      std::vector<MPI_Aint> theAddrs1 (numAddrs1);
      std::vector<MPI_Aint> theAddrs2 (numAddrs2);

      // This does not exist in MPI 1.1, but does exist in MPI 2.0.
      err = MPI_Type_get_contents (t1, numInts1, numAddrs1, numTypes1,
                                   numInts1 == 0 ? NULL : &theInts1[0],
                                   numAddrs1 == 0 ? NULL : &theAddrs1[0],
                                   numTypes1 == 0 ? NULL : &theTypes1[0]);
      REPORT_MPI_ERR(err, "MPI_Type_get_contents");
      // This does not exist in MPI 1.1, but does exist in MPI 2.0.
      err = MPI_Type_get_contents (t2, numInts2, numAddrs2, numTypes2,
                                   numInts2 == 0 ? NULL : &theInts2[0],
                                   numAddrs2 == 0 ? NULL : &theAddrs2[0],
                                   numTypes2 == 0 ? NULL : &theTypes2[0]);
      REPORT_MPI_ERR(err, "MPI_Type_get_contents");
      TEUCHOS_TEST_FOR_EXCEPTION
        (err != MPI_SUCCESS, std::logic_error, "MPI_Type_get_contents "
         "returned err != MPI_SUCCESS.");

      for (int k = 0; k < numInts1; ++k) {
        if (theInts1[k] != theInts2[k]) {
          return false;
        }
      }
      for (int k = 0; k < numAddrs1; ++k) {
        if (theAddrs1[k] != theAddrs2[k]) {
          return false;
        }
      }
    }

    // Compare all the MPI_Datatype instances, recursively.
    for (int k = 0; k < numTypes1; ++k) {
      const bool same = mpiDatatypeIsSame (theTypes1[k], theTypes2[k]);
      // For non-built-in (custom) MPI_Datatype instances, the
      // instance returned from MPI_Type_get_contents is a "new"
      // instance and must therefore be freed.  It is illegal to
      // call MPI_Type_free on a built-in MPI_Datatype instance.
      if (mpiDatatypeIsCustom (theTypes1[k])) {
        err = MPI_Type_free (&theTypes1[k]);
        TEUCHOS_TEST_FOR_EXCEPTION
          (err != MPI_SUCCESS, std::logic_error, "MPI_Type_free "
           "returned err != MPI_SUCCESS.");
      }
      if (mpiDatatypeIsCustom (theTypes2[k])) {
        err = MPI_Type_free (&theTypes2[k]);
        TEUCHOS_TEST_FOR_EXCEPTION
          (err != MPI_SUCCESS, std::logic_error, "MPI_Type_free "
           "returned err != MPI_SUCCESS.");
      }
      if (! same) {
        return false;
      }
    }
    return true;
  }
#endif // MPI_VERSION < 2
}

TEUCHOS_UNIT_TEST( MpiTypeTraits, TestMpiDatatypeIsCustom )
{
  using Teuchos::Details::MpiTypeTraits;
  using std::endl;
  int err = MPI_SUCCESS;

  out << "Test mpiDatatypeIsCustom" << endl;
  Teuchos::OSTab tab1 (out);

  // In order to debug any MPI errors, tell MPI to let its functions
  // return error codes, instead of aborting right away.
  //
  // This function does not exist in MPI 1.1, but does exist in MPI
  // 2.0.  MPI 1 has MPI_Errhandler_set, which is deprecated as of MPI
  // 2.0.
  err = MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Comm_set_errhandler "
     "returned err != MPI_SUCCESS.");

  TEST_ASSERT( ! mpiDatatypeIsCustom (MPI_INT) );
  TEST_ASSERT( ! mpiDatatypeIsCustom (MPI_DOUBLE) );

  MPI_Datatype mpi3doubles;
  err = MPI_Type_contiguous (3, MPI_DOUBLE, &mpi3doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_contiguous "
     "returned err != MPI_SUCCESS.");

  TEST_ASSERT( mpiDatatypeIsCustom (mpi3doubles) );

  err = MPI_Type_free (&mpi3doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_free "
     "returned err != MPI_SUCCESS.");
}

TEUCHOS_UNIT_TEST( MpiTypeTraits, TestMpiDatatypeIsSame )
{
  using Teuchos::Details::MpiTypeTraits;
  using std::endl;
  int err = MPI_SUCCESS;

  out << "Test mpiDatatypeIsSame" << endl;
  Teuchos::OSTab tab1 (out);

  TEST_ASSERT( mpiDatatypeIsSame (MPI_INT, MPI_INT) );
  TEST_ASSERT( mpiDatatypeIsSame (MPI_DOUBLE, MPI_DOUBLE) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_DOUBLE, MPI_INT) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_INT, MPI_DOUBLE) );

  MPI_Datatype mpi3doubles;
  err = MPI_Type_contiguous (3, MPI_DOUBLE, &mpi3doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_contiguous "
     "returned err != MPI_SUCCESS.");
  MPI_Datatype mpi4doubles;
  err = MPI_Type_contiguous (4, MPI_DOUBLE, &mpi4doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_contiguous "
     "returned err != MPI_SUCCESS.");

  TEST_ASSERT( ! mpiDatatypeIsSame (mpi3doubles, MPI_DOUBLE) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_DOUBLE, mpi3doubles) );
  TEST_ASSERT( ! mpiDatatypeIsSame (mpi3doubles, MPI_INT) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_INT, mpi3doubles) );
  TEST_ASSERT( ! mpiDatatypeIsSame (mpi4doubles, MPI_DOUBLE) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_DOUBLE, mpi4doubles) );
  TEST_ASSERT( ! mpiDatatypeIsSame (mpi4doubles, MPI_INT) );
  TEST_ASSERT( ! mpiDatatypeIsSame (MPI_INT, mpi4doubles) );

  TEST_ASSERT( mpiDatatypeIsSame (mpi3doubles, mpi3doubles) );
  TEST_ASSERT( mpiDatatypeIsSame (mpi4doubles, mpi4doubles) );
  TEST_ASSERT( ! mpiDatatypeIsSame (mpi3doubles, mpi4doubles) );
  TEST_ASSERT( ! mpiDatatypeIsSame (mpi4doubles, mpi3doubles) );

  err = MPI_Type_free (&mpi3doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_free "
     "returned err != MPI_SUCCESS.");
  err = MPI_Type_free (&mpi4doubles);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::logic_error, "MPI_Type_free "
     "returned err != MPI_SUCCESS.");
}


TEUCHOS_UNIT_TEST( MpiTypeTraits, CompareWithRawMpi )
{
  using Teuchos::Details::MpiTypeTraits;
  using std::endl;

  out << "Compare MpiTypeTraits<T>::getType result with known result" << endl;
  Teuchos::OSTab tab1 (out);
  out << "Test one-argument version of MpiTypeTraits<T>::getType" << endl;

  MPI_Datatype dt_char = MpiTypeTraits<char>::getType (static_cast<char> ('a'));

  MPI_Datatype dt_short =
    MpiTypeTraits<short>::getType (static_cast<short> (42));
  MPI_Datatype dt_unsigned_short =
    MpiTypeTraits<unsigned short>::getType (static_cast<unsigned short> (42));

  MPI_Datatype dt_int =
    MpiTypeTraits<int>::getType (static_cast<int> (42));
  MPI_Datatype dt_unsigned_int =
    MpiTypeTraits<unsigned int>::getType (static_cast<unsigned int> (42));

  MPI_Datatype dt_long =
    MpiTypeTraits<long>::getType (static_cast<long> (42));
  MPI_Datatype dt_unsigned_long =
    MpiTypeTraits<unsigned long>::getType (static_cast<unsigned long> (42));

  MPI_Datatype dt_float =
    MpiTypeTraits<float>::getType (static_cast<float> (4.2));
  MPI_Datatype dt_double =
    MpiTypeTraits<double>::getType (static_cast<double> (4.2));

  // Sanity check.
  TEST_ASSERT( ! mpiDatatypeIsSame (dt_double, dt_int) );

  TEST_ASSERT( mpiDatatypeIsSame (dt_char, MPI_CHAR) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_short, MPI_SHORT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_short, MPI_UNSIGNED_SHORT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_int, MPI_INT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_int, MPI_UNSIGNED) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_long, MPI_LONG) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_long, MPI_UNSIGNED_LONG) );

  TEST_ASSERT( mpiDatatypeIsSame (dt_float, MPI_FLOAT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_double, MPI_DOUBLE) );

  out << "Test zero-argument version of MpiTypeTraits<T>::getType, "
    "for types T where that version is known to exist" << endl;

  dt_char = MpiTypeTraits<char>::getType ();
  dt_short = MpiTypeTraits<short>::getType ();
  dt_unsigned_short = MpiTypeTraits<unsigned short>::getType ();
  dt_int = MpiTypeTraits<int>::getType ();
  dt_unsigned_int = MpiTypeTraits<unsigned int>::getType ();
  dt_long = MpiTypeTraits<long>::getType ();
  dt_unsigned_long = MpiTypeTraits<unsigned long>::getType ();
  dt_float = MpiTypeTraits<float>::getType ();
  dt_double = MpiTypeTraits<double>::getType ();

  // Sanity check.
  TEST_ASSERT( ! mpiDatatypeIsSame (dt_double, dt_int) );

  TEST_ASSERT( mpiDatatypeIsSame (dt_char, MPI_CHAR) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_short, MPI_SHORT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_short, MPI_UNSIGNED_SHORT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_int, MPI_INT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_int, MPI_UNSIGNED) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_long, MPI_LONG) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_unsigned_long, MPI_UNSIGNED_LONG) );

  TEST_ASSERT( mpiDatatypeIsSame (dt_float, MPI_FLOAT) );
  TEST_ASSERT( mpiDatatypeIsSame (dt_double, MPI_DOUBLE) );

  TEST_ASSERT( MpiTypeTraits<char>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<short>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<unsigned short>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<int>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<unsigned int>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<long>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<unsigned long>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<float>::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<double>::isSpecialized );

  TEST_ASSERT( ! MpiTypeTraits<char>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<short>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<unsigned short>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<int>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<unsigned int>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<long>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<unsigned long>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<float>::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<double>::needsFree );

#ifdef HAVE_TEUCHOS_COMPLEX

  TEST_ASSERT( MpiTypeTraits<std::complex<float> >::isSpecialized );
  TEST_ASSERT( MpiTypeTraits<std::complex<double> >::isSpecialized );

  MPI_Datatype dt_complex_float =
    MpiTypeTraits<std::complex<float> >::getType ();
  MPI_Datatype dt_complex_double =
    MpiTypeTraits<std::complex<float> >::getType ();

// We make no promises about this for MPI_VERSION < 3.
#if MPI_VERSION >= 3
  TEST_ASSERT( ! MpiTypeTraits<std::complex<float> >::needsFree );
  TEST_ASSERT( ! MpiTypeTraits<std::complex<double> >::needsFree );
#endif // MPI_VERSION >= 3

  if (MpiTypeTraits<std::complex<float> >::needsFree) {
    (void) MPI_Type_free (&dt_complex_float);
  }
  if (MpiTypeTraits<std::complex<double> >::needsFree) {
    (void) MPI_Type_free (&dt_complex_double);
  }
#endif // HAVE_TEUCHOS_COMPLEX
}

#endif // MPI_VERSION >= 2

} // namespace (anonymous)
