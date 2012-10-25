/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_mpiTraits_hpp_
#define _fei_mpiTraits_hpp_

#include <fei_macros.hpp>
#include <fei_defs.h>

#ifndef FEI_SER

namespace fei {

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<class T>
    struct mpiTraits {};

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<char> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_CHAR);};
    };

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<float> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_FLOAT);};
    };

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<double> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_DOUBLE);};
    };

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<int> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_INT);};
    };

  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<long> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_LONG);};
    };

#ifdef EIGHT_BYTE_GLOBAL_ID
  /** Simple trait struct for querying the mpi-data-type of
   an intrinsic type. */
  template<>
    struct mpiTraits<GlobalID> {
      /** mpi type query */
      static MPI_Datatype mpi_type() {return(MPI_LONG_LONG);};
    };
#endif

} //namespace fei

#endif //FEI_SER

#endif // _fei_mpiTraits_hpp_

