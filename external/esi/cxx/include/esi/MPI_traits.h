#ifndef __ESI_MPI_traits_h
#define __ESI_MPI_traits_h

namespace esi {

/** The ESI MPI_traits file. */

template<class T>
struct MPI_traits {};

template<>
struct MPI_traits<float> {
  static MPI_Datatype mpi_type() {return(MPI_FLOAT);};
};

template<>
struct MPI_traits<double> {
  static MPI_Datatype mpi_type() {return(MPI_DOUBLE);};
};

template<>
struct MPI_traits<int> {
  static MPI_Datatype mpi_type() {return(MPI_INT);};
};

template<>
struct MPI_traits<long> {
  static MPI_Datatype mpi_type() {return(MPI_LONG);};
};

#ifndef ESI_NO_COMPLEX
//temporarily close the esi namespace:
};
#include <complex>
//reopen the esi namespace:
namespace esi {

template<>
struct MPI_traits<std::complex<float> > {
  static MPI_Datatype mpi_type() {return(MPI_COMPLEX);};
};

template<>
struct MPI_traits<std::complex<double> > {
  static MPI_Datatype mpi_type() {return(MPI_DOUBLE_COMPLEX);};
};

#endif /* ESI_NO_COMPLEX */

};     // esi namespace

#endif
