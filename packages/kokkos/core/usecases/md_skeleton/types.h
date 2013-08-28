/*
 * types.h
 *
 *  Created on: Aug 27, 2013
 *      Author: crtrott
 */

#ifndef TYPES_H_
#define TYPES_H_

/* Determine default device type and necessary includes */

#ifdef _OPENMP
  #include <Kokkos_OpenMP.hpp>
#else
  #include <Kokkos_Threads.hpp>
#endif

#ifdef KOKKOS_HAVE_CUDA
  #include <Kokkos_Cuda.hpp>
  #include <cuda.h>
  #include <cuda_runtime.h>
  typedef Kokkos::Cuda device_type;
#else
  #ifdef _OPENMP
    typedef Kokkos::OpenMP device_type;
  #else
    typedef Kokkos::Threads device_type;
  #endif

  struct double2 {
    double x, y;
  };
#endif


#include <Kokkos_View.hpp>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Timer.hpp>

/* Define types used throughout the code */

//Position arrays
typedef Kokkos::View<double*[3], Kokkos::LayoutRight, device_type>                                   t_x_array ;
typedef t_x_array::HostMirror                                                                        t_x_array_host ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, device_type>                             t_x_array_const ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, device_type, Kokkos::MemoryRandomRead >  t_x_array_randomread ;

//Force array
typedef Kokkos::View<double*[3],  device_type>                                                       t_f_array ;


//Neighborlist
typedef Kokkos::View<int**, device_type >                                                            t_neighbors ;
typedef Kokkos::View<const int**, device_type >                                                      t_neighbors_const ;
typedef Kokkos::View<int*, device_type, Kokkos::MemoryUnmanaged >                                    t_neighbors_sub ;
typedef Kokkos::View<const int*, device_type, Kokkos::MemoryUnmanaged >                              t_neighbors_const_sub ;

//1d int array
typedef Kokkos::View<int*, device_type >                                                             t_int_1d ;
typedef t_int_1d::HostMirror                                                                         t_int_1d_host ;
typedef Kokkos::View<const int*, device_type >                                                       t_int_1d_const ;
typedef Kokkos::View<int*, device_type , Kokkos::MemoryUnmanaged>                                    t_int_1d_um ;
typedef Kokkos::View<const int* , device_type , Kokkos::MemoryUnmanaged>                             t_int_1d_const_um ;

//2d int array
typedef Kokkos::View<int**, Kokkos::LayoutRight, device_type >                                       t_int_2d ;
typedef t_int_2d::HostMirror                                                                         t_int_2d_host ;

//Scalar ints
typedef Kokkos::View<int[1], Kokkos::LayoutLeft, device_type>                                        t_int_scalar ;
typedef t_int_scalar::HostMirror                                                                     t_int_scalar_host ;

#endif /* TYPES_H_ */
