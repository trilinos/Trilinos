#ifndef TPETRA_SPACE_HPP
#define TPETRA_SPACE_HPP

#include <Kokkos_Core.hpp>

namespace Tpetra {

    template <typename Space>
    Space &get_space(int i);

    void drop_exec_spaces();


#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    Kokkos::Serial &get_space(int i);
#endif


#ifdef KOKKOS_ENABLE_CUDA
    template<>
    Kokkos::Cuda &get_space(int i);
#endif


#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    Kokkos::OpenMP &get_space(int i);
#endif
}

#endif // TPETRA_SPACE_HPP
