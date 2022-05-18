#pragma once

#include <Kokkos_Core.hpp>

namespace Tpetra {

    template <typename Space>
    Space &get_space(int i);

    // Kokkos::DefaultExecutionSpace &get_exec_space(int i);
    void drop_exec_spaces();





    template<>
    Kokkos::Serial &get_space(int i);

    template<>
    Kokkos::Cuda &get_space(int i);

}