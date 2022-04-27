#pragma once

#include <Kokkos_Core.hpp>

namespace Tpetra {

    Kokkos::DefaultExecutionSpace &get_exec_space(int i);
    void drop_exec_spaces();

}