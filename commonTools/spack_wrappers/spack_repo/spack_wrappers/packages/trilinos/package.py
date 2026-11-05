from spack_repo.builtin.packages.trilinos.package import Trilinos as BuiltinTrilinos

from spack.package import *


class Trilinos(BuiltinTrilinos):

    depends_on("kokkos@5.2.2", when='@develop')
    depends_on("kokkos-kernels@5.2.2", when='@develop')
