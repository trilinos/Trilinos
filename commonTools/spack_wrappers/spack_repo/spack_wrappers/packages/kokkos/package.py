from spack_repo.builtin.packages.kokkos.package import Kokkos as BuiltinKokkos

from spack.package import *


class Kokkos(BuiltinKokkos):

    version("5.2.2", sha256="d6557aaef39302282a15f9c770433d1fcdf4e961dfd6d9ed726b9d0d0f546b9f")
    version("5.2.1", sha256="3f754c99aa6130b1dd6520d904db7b2fd44ed618cd91e0dfd921956f23f6812d")
    version("5.2.0", sha256="54993e0682d80b78939bbf260490f8cf31428bb883c0309961369997f15d94df")
    
    #depends_on("googletest@1.10.0:", when="+tests")

    patch("fix_mdspan_debug_logic.patch", when="@5.0:5.2.0")
    patch("gnu12nullptr.patch", when="@:5.1.1 %cxx=gcc@12.5")
