from spack_repo.builtin.packages.kokkos.package import Kokkos as BuiltinKokkos

from spack.package import *


class Kokkos(BuiltinKokkos):

    version("5.2.2", sha256="d6557aaef39302282a15f9c770433d1fcdf4e961dfd6d9ed726b9d0d0f546b9f")
    version("5.2.1", sha256="3f754c99aa6130b1dd6520d904db7b2fd44ed618cd91e0dfd921956f23f6812d")
    version("5.2.0", sha256="54993e0682d80b78939bbf260490f8cf31428bb883c0309961369997f15d94df")
    version("5.2.2", sha256="d6557aaef39302282a15f9c770433d1fcdf4e961dfd6d9ed726b9d0d0f546b9f")
    version("5.2.1", sha256="3f754c99aa6130b1dd6520d904db7b2fd44ed618cd91e0dfd921956f23f6812d")
    version("5.2.0", sha256="54993e0682d80b78939bbf260490f8cf31428bb883c0309961369997f15d94df")
    version("5.1.1", sha256="8bdbee0f0ac383436743ad8a9e3e928705b34b31a25a92dc5179c52a3aa98519")
    version("5.1.0", sha256="7bdbdfc88033ed7d940c7940ed8919e1f2b78a9656c69276beb76ad45c41ec4e")
    version("5.0.2", sha256="188817bb452ca805ee8701f1c5adbbb4fb83dc8d1c50624566a18a719ba0fa5e")
    version("5.0.1", sha256="cf7d8515ca993229929be9f051aecd8f93cde325adac8a4f82ed6848adace218")
    version("5.0.0", sha256="c45f3e19c3eb71fc8b7210cb04cac658015fc1839e7cc0571f7406588ff9bcef")
    version("4.7.04", sha256="4213b248c39e112299fa94ee08817e51126fc02996ed6e2ab56aec4cdb80ee1f")
    version("4.7.03", sha256="969e7933b9426219b220f08036e489b3226e6d8cd24eecf2c5b80df8c37443c0")
    version("4.7.02", sha256="a81826ac0a167933d13506bc2a986fb5517038df9abb780fe9bb2c1d4e80803b")
    
    depends_on("googletest@1.10.0:", type="test")

    patch("fix_mdspan_debug_logic.patch", when="@5.0:5.2.0")
    patch("gnu12nullptr.patch", when="@:5.1.1 %cxx=gcc@12.5")
