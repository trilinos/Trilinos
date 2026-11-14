from spack_repo.builtin.packages.trilinos.package import Trilinos as BuiltinTrilinos

from spack.package import *


class Trilinos(BuiltinTrilinos):


    version("17.3.0", sha256="8094740da60bf75b675ad556ed00e337c8547be0dc9603893c4894f4abf0396b")
    version("17.2.1", sha256="b1d076c9c474893ddb340f0e2418a2dec50f52edaf134b54e4c01c5c53032eae")
    version("17.2.0", sha256="3df928f8d600fd8952593565c507112bc562ac32f3ef9ba6217d9a86072c3e72")
    version("17.1.2", sha256="f6ebacae0fd09c096eb7ff35c45c455525e22e4d317e98ceee5bbcfd880dc353")
    version("17.1.1", sha256="482bb76c9a993767a9e649261462f70620254a631d34fb089482d8474dfb2766")
    version("17.0.0", sha256="7afa68fc6bf1dfdcd0c07f7b61055b03509e62cee1a835d570201b46aa440a6b")
    version("16.2.2", sha256="211810d0f642736ac93ca500c3e28c72a9e8fa9991a757d60c0ecc85108e35d4")
    version("16.2.1", sha256="c68a9d28fc9e7b06f33804c1f5e998234820878c21ae075649483b56992cea05")
    version("16.2.0", sha256="543aa56232d7c0cbe73705fab2d3b5524f11b15fef8917aa14de02d23a5ca418")
    
    depends_on("kokkos@=5.2.2", when="@master:")
    depends_on("kokkos-kernels@=5.2.2", when="@master:")
    depends_on("kokkos@5.2.2", when='@develop')
    depends_on("kokkos-kernels@5.2.2", when='@develop')
    depends_on("kokkos@=5.2.2", when="@17.3.0")
    depends_on("kokkos-kernels@=5.2.2", when="@17.3.0")
    depends_on("kokkos@=5.2.1", when="@17.2.1")
    depends_on("kokkos-kernels@=5.2.1", when="@17.2.1")
    depends_on("kokkos@=5.2.1", when="@17.2.0")
    depends_on("kokkos-kernels@=5.2.1", when="@17.2.0")
    depends_on("kokkos@=5.1.1", when="@17.1.2")
    depends_on("kokkos-kernels@=5.1.1", when="@17.1.2")
    depends_on("kokkos@=5.1.1", when="@17.1.1")
    depends_on("kokkos-kernels@=5.1.1", when="@17.1.1")
    depends_on("kokkos@=5.1.1", when="@17.1.0")
    depends_on("kokkos-kernels@=5.1.1", when="@17.1.0")
    depends_on("kokkos@=5.0.2", when="@17.0")
    depends_on("kokkos-kernels@=5.0.2", when="@17.0")
    depends_on("kokkos@=4.7.04", when="@16.2.2")
    depends_on("kokkos-kernels@=4.7.04", when="@16.2.2")
    depends_on("kokkos@=4.7.04", when="@16.2.1")
    depends_on("kokkos-kernels@=4.7.04", when="@16.2.1")
    depends_on("kokkos@=4.7.01", when="@16.2.0")
    depends_on("kokkos-kernels@=4.7.01", when="@16.2.0")
