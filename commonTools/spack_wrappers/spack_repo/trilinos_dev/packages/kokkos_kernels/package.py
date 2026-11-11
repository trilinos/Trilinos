from spack_repo.builtin.packages.kokkos_kernels.package import KokkosKernels as BuiltinKokkosKernels

from spack.package import *


class KokkosKernels(BuiltinKokkosKernels):

    version("5.2.2", sha256="4b95ac27319714a393b9a9ddf3a77fb7ffe751fe215a02084dd781856814d55f")
    version("5.2.1", sha256="346068d9c3bf25293f0460a8edc16cb02edf327d2c4197840ecd26ac37844df7")
    version("5.2.0", sha256="57f278e7a25940f9a1237d8b370a1de7f09ea282626b4594e4880042181acdcd")
    depends_on("kokkos@5.2.2", when="@5.2.2")
    depends_on("kokkos@5.2.1", when="@5.2.1")
    depends_on("kokkos@5.2.0", when="@5.2.0")
