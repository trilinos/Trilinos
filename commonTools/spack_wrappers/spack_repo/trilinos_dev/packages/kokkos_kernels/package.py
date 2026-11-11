from spack_repo.builtin.packages.kokkos_kernels.package import KokkosKernels as BuiltinKokkosKernels

from spack.package import *


class KokkosKernels(BuiltinKokkosKernels):

    version("5.2.2", sha256="4b95ac27319714a393b9a9ddf3a77fb7ffe751fe215a02084dd781856814d55f")
    version("5.2.1", sha256="346068d9c3bf25293f0460a8edc16cb02edf327d2c4197840ecd26ac37844df7")
    version("5.2.0", sha256="57f278e7a25940f9a1237d8b370a1de7f09ea282626b4594e4880042181acdcd")
    version("5.1.1", sha256="4415c2a6e14e2bba9aa978917d2ffbcbe32760d3aba3a33bf7a267d50e7e20c9")
    version("5.1.0", sha256="c003cd53126dee651f41a3b003e443950c3030246785ae968055ce015c89e0d5")
    version("5.0.2", sha256="7c7af2c3659ecc620cc7b7016876330d9f288e8c6fd7b70b70907687df823b43")
    version("5.0.1", sha256="c9d0b507ab754b347b71b530683e9dea8dbf4d2e3cdadb863dcb494b07bbf0b2")
    version("5.0.0", sha256="e1d7e7040b86f141004126c3fa5711f28697803d921c3558d82192a32156b1b2")
    version("4.7.04", sha256="dc357de6f57c42ad6dd45228d38da65dc835788f66309aaed6feeb65cae62210")
    version("4.7.03", sha256="902701e3481b2c535925a02e0918baed6f82186a25297e45314712b0905d9780")
    version("4.7.02", sha256="2d3b3e10ac112e382b88c50d66e4222ce543fca5d42be6d8376a684b82b8b238")

    depends_on("kokkos@5.2.2", when="@5.2.2")
    depends_on("kokkos@5.2.1", when="@5.2.1")
    depends_on("kokkos@5.2.0", when="@5.2.0")
    depends_on("kokkos@5.1.1", when="@5.1.1")
    depends_on("kokkos@5.0.2", when="@5.0.2")
    depends_on("kokkos@5.0.1", when="@5.0.1")
    depends_on("kokkos@5.0.0", when="@5.0.0")
    depends_on("kokkos@4.7.04", when="@4.7.04")
    depends_on("kokkos@4.7.03", when="@4.7.03")
    depends_on("kokkos@4.7.02", when="@4.7.02")
