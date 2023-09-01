#
# Define Some CMake and TriBiTS tweaks just for Trilinos (a.k.a. hacks)
#

# Define empty macro to allow Kokkos and KokkosKernels to keep working until
# updated versions after Kokkos 4.1 are snapshotted in.
macro(tribits_exclude_autotools_files)
endmacro()
