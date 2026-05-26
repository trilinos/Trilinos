source /home/runner/.bashrc
export TRILINOS_DIR=$(realpath ../..)
source ${TRILINOS_DIR}/packages/framework/GenConfig/gen-config.sh \
  rhel_clang-openmpi_release-debug_shared_no-kokkos-arch_no-asan_no-complex_no-fpic_mpi_no-pt_no-rdc_no-uvm_deprecated-on_no-package-enables \
  --cmake-fragment GenConfigSettings.cmake \
  --force -y \
  "$@" \
  ${TRILINOS_DIR}
