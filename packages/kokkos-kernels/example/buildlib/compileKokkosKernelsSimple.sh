KOKKOS_PATH="${HOME}/Kokkos/kokkos" #path to kokkos source
KOKKOSKERNELS_PATH="../.."          #path to kokkos-kernels top directory.

KOKKOSKERNELS_SCALARS=double #the scalar types to instantiate =double,float...
KOKKOSKERNELS_LAYOUTS=LayoutLeft #the layout types to instantiate.
KOKKOSKERNELS_ORDINALS=int #ordinal types to instantiate
KOKKOSKERNELS_OFFSETS=int #offset types to instantiate
CXX=${KOKKOS_PATH}/bin/nvcc_wrapper
KOKKOSKERNELS_OPTIONS=eti-only #options for kokkoskernels  
KOKKOS_DEVICES=Cuda
KOKKOS_ARCHS=SKX,Volta70
KOKKOS_CUDA_OPTIONS=enable_lambda
CXXFLAGS="-Wall -pedantic -Werror -O3 -g -Wshadow -Wsign-compare -Wtype-limits -Wuninitialized"

../../cm_generate_makefile.bash --kokkoskernels-path=${KOKKOSKERNELS_PATH} --with-scalars=${KOKKOSKERNELS_SCALARS} --with-ordinals=${KOKKOSKERNELS_ORDINALS} --with-offsets=${KOKKOSKERNELS_OFFSETS} --kokkos-path=${KOKKOS_PATH} --with-devices=${KOKKOS_DEVICES} --arch=${KOKKOS_ARCHS} --compiler=${CXX} --with-cuda-options=${KOKKOS_CUDA_OPTIONS} --with-options=${KOKKOSKERNELS_OPTIONS} --cxxflags="${CXXFLAGS}"

# Call "../../scripts/cm_generate_makefile.bash --help" for options
