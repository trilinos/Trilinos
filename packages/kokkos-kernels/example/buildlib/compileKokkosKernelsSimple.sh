KOKKOS_PATH=${HOME}/proj/kokkos #path to kokkos source
KOKKOSKERNELS_SCALARS=double #the scalar types to instantiate =double,float...
KOKKOSKERNELS_LAYOUTS=LayoutLeft #the layout types to instantiate.
KOKKOSKERNELS_ORDINALS=int #ordinal types to instantiate
KOKKOSKERNELS_OFFSETS=int #offset types to instantiate
KOKKOSKERNELS_PATH=../.. #path to kokkos-kernels top directory.
CXX=icpc #${KOKKOS_PATH}/config/nvcc_wrapper #icpc #
KOKKOSKERNELS_OPTIONS=eti-only #options for kokkoskernels  
KOKKOS_DEVICES=OpenMP # other devices Cuda,Serial ..
KOKKOS_ARCHS=KNL
CXXFLAGS="-Wall -pedantic -Werror -O3 -g -Wshadow -Wsign-compare -Wtype-limits -Wuninitialized"

../../scripts/generate_makefile.bash --kokkoskernels-path=${KOKKOSKERNELS_PATH} --with-scalars=${KOKKOSKERNELS_SCALARS} --with-ordinals=${KOKKOSKERNELS_ORDINALS} --with-offsets=${KOKKOSKERNELS_OFFSETS} --kokkos-path=${KOKKOS_PATH} --with-devices=${KOKKOS_DEVICES} --arch=${KOKKOS_ARCHS} --compiler=${CXX} --with-options=${KOKKOSKERNELS_OPTIONS}  --cxxflags="${CXXFLAGS}"
