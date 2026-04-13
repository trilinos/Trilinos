set(RIGHT_LAYOUTS         LAYOUTRIGHT)
set(LEFT_LAYOUTS          LAYOUTLEFT)
set(LAYOUTS               LAYOUTLEFT LAYOUTRIGHT)
set(LAYOUTLEFT_CPP_TYPE   Kokkos::LayoutLeft)
set(LAYOUTRIGHT_CPP_TYPE  Kokkos::LayoutRight)

kokkoskernels_add_option("INST_LAYOUTLEFT" ${KOKKOSKERNELS_ADD_DEFAULT_ETI} BOOL
  "Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is KokkosKernels_INST_LAYOUTLEFT=ON by default.  Disabling this may increase build times. Default: ON")

kokkoskernels_add_option("INST_LAYOUTRIGHT" OFF BOOL
  "Whether to pre instantiate kernels for the view layout LayoutRight.  This option is KokkosKernels_INST_LAYOUTRIGHT=OFF by default.  Disabling this may increase build times. Default: OFF")

if(KOKKOSKERNELS_INST_LAYOUTLEFT)
  list(APPEND LAYOUT_LIST "LayoutLeft")
endif()

if(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  list(APPEND LAYOUT_LIST "LayoutRight")
endif()
