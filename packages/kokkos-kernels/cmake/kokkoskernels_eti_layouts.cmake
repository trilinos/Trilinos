SET(RIGHT_LAYOUTS
  LAYOUTRIGHT)
SET(LEFT_LAYOUTS
  LAYOUTLEFT)
SET(LAYOUTS
  LAYOUTLEFT
  LAYOUTRIGHT)
SET(LAYOUTLEFT_CPP_TYPE Kokkos::LayoutLeft)
SET(LAYOUTRIGHT_CPP_TYPE Kokkos::LayoutRight)

KOKKOSKERNELS_ADD_OPTION(
  INST_LAYOUTLEFT
  ${KOKKOSKERNELS_ADD_DEFAULT_ETI}
  BOOL
  "Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is KokkosKernels_INST_LAYOUTLEFT=ON by default.  Disabling this may increase build times. Default: ON"
  )

KOKKOSKERNELS_ADD_OPTION(
  INST_LAYOUTRIGHT
  OFF
  BOOL
  "Whether to pre instantiate kernels for the view layout LayoutRight.  This option is KokkosKernels_INST_LAYOUTRIGHT=OFF by default.  Disabling this may increase build times. Default: OFF"
  )

IF (KOKKOSKERNELS_INST_LAYOUTLEFT)
  LIST(APPEND LAYOUT_LIST "LayoutLeft")
ENDIF()

IF (KOKKOSKERNELS_INST_LAYOUTRIGHT)
  LIST(APPEND LAYOUT_LIST "LayoutRight")
ENDIF()
