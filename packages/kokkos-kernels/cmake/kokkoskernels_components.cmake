# Define component dependencies and enable
# them selectively based on what the user
# requests.

KOKKOSKERNELS_ADD_OPTION(
        "ENABLE_ALL_COMPONENTS"
        ON
        BOOL
        "Whether to build all the library's components. Default: ON"
)

# BATCHED only depends on COMMON which
# is always enabled so nothing more needs
# to be enabled for this component.
KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_BATCHED"
        OFF
        BOOL
        "Whether to build the batched component. Default: OFF"
)

# BLAS only depends on COMMON which
# is always enabled so nothing more needs
# to be enabled for this component.
KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_BLAS"
        OFF
        BOOL
        "Whether to build the blas component. Default: OFF"
)

# SPARSE depends on everything else at the moment.
KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_SPARSE"
        OFF
        BOOL
        "Whether to build the sparse component. Default: OFF"
)

# GRAPH depends on everything else at the moment.
KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_GRAPH"
        OFF
        BOOL
        "Whether to build the graph component. Default: OFF"
)

# The user requested individual components,
# the assumption is that a full build is not
# desired and ENABLE_ALL_COMPONENETS is turned
# off.
IF (KokkosKernels_ENABLE_COMPONENT_BATCHED OR KokkosKernels_ENABLE_COMPONENT_BLAS
    OR KokkosKernels_ENABLE_COMPONENT_GRAPH OR KokkosKernels_ENABLE_COMPONENT_SPARSE)
  SET(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "" FORCE)
ENDIF()

# Graph depends on everything else because it depends
# on Sparse at the moment, breaking that dependency will
# allow for Graph to build pretty much as a standalone
# component.
IF (KokkosKernels_ENABLE_COMPONENT_GRAPH)
  SET(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_SPARSE ON CACHE BOOL "" FORCE)
ENDIF()

# Sparse depends on everything else so no real benefit
# here unfortunately...
IF (KokkosKernels_ENABLE_COMPONENT_SPARSE)
  SET(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_GRAPH ON CACHE BOOL "" FORCE)
ENDIF()

# At this point, if ENABLE_ALL_COMPONENTS is
# still ON we need to enable all individual
# components as they are required for this
# build.
IF (KokkosKernels_ENABLE_ALL_COMPONENTS)
  SET(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_SPARSE ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_GRAPH ON CACHE BOOL "" FORCE)
ENDIF()
