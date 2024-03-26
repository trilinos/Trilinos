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

KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_LAPACK"
        OFF
        BOOL
        "Whether to build the lapack component. Default: OFF"
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
KOKKOSKERNELS_ADD_OPTION(
       "ENABLE_COMPONENT_ODE"
        OFF
        BOOL
        "Whether to build the ode component. Default: OFF"
)


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
  SET(KokkosKernels_ENABLE_COMPONENT_LAPACK ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_GRAPH ON CACHE BOOL "" FORCE)
ENDIF()

# If user requested to enable all components, enable all components
IF (KokkosKernels_ENABLE_ALL_COMPONENTS)
  SET(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_LAPACK ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_SPARSE ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_GRAPH ON CACHE BOOL "" FORCE)
  SET(KokkosKernels_ENABLE_COMPONENT_ODE ON CACHE BOOL "" FORCE)
ENDIF()

# KOKKOSKERNELS_ALL_COMPONENTS_ENABLED says whether all components are on,
# regardless of how this came to be
# this is in the cache so we can use it as a global variable,
# but marking it as advanced should hide it from GUIs
IF (    KokkosKernels_ENABLE_COMPONENT_BATCHED 
    AND KokkosKernels_ENABLE_COMPONENT_BLAS
    AND KokkosKernels_ENABLE_COMPONENT_LAPACK
    AND KokkosKernels_ENABLE_COMPONENT_GRAPH 
    AND KokkosKernels_ENABLE_COMPONENT_SPARSE
    AND KokkosKernels_ENABLE_COMPONENT_ODE)
  SET(KOKKOSKERNELS_ALL_COMPONENTS_ENABLED ON CACHE BOOL "" FORCE)
ELSE()
  SET(KOKKOSKERNELS_ALL_COMPONENTS_ENABLED OFF CACHE BOOL "" FORCE)
ENDIF()
mark_as_advanced(FORCE KOKKOSKERNELS_ALL_COMPONENTS_ENABLED)