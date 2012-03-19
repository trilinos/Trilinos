SET(${PROJECT_NAME}_PACKAGES
  CTeuchos
  ForTeuchos
  Coupler
  DataTransferKit
  STARCCM
  DeCART
  Nemesis
  Denovo
  CASLBOA
  CASLRAVE
  RELAP5
  LIME
  VRIPSS
  )

# Force Panzer, Drekar and other MPI-dependent code to off
SET(TPL_ENABLE_MPI OFF)
