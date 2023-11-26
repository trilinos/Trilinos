ADD_SUBDIRECTORY(CG)
ADD_SUBDIRECTORY(cusolver)
ADD_SUBDIRECTORY(GMRES)
ADD_SUBDIRECTORY(SPMV)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_CG.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/run_CG.sh
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_GMRES.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/run_GMRES.sh
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_SPMV.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/run_SPMV.sh
)
