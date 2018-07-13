set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Enable/Disable output of compile commands during generation." FORCE)

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/.ycm_extra_conf.py.in
    ${CMAKE_SOURCE_DIR}/.ycm_extra_conf.py
    @ONLY
)
