# Disable Sacado GTest tests with this intel-18 compiler to avoid build error
# (#7778)
ATDM_SET_ENABLE(Sacado_ENABLE_GTest OFF)
