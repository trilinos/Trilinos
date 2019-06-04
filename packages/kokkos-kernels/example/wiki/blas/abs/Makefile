MAKEFILE_PATH := $(subst Makefile,,$(abspath $(lastword $(MAKEFILE_LIST))))

ifndef KOKKOSKERNELS_PATH
  KOKKOSKERNELS_PATH = $(MAKEFILE_PATH)../../../..
endif
ifndef KOKKOS_PATH
  KOKKOS_PATH = $(KOKKOSKERNELS_PATH)/../kokkos
endif

SRC = $(wildcard $(MAKEFILE_PATH)*.cpp)
HDR = $(wildcard $(MAKEFILE_PATH)*.hpp)

include $(KOKKOSKERNELS_PATH)/example/wiki/Makefile

