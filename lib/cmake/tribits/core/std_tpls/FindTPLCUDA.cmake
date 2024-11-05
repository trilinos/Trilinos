# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

find_package(CUDAToolkit REQUIRED)  # Will abort if not found!
tribits_extpkg_create_imported_all_libs_target_and_config_file( CUDA
  INNER_FIND_PACKAGE_NAME  CUDAToolkit
  IMPORTED_TARGETS_FOR_ALL_LIBS  CUDA::cufft  CUDA::cublas  CUDA::cudart CUDA::cuda_driver )
# Above, we could add more dependencies if we need
