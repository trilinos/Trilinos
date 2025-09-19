#@HEADER
# ************************************************************************
#
#                        Kokkos v. 4.0
#       Copyright (2022) National Technology & Engineering
#               Solutions of Sandia, LLC (NTESS).
#
# Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
# See https://kokkos.org/LICENSE for license information.
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
#@HEADER

"""Check if an API change in Kokkos Kernels might be undocumented

This script is intended to run in a github action when a Pull Request is opened against 
Kokkos Kernels' develop branch. It will read the output of git diff --name-only, stored locally 
in 'modified_files.txt' to get a list of modified files. If one of the modified files is a key 
in src_doc_mapping it checks if one of the corresponding documentation files has also been 
modified which indicates some amount of documentation was done. Otherwise the script will
return with exit core 1. For API files stored in one of the 'src_directories' but not listed
in 'src_doc_mapping' an exit code 1 is also returned.
"""

import sys

src_directories = ["batched/dense/src",
                   "batched/sparse/src",
                   "blas/src",
                   "common/src",
                   "graph/src",
                   "lapack/src",
                   "ode/src",
                   "sparse/src"]

src_doc_mapping = dict([
                        ('blas/src/KokkosBlas1_abs.hpp', ['blas1_abs.rst']),
                        ('blas/src/KokkosBlas1_axpby.hpp', ['blas1_axpy.rst']),
                        ('blas/src/KokkosBlas1_dot.hpp', ['blas1_dot.rst']),
                        # ('blas/src/KokkosBlas1_fill.hpp', ),
                        ('blas/src/KokkosBlas1_iamax.hpp', ['blas1_iamax.rst']),
                        # ('blas/src/KokkosBlas1_mult.hpp', ),
                        ('blas/src/KokkosBlas1_nrm1.hpp', ['blas1_nrm1.rst']),
                        ('blas/src/KokkosBlas1_nrm2.hpp', ['blas1_nrm2.rst']),
                        # ('blas/src/KokkosBlas1_nrm2_squared.hpp', ),
                        # ('blas/src/KokkosBlas1_nrm2w.hpp', ),
                        # ('blas/src/KokkosBlas1_nrm2w_squared.hpp', ),
                        # ('blas/src/KokkosBlas1_nrminf.hpp', ['']),
                        # ('blas/src/KokkosBlas1_reciprocal.hpp', ['']),
                        ('blas/src/KokkosBlas1_rot.hpp', ['blas1_rot.rst']),
                        ('blas/src/KokkosBlas1_rotg.hpp', ['blas1_rotg.rst']),
                        ('blas/src/KokkosBlas1_rotm.hpp', ['blas1_rotm.rst']),
                        ('blas/src/KokkosBlas1_rotmg.hpp', ['blas1_rotmg.rst']),
                        ('blas/src/KokkosBlas1_scal.hpp', ['blas1_scal.rst']),
                        # ('blas/src/KokkosBlas1_set.hpp', ['']),
                        # ('blas/src/KokkosBlas1_sum.hpp', ['']),
                        ('blas/src/KokkosBlas1_swap.hpp', ['blas1_swap.rst']),
                        # ('blas/src/KokkosBlas1_update.hpp', ['']),
                        ('blas/src/KokkosBlas2_gemv.hpp', ['blas2_gemv.rst']),
                        ('blas/src/KokkosBlas2_ger.hpp', ['blas2_ger.rst']),
                        ('blas/src/KokkosBlas2_syr.hpp', ['blas2_syr.rst']),
                        ('blas/src/KokkosBlas2_syr2.hpp', ['blas2_syr2.rst']),
                        ('blas/src/KokkosBlas3_gemm.hpp', ['blas3_gemm.rst']),
                        ('blas/src/KokkosBlas3_trmm.hpp', ['blas3_trmm.rst']),
                        ('blas/src/KokkosBlas3_trsm.hpp', ['blas3_trsm.rst']),
                        ('lapack/src/KokkosLapack_gesv.hpp', ['docs/source/API/lapack/gesv.rst']),
                        ('lapack/src/KokkosLapack_svd.hpp', ['docs/source/API/lapack/gesvd.rst']),
                        ('lapack/src/KokkosLapack_trtri.hpp', ['docs/source/API/lapack/trtri.rst']),
                        ('graph/src/KokkosGraph_Distance1Color.hpp', ['docs/source/API/graph/distance1_color.rst']),
                        ('graph/src/KokkosGraph_Distance2Color.hpp', ['docs/source/API/graph/distance2_color.rst']),
                        ('graph/src/KokkosGraph_RCB.hpp', ['docs/source/API/graph/rcb.rst']),
                        ('sparse/src/KokkosKernels_Handle.hpp', ['docs/source/API/sparse/kokkoskernelshandle.rst',
                                                                 'docs/source/API/sparse/handle_get_create_destroy.rst']),
                        ('sparse/src/KokkosSparse_BsrMatrix.hpp', ['docs/source/API/sparse/bsr_matrix.rst',
                                                                   'docs/source/API/sparse/BsrMatrix_constructors.rst',
                                                                   'bsr_row_view.rst']),
                        ('sparse/src/KokkosSparse_CrsMatrix.hpp', ['docs/source/API/sparse/crs_matrix.rst',
                                                                   'docs/source/API/sparse/CrsMatrix_constructors.rst',
                                                                   'docs/source/API/sparse/sparse_row_view.rst']),
                        ('sparse/src/KokkosSparse_CcsMatrix.hpp', ['docs/source/API/sparse/ccs_matrix.rst',
                                                                   'docs/source/API/sparse/CcsMatrix_constructors.rst']),
                        ('sparse/src/KokkosSparse_CooMatrix.hpp', ['docs/source/API/sparse/coo_matrix.rst',
                                                                   'docs/source/API/sparse/CooMatrix_constructors.rst']),
                        # ('sparse/src/KokkosSparse_CcsMatrix.hpp', ['']),
                        # ('sparse/src/KokkosSparse_CooMatrix.hpp', ['']),
                        # ('sparse/src/KokkosSparse_IOUtils.hpp', []),
                        # ('sparse/src/KokkosSparse_LUPrec.hpp', []),
                        # ('sparse/src/KokkosSparse_MatrixPrec.hpp', []),
                        # ('sparse/src/KokkosSparse_OrdinalTraits.hpp', []),
                        # ('sparse/src/KokkosSparse_Preconditioner.hpp', []),
                        ('sparse/src/KokkosSparse_SortCrs.hpp', ['docs/source/API/sparse/sort_crs.rst']),
                        # ('sparse/src/KokkosSparse_StaticCcsGraph.hpp', []),
                        # ('sparse/src/KokkosSparse_StaticCrsGraph.hpp', []),
                        ('sparse/src/KokkosSparse_Utils.hpp', ['docs/source/API/sparse/extract_diagonal_blocks_rcb.rst']),
                        # ('sparse/src/KokkosSparse_Utils_cusparse.hpp', []),
                        # ('sparse/src/KokkosSparse_Utils_mkl.hpp', []),
                        # ('sparse/src/KokkosSparse_Utils_rocsparse.hpp', []),
                        # ('sparse/src/KokkosSparse_ccs2crs.hpp', []),
                        # ('sparse/src/KokkosSparse_coo2crs.hpp', []),
                        # ('sparse/src/KokkosSparse_crs2ccs.hpp', []),
                        # ('sparse/src/KokkosSparse_crs2coo.hpp', []),
                        # ('sparse/src/KokkosSparse_findRelOffset.hpp', []),
                        ('sparse/src/KokkosSparse_gauss_seidel.hpp', ['docs/source/API/sparse/gauss_seidel_apply_backward.rst',
                                                                      'docs/source/API/sparse/gauss_seidel_apply_forward.rst',
                                                                      'docs/source/API/sparse/gauss_seidel_apply_symmetric.rst',
                                                                      'docs/source/API/sparse/gauss_seidel_numeric.rst',
                                                                      'docs/source/API/sparse/gauss_seidel_symbolic.rst']),
                        ('sparse/src/KokkosSparse_gauss_seidel_handle.hpp', []),
                        # ('sparse/src/KokkosSparse_getDiagCopy.hpp', []),
                        # ('sparse/src/KokkosSparse_gmres.hpp', []),
                        # ('sparse/src/KokkosSparse_gmres_handle.hpp', []),
                        # ('sparse/src/KokkosSparse_mdf.hpp', []),
                        # ('sparse/src/KokkosSparse_mdf_handle.hpp', []),
                        ('sparse/src/KokkosSparse_par_ilut.hpp', ['docs/source/API/sparse/par_ilut.rst']),
                        ('sparse/src/KokkosSparse_par_ilut_handle.hpp', ['docs/source/API/sparse/par_ilut.rst']),
                        ('sparse/src/KokkosSparse_spadd.hpp', ['docs/source/API/sparse/spadd_numeric.rst',
                                                               'docs/source/API/sparse/spadd_symbolic.rst']),
                        # ('sparse/src/KokkosSparse_spadd_handle.hpp', [])
                        ('sparse/src/KokkosSparse_spgemm.hpp', ['docs/source/API/sparse/spgemm_numeric.rst',
                                                                'docs/source/API/sparse/spgemm_symbolic.rst']),
                        # ('sparse/src/KokkosSparse_spgemm_handle.hpp', [])
                        # ('sparse/src/KokkosSparse_spgemm_jacobi.hpp', [])
                        ('KokkosSparse_spgemm_numeric.hpp', ['docs/source/API/sparse/spgemm_numeric.rst']),
                        ('sparse/src/KokkosSparse_spgemm_symbolic.hpp', ['docs/source/API/sparse/spgemm_symbolic.rst']),
                        ('sparse/src/KokkosSparse_spiluk.hpp', ['docs/source/API/sparse/spiluk_numeric.rst',
                                                                'docs/source/API/sparse/spiluk_symbolic.rst']),
                        # ('sparse/src/KokkosSparse_spiluk_handle.hpp', []),
                        ('sparse/src/KokkosSparse_sptrsv.hpp', ['docs/source/API/sparse/sptrsv_solve.rst',
                                                                'docs/source/API/sparse/sptrsv_symbolic.rst']),
                        ('sparse/src/KokkosSparse_trsv.hpp', ['docs/source/API/sparse/sptrsv_solve.rst',
                                                              'docs/source/API/sparse/sptrsv_symbolic.rst']),
                        # ('sparse/src/KokkosSparse_sptrsv_cholmod.hpp', ['']),
                        # ('sparse/src/KokkosSparse_sptrsv_handle.hpp', [])
                        # ('sparse/src/KokkosSparse_sptrsv_superlu.hpp', [])
                        # ('sparse/src/KokkosSparse_sptrsv_supernode.hpp', [])
                        ('sparse/src/KokkosSparse_spmv.hpp', ['docs/source/API/sparse/spmv.rst'])
                        # ('sparse/src/KokkosSparse_spmv_deprecated.hpp', []),
                        # ('sparse/src/KokkosSparse_spmv_handle.hpp', [])
                        ])

with open('./modified_files.txt', 'r', encoding='utf-8') as diffs_list:
    modified_files = [line.strip() for line in diffs_list.readlines()]

modified_public_files = []
new_apis_to_document  = []
undocumented_changes  = []

# Loop over changed files
for modified_file in modified_files:
    # Check if file belongs to one of our source directories
    if any(src_dir in modified_file for src_dir in src_directories):
        modified_public_files.append(modified_file)
        # Look for documentation associated with modified file
        doc_files = src_doc_mapping.get(modified_file)
        if doc_files is None:
            new_apis_to_document.append(modified_file)
        else:
            # Construct the intersection of modified files
            # and associated documentation files that should
            # be modified. If the intersection is empty some
            # documentation is likely missing!
            intersection = set.intersection(set(doc_files), set(modified_files))
            if len(intersection) == 0:
                undocumented_changes.append(modified_file)

RETURN_VALUE = 0

print("Modified public files:")
for modified_file in modified_public_files:
    print("   "+str(modified_file))
print("")

if new_apis_to_document:
    print("New undocumented public files:")
    for new_api in new_apis_to_document:
        print("   "+str(new_api))
    print("Note: you will need to update the src_doc_mapping dictionary in check_api_updates.py")
    print("")
    RETURN_VALUE = 1

if undocumented_changes:
    print("Likely undocumented public files:")
    for undoc_change in undocumented_changes:
        print("   "+str(undoc_change))
    RETURN_VALUE = 1

sys.exit(RETURN_VALUE)
