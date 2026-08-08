#! /usr/bin/env python3

# SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

"""Check if an API change in Kokkos Kernels might be undocumented

This script is intended to run in a github action when a Pull Request
is opened against Kokkos Kernels' develop branch. It will read the a
pipe or redirected list of files to get a list of modified files. If
one of the modified files is a key in SRC_DOC_MAPPING it checks if one
of the corresponding documentation files has also been modified which
indicates some amount of documentation was done. Otherwise the script
will return with exit core 1. For API files stored in one of the
'SRC_DIRECTORIES' but not listed in 'SRC_DOC_MAPPING' an exit code 1
is also returned.
"""

import sys, argparse, io, difflib, pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

SRC_DIRECTORIES = [
    "batched/dense/src",
    "batched/sparse/src",
    "blas/src",
    "common/src",
    "graph/src",
    "lapack/src",
    "ode/src",
    "sparse/src"
]

SRC_DOC_MAPPING = dict([
    ('batched/dense/src/KokkosBatched_ApplyHouseholder_Decl.hpp', ['docs/source/API/batched/dense/batched_apply_householder.rst']),
    ('batched/dense/src/KokkosBatched_Householder_Decl.hpp', ['docs/source/API/batched/dense/batched_householder.rst']),
    ('blas/src/KokkosBlas1_abs.hpp', ['docs/source/API/blas/blas1_abs.rst']),
    ('blas/src/KokkosBlas1_axpby.hpp', ['docs/source/API/blas/blas1_axpy.rst', 'docs/source/API/blas/blas1_axpby.rst']),
    ('blas/src/KokkosBlas1_dot.hpp', ['docs/source/API/blas/blas1_dot.rst']),
    # ('blas/src/KokkosBlas1_fill.hpp', ),
    ('blas/src/KokkosBlas1_iamax.hpp', ['docs/source/API/blas/blas1_iamax.rst']),
    # ('blas/src/KokkosBlas1_mult.hpp', ),
    ('blas/src/KokkosBlas1_nrm1.hpp', ['docs/source/API/blas/blas1_nrm1.rst']),
    ('blas/src/KokkosBlas1_nrm2.hpp', ['docs/source/API/blas/blas1_nrm2.rst']),
    # ('blas/src/KokkosBlas1_nrm2_squared.hpp', ),
    # ('blas/src/KokkosBlas1_nrm2w.hpp', ),
    # ('blas/src/KokkosBlas1_nrm2w_squared.hpp', ),
    # ('blas/src/KokkosBlas1_nrminf.hpp', ['']),
    # ('blas/src/KokkosBlas1_reciprocal.hpp', ['']),
    ('blas/src/KokkosBlas1_rot.hpp', ['docs/source/API/blas/blas1_rot.rst']),
    ('blas/src/KokkosBlas1_rotg.hpp', ['docs/source/API/blas/blas1_rotg.rst']),
    ('blas/src/KokkosBlas1_rotm.hpp', ['docs/source/API/blas/blas1_rotm.rst']),
    ('blas/src/KokkosBlas1_rotmg.hpp', ['docs/source/API/blas/blas1_rotmg.rst']),
    ('blas/src/KokkosBlas1_scal.hpp', ['docs/source/API/blas/blas1_scal.rst']),
    ('blas/src/KokkosBlas1_set.hpp', ['docs/source/API/blas/blas1_set.rst']),
    # ('blas/src/KokkosBlas1_sum.hpp', ['']),
    ('blas/src/KokkosBlas1_swap.hpp', ['docs/source/API/blas/blas1_swap.rst']),
    # ('blas/src/KokkosBlas1_update.hpp', ['']),
    ('blas/src/KokkosBlas2_gemv.hpp', ['docs/source/API/blas/blas2_gemv.rst']),
    ('blas/src/KokkosBlas2_ger.hpp', ['docs/source/API/blas/blas2_ger.rst']),
    ('blas/src/KokkosBlas2_syr.hpp', ['docs/source/API/blas/blas2_syr.rst']),
    ('blas/src/KokkosBlas2_syr2.hpp', ['docs/source/API/blas/blas2_syr2.rst']),
    ('blas/src/KokkosBlas3_gemm.hpp', ['docs/source/API/blas/blas3_gemm.rst']),
    ('blas/src/KokkosBlas3_trmm.hpp', ['docs/source/API/blas/blas3_trmm.rst']),
    ('blas/src/KokkosBlas3_trsm.hpp', ['docs/source/API/blas/blas3_trsm.rst']),
    ('common/src/KokkosKernels_LowerBound.hpp', ['docs/source/API/common/lower_bound.rst']),
    ('common/src/KokkosKernels_UpperBound.hpp', ['docs/source/API/common/upper_bound.rst']),
    ('lapack/src/KokkosLapack_geqrf.hpp', ['docs/source/API/lapack/geqrf.rst']),
    ('lapack/src/KokkosLapack_gemqr.hpp', ['docs/source/API/lapack/gemqr.rst']),
    ('lapack/src/KokkosLapack_gegqr.hpp', ['docs/source/API/lapack/gegqr.rst']),
    ('lapack/src/KokkosLapack_potrf.hpp', ['docs/source/API/lapack/potrf.rst']),
    ('lapack/src/KokkosLapack_potrs.hpp', ['docs/source/API/lapack/potrs.rst']),
    ('lapack/src/KokkosLapack_gesv.hpp', ['docs/source/API/lapack/gesv.rst']),
    ('lapack/src/KokkosLapack_svd.hpp', ['docs/source/API/lapack/gesvd.rst']),
    ('lapack/src/KokkosLapack_trtri.hpp', ['docs/source/API/lapack/trtri.rst']),
    ('graph/src/KokkosGraph_Distance1Color.hpp', ['docs/source/API/graph/distance1_color.rst']),
    ('graph/src/KokkosGraph_Distance2Color.hpp', ['docs/source/API/graph/distance2_color.rst']),
    ('graph/src/KokkosGraph_LoadBalance.hpp', ['docs/source/API/graph/load_balance.rst']),
    ('graph/src/KokkosGraph_Merge.hpp', ['docs/source/API/graph/merge.rst']),
    ('graph/src/KokkosGraph_MergePath.hpp', ['docs/source/API/graph/merge_path.rst']),
    ('graph/src/KokkosGraph_RCB.hpp', ['docs/source/API/graph/rcb.rst']),
    ('sparse/src/KokkosKernels_Handle.hpp', [
        'docs/source/API/sparse/kokkoskernelshandle.rst',
        'docs/source/API/sparse/handle_get_create_destroy.rst']),
    ('sparse/src/KokkosSparse_BsrMatrix.hpp', [
        'docs/source/API/sparse/bsr_matrix.rst',
        'docs/source/API/sparse/BsrMatrix_constructors.rst',
        'docs/source/API/sparse/bsr_row_view.rst']),
    ('sparse/src/KokkosSparse_CrsMatrix.hpp', [
        'docs/source/API/sparse/crs_matrix.rst',
        'docs/source/API/sparse/CrsMatrix_constructors.rst',
        'docs/source/API/sparse/sparse_row_view.rst']),
    ('sparse/src/KokkosSparse_CcsMatrix.hpp', [
        'docs/source/API/sparse/ccs_matrix.rst',
        'docs/source/API/sparse/CcsMatrix_constructors.rst']),
    ('sparse/src/KokkosSparse_CooMatrix.hpp', [
        'docs/source/API/sparse/coo_matrix.rst',
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
    ('sparse/src/KokkosSparse_Utils.hpp', [
        'docs/source/API/sparse/extract_diagonal_blocks_rcb_deprecated.rst',
        'docs/source/API/sparse/extract_diagonal_blocks_rcb.rst']),
    # ('sparse/src/KokkosSparse_Utils_cusparse.hpp', []),
    # ('sparse/src/KokkosSparse_Utils_mkl.hpp', []),
    # ('sparse/src/KokkosSparse_Utils_rocsparse.hpp', []),
    # ('sparse/src/KokkosSparse_ccs2crs.hpp', []),
    # ('sparse/src/KokkosSparse_coo2crs.hpp', []),
    # ('sparse/src/KokkosSparse_crs2ccs.hpp', []),
    # ('sparse/src/KokkosSparse_crs2coo.hpp', []),
    # ('sparse/src/KokkosSparse_findRelOffset.hpp', []),
    ('sparse/src/KokkosSparse_gauss_seidel.hpp', [
        'docs/source/API/sparse/gauss_seidel_apply_backward.rst',
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
    ('sparse/src/KokkosSparse_spadd.hpp', [
        'docs/source/API/sparse/spadd_numeric.rst',
        'docs/source/API/sparse/spadd_symbolic.rst']),
    # ('sparse/src/KokkosSparse_spadd_handle.hpp', [])
    ('sparse/src/KokkosSparse_spgemm.hpp', [
        'docs/source/API/sparse/spgemm_numeric.rst',
        'docs/source/API/sparse/spgemm_symbolic.rst']),
    # ('sparse/src/KokkosSparse_spgemm_handle.hpp', [])
    # ('sparse/src/KokkosSparse_spgemm_jacobi.hpp', [])
    ('sparse/src/KokkosSparse_spgemm_numeric.hpp', ['docs/source/API/sparse/spgemm_numeric.rst']),
    ('sparse/src/KokkosSparse_spgemm_symbolic.hpp', ['docs/source/API/sparse/spgemm_symbolic.rst']),
    ('sparse/src/KokkosSparse_spiluk.hpp', [
        'docs/source/API/sparse/spiluk_numeric.rst',
        'docs/source/API/sparse/spiluk_symbolic.rst']),
    # ('sparse/src/KokkosSparse_spiluk_handle.hpp', []),
    ('sparse/src/KokkosSparse_sptrsv.hpp', [
        'docs/source/API/sparse/sptrsv_solve.rst',
        'docs/source/API/sparse/sptrsv_symbolic.rst']),
    ('sparse/src/KokkosSparse_trsv.hpp', [
        'docs/source/API/sparse/sptrsv_solve.rst',
        'docs/source/API/sparse/sptrsv_symbolic.rst']),
    # ('sparse/src/KokkosSparse_sptrsv_cholmod.hpp', ['']),
    # ('sparse/src/KokkosSparse_sptrsv_handle.hpp', [])
    # ('sparse/src/KokkosSparse_sptrsv_superlu.hpp', [])
    # ('sparse/src/KokkosSparse_sptrsv_supernode.hpp', [])
    ('sparse/src/KokkosSparse_spmv.hpp', ['docs/source/API/sparse/spmv.rst'])
    # ('sparse/src/KokkosSparse_spmv_deprecated.hpp', []),
    # ('sparse/src/KokkosSparse_spmv_handle.hpp', [])
])

###############################################################################
def check_api_updates_impl(verbose=False, input_files=None):
###############################################################################
    modified_files = set(input_files)

    modified_public_files = []
    new_apis_to_document  = []
    undocumented_changes  = []
    all_known_doc_files   = {x for val in SRC_DOC_MAPPING.values() for x in val}
    indent = "    "
    fail_prefix = "FAILED CHECK: "
    pass_prefix = "PASSED CHECK: "

    if verbose:
        print("All modified files")
        for modified_file in sorted(modified_files):
            print(f"{indent}{modified_file}")
        print()

        modified_doc_files = modified_files & all_known_doc_files
        print("All modified documentation files")
        for modified_doc_file in sorted(modified_doc_files):
            print(f"{indent}{modified_doc_file}")
        print()

        print("Running checks ...")

    # Loop over changed files
    for modified_file in sorted(modified_files):
        if verbose:
            print(f"{indent}Checking modified file {modified_file}")

        # Check if file belongs to one of our source directories
        if any(src_dir in modified_file for src_dir in SRC_DIRECTORIES):
            modified_public_files.append(modified_file)
            # Look for documentation associated with modified file
            if modified_file not in SRC_DOC_MAPPING:
                new_apis_to_document.append(modified_file)
                if verbose:
                    print(f"{indent*2}{fail_prefix}This file has no src doc mapping!")
            else:
                # Construct the intersection of modified files
                # and associated documentation files that should
                # be modified. If the intersection is empty some
                # documentation is likely missing!
                doc_files = set(SRC_DOC_MAPPING[modified_file])
                intersection = doc_files & modified_files
                if not intersection:
                    if verbose:
                        print(f"{indent*2}{fail_prefix}This file has potential undocumented changes due to unchanged {doc_files}!")
                    undocumented_changes.append(modified_file)
                else:
                    if verbose:
                        print(f"{indent*2}{pass_prefix}This file has the expected documentation updates, good job")

        elif modified_file in modified_doc_files:
            if verbose:
                print(f"{indent*2}This file appears to be a documentation file")
        else:
            if verbose:
                print(f"{indent*2}This file does not appear to be a public file, skipping it")

    if verbose:
        print()
        print("Modified public files:")
        for modified_file in sorted(modified_public_files):
            print(f"{indent}{modified_file}")
        print()

    result = True
    if new_apis_to_document:
        print(f"{fail_prefix}New undocumented public files:")
        for new_api in sorted(new_apis_to_document):
            print(f"{indent}{new_api}")
        print(f"{indent}Note: you will need to update the SRC_DOC_MAPPING dictionary in check_api_updates.py")
        print()
        result = False

    if undocumented_changes:
        print(f"{fail_prefix}Likely undocumented public files:")
        for undoc_change in sorted(undocumented_changes):
            print(f"{indent}{undoc_change}")
        print()
        result = False

    return result

###############################################################################
# Test data: covers all three scenarios checked by the script:
#   1. A file whose associated documentation was also modified    -> PASS
#   2. A file in the mapping whose documentation was NOT modified -> FAIL undocumented
#   3. A file in SRC_DIRECTORIES with no mapping entry            -> FAIL new API
#   4. A file outside SRC_DIRECTORIES                             -> ignored
_TEST_MODIFIED_FILES = [
    "blas/src/KokkosBlas1_dot.hpp",               # mapped; doc also present -> pass
    "docs/source/API/blas/blas1_dot.rst",         # doc for dot
    "blas/src/KokkosBlas1_abs.hpp",               # mapped; doc NOT present  -> fail undocumented
    "blas/src/KokkosBlas_NewUndocumentedAPI.hpp", # in SRC_DIRECTORIES, unmapped -> fail new API
    "README.md",                                  # outside SRC_DIRECTORIES -> ignored
]

_TEST_EXPECTED_OUTPUT = \
"""All modified files
    README.md
    blas/src/KokkosBlas1_abs.hpp
    blas/src/KokkosBlas1_dot.hpp
    blas/src/KokkosBlas_NewUndocumentedAPI.hpp
    docs/source/API/blas/blas1_dot.rst

All modified documentation files
    docs/source/API/blas/blas1_dot.rst

Running checks ...
    Checking modified file README.md
        This file does not appear to be a public file, skipping it
    Checking modified file blas/src/KokkosBlas1_abs.hpp
        FAILED CHECK: This file has potential undocumented changes due to unchanged {'docs/source/API/blas/blas1_abs.rst'}!
    Checking modified file blas/src/KokkosBlas1_dot.hpp
        PASSED CHECK: This file has the expected documentation updates, good job
    Checking modified file blas/src/KokkosBlas_NewUndocumentedAPI.hpp
        FAILED CHECK: This file has no src doc mapping!
    Checking modified file docs/source/API/blas/blas1_dot.rst
        This file appears to be a documentation file

Modified public files:
    blas/src/KokkosBlas1_abs.hpp
    blas/src/KokkosBlas1_dot.hpp
    blas/src/KokkosBlas_NewUndocumentedAPI.hpp

FAILED CHECK: New undocumented public files:
    blas/src/KokkosBlas_NewUndocumentedAPI.hpp
    Note: you will need to update the SRC_DOC_MAPPING dictionary in check_api_updates.py

FAILED CHECK: Likely undocumented public files:
    blas/src/KokkosBlas1_abs.hpp

"""

_TEST_EXPECTED_RESULT = False

###############################################################################
def check_api_updates(verbose=False, test=False):
###############################################################################
    if test:
        passed = True

        # Validate every key and value in SRC_DOC_MAPPING points to a real file
        for src_file, doc_files in SRC_DOC_MAPPING.items():
            if not (REPO_ROOT / src_file).is_file():
                print(f"FAILED: SRC_DOC_MAPPING key does not exist: {src_file}")
                passed = False
            for doc_file in doc_files:
                if not (REPO_ROOT / doc_file).is_file():
                    print(f"FAILED: SRC_DOC_MAPPING value does not exist: {doc_file} (mapped from {src_file})")
                    passed = False

        buf = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buf
        try:
            # We always want to test verbose mode
            result = check_api_updates_impl(verbose=True, input_files=_TEST_MODIFIED_FILES)
        finally:
            sys.stdout = old_stdout

        actual_output = buf.getvalue()

        if actual_output.strip() != _TEST_EXPECTED_OUTPUT.strip():
            print("FAILED: Output mismatch")
            diff = difflib.unified_diff(
                _TEST_EXPECTED_OUTPUT.splitlines(keepends=True),
                actual_output.splitlines(keepends=True),
                fromfile="expected",
                tofile="actual",
            )
            print("".join(diff))
            passed = False

        if result != _TEST_EXPECTED_RESULT:
            print(f"FAILED: Expected return value {_TEST_EXPECTED_RESULT}, got {result}")
            passed = False

        if passed:
            print("PASSED: test output matches expected output")

    else:
        return check_api_updates_impl(verbose=True, input_files=[line.strip() for line in sys.stdin])

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        description=description,
        epilog="""
\033[1mEXAMPLES:\033[0m
   \033[1;32m# Check for documentation for changes against origin/develop \033[0m
  > git diff --name-only $(git merge-base origin/develop HEAD) | ./scripts/%(prog)s

   \033[1;32m# Check for documentation in a github action \033[0m
  > git diff --name-only origin/$GITHUB_BASE_REF | ./scripts/%(prog)s
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Produce additional output"
    )

    parser.add_argument(
        "--test",
        action="store_true",
        help="Run built-in self-test with fake data instead of reading from stdin"
    )

    parsed_args = parser.parse_args(args[1:])

    if not parsed_args.test and sys.stdin.isatty():
        raise SystemExit("check_api_updates requires pipe or input redirect")

    return parsed_args

###############################################################################
def main(description):
###############################################################################
    parsed_args = parse_command_line(sys.argv, description)
    return check_api_updates(**vars(parse_command_line(sys.argv, description)))

###############################################################################
if __name__ == "__main__":
    sys.exit(0 if main(__doc__) else 1)
