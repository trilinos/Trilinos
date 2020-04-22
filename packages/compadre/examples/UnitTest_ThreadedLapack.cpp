#include <Compadre_Config.h>
#include "Compadre_LinearAlgebra_Definitions.hpp"

#include <assert.h>
#include <Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>

#ifdef COMPADRE_USE_MPI
#include <mpi.h>
#endif

using namespace Compadre;

//! This tests whether or not the LAPACK+BLAS combination is safe to have a parallel_for wrapping them
//! If they are not threadsafe, then nothing in this toolkit will work and another library will have
//! be provided that is threadsafe.
//!
//! Creates num_matrices copies of a matrix P and RHS
//! P has dimensions P_rows x P_cols 
//! RHS has dimensions RHS_rows x RHS_cols
//!
//! P*X=RHS is an overdetermined system from the perspective of rows to columns
//! But it is underdetermined in the sense that it actually has a large nullspace
//! because we fill P with ones, and RHS with ones. 
//!
//! The result is that solution should be 1 / P_cols for each entry in RHS
//!
//! This is a fairly easy exercise for dgelsd in LAPACK, but if LAPACK+BLAS are not
//! thread safe, then this will fail to get the correct answer.
//!
//! If this test fails, then either a different LAPACK library should be provided or
//! compiled. If you are compiling, be sure to add the -frecursive 
//! see https://github.com/xianyi/OpenBLAS/issues/477

// called from command line
int main (int argc, char* args[]) {

// initializes MPI (if available) with command line arguments given
#ifdef COMPADRE_USE_MPI
MPI_Init(&argc, &args);
#endif

int number_wrong = 0;

// initializes Kokkos with command line arguments given
Kokkos::initialize(argc, args);
{
    const int P_rows   = 100;
    const int P_cols   = 50;
    assert((P_rows >= P_cols) && "P must not be underdetermined.");

    const int RHS_rows = P_rows;

    const int num_matrices = 20;

    Kokkos::Profiling::pushRegion("Instantiating Data");
    auto all_P   = Kokkos::View<double*>("P", num_matrices*P_cols*P_rows);
    auto all_RHS = Kokkos::View<double*>("RHS", num_matrices*RHS_rows*RHS_rows);
    Kokkos::Profiling::popRegion();

    Kokkos::parallel_for("Fill Matrices", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,num_matrices), KOKKOS_LAMBDA(const int i) {
        Kokkos::View<double**, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
            P(all_P.data() + i*P_cols*P_rows, P_rows, P_cols);

        Kokkos::View<double**, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
            RHS(all_RHS.data() + i*RHS_rows*RHS_rows, RHS_rows, RHS_rows);

        for (int j=0; j<P_rows; ++j) {
            for (int k=0; k<P_cols; ++k) {
                P(j,k) = 1.0;
            }
        }

        for (int j=0; j<RHS_rows; ++j) {
            for (int k=0; k<RHS_rows; ++k) {
                RHS(j,k) = 1.0;
            }
        }
    });
    Kokkos::fence();


    // call SVD on all num_matrices
    ParallelManager pm;
    GMLS_LinearAlgebra::batchSVDFactorize(pm, true, all_P.data(), P_rows, P_cols, true, all_RHS.data(), RHS_rows, RHS_rows, P_rows, P_cols, P_rows, num_matrices);


    const double tol = 1e-10;
    Kokkos::parallel_reduce("Check Solution", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,num_matrices), KOKKOS_LAMBDA(const int i, int& t_wrong) {
        Kokkos::View<double**, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
            RHS(all_RHS.data() + i*RHS_rows*RHS_rows, RHS_rows, RHS_rows);

        // check the diagonals for the true solution
        for (int j=0; j<P_cols; ++j) {
            if (std::abs(RHS(j,j)-1./P_cols) > tol) {
                t_wrong++;
            }
        }

    }, number_wrong);

}
// finalize Kokkos
Kokkos::finalize();
#ifdef COMPADRE_USE_MPI
MPI_Finalize();
#endif

#ifndef LAPACK_DECLARED_THREADSAFE
    printf("LAPACK_DECLARED_THREADSAFE=OFF. Massive performance loss due to serial LAPACK implementation provided at configure.\n");
    return 77;
#endif

if (number_wrong > 0) {
    printf("Incorrect result. LAPACK IS NOT THREADSAFE AND CANNOT BE USED WITH THIS TOOLKIT! Either provide a thread safe LAPACK+BLAS combination or set -DLAPACK_DECLARED_THREADSAFE:BOOL=OFF in CMake and take a MASSIVE performance hit.\n");
    return -1;
}
return 0;
}
