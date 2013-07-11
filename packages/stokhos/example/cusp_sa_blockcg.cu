#include <cusp/precond/block_smoothed_aggregation.h>
#include <cusp/gallery/poisson.h>
#include <cusp/csr_matrix.h>
#include <cusp/detail/timer.h>
#include <iostream>
#include <cusp/krylov/blockcg.h>
#include <cusp/block_monitor.h>

int main(void)
{

    cudaSetDevice(1);
    cudaDeviceReset();

    typedef int                 IndexType;
    typedef float               ValueType;
    typedef cusp::device_memory MemorySpace;

    // create an empty sparse matrix structure
    cusp::csr_matrix<IndexType, ValueType, MemorySpace> A;

    // create 2D Poisson problem
    cusp::gallery::poisson5pt(A, 128, 256);
    std::cout << "N =" << A.num_rows<< std::endl;
    std::cout << "nnz of A = " << A.num_entries << std::endl;

    // solve with multiple RHS 
    {
        std::cout << "\nSolving Ax = b with multiple RHS..." << std::endl;
        // allocate storage for solution (x) and right hand side (b)
        cusp::array2d<ValueType, MemorySpace> x(A.num_rows, 128, 0);
        cusp::array2d<ValueType, MemorySpace> b(A.num_rows, x.num_cols, 1);
	IndexType numRHS = x.num_cols;

	std::cout << "numRHS = " << x.num_cols << std::endl;

	// set stopping criteria (iteration_limit = 1000, absolute_tolerance = 1e-6)
        cusp::default_block_monitor<ValueType> monitor(b, 100, 1e-6, 0);
        cusp::detail::timer t2;
	// setup preconditioner
        cusp::precond::aggregation::block_smoothed_aggregation<IndexType, ValueType, MemorySpace> M(A, numRHS);
        // solve
	t2.start();	
        cusp::krylov::blockcg(A, x, b, monitor, M);
        std::cout << "Solve Time = " << t2.seconds_elapsed() << std::endl;
        // print hierarchy information
        std::cout << "\nPreconditioner statistics" << std::endl;
        M.print();
    }

    return 0;
}

