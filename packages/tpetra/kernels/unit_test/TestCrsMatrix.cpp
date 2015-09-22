//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <stdexcept>

namespace { // anonymous

  using std::cerr;
  using std::endl;

  // Create a test sparse matrix A.
  //
  // Identify the matrix to create by number (whichMatrix).  The
  // following lists the valid options for whichMatrix:
  //
  // 0: A square 10 x 10 nonsymmetric diagonally dominant sparse matrix.
  //
  // \param ptr [out] Array of row offsets, of length numRows+1.
  // \param ind [out] Array of column indices, of length nnz.
  // \param val [out] Array of entries (values), of length nnz.
  // \param numRows [out] The number of rows in the matrix.
  // \param numCols [out] The number of columns in the matrix.
  // \param nnz [out] The number of stored entries in the matrix.
  // \param whichMatrix [in] The index of the matrix to create.
  template<typename MemorySpace>
  void
  makeSparseMatrix (Kokkos::View<typename MemorySpace::size_type*, MemorySpace>& ptr,
                    Kokkos::View<int*, MemorySpace>& ind,
                    Kokkos::View<double*, MemorySpace>& val,
                    int& numRows,
                    int& numCols,
                    int& nnz,
                    const int whichMatrix)
  {
    using Kokkos::HostSpace;
    using Kokkos::MemoryUnmanaged;
    using Kokkos::View;
    typedef typename MemorySpace::size_type size_type;

    if (whichMatrix == 0) {
      numRows = 10;
      numCols = 10;
      nnz = 21;
      const size_type ptrRaw[] = {0, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21};
      const int indRaw[] = {0, 1, 9,
                            1, 2,
                            2, 3,
                            3, 4,
                            4, 5,
                            5, 6,
                            6, 7,
                            7, 8,
                            8, 9,
                            1, 9};
      const double valRaw[] = {1.0, 4.0, 0.5,
                               0.5, 5.0,
                               1.0, 6.0,
                               1.5, 7.0,
                               2.0, 8.0,
                               2.5, 9.0,
                               3.0, 10.0,
                               3.5, 11.0,
                               4.0, 12.0,
                               4.5, 13.0};

      typedef View<size_type*,MemorySpace> ptr_type ;
      typedef View<int*,   MemorySpace> ind_type ;
      typedef View<double*,MemorySpace> val_type ;

      // Create the output Views.
      ptr = ptr_type("ptr", numRows + 1);
      ind = ind_type("ind", nnz);
      val = val_type("val", nnz);

      // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
      typename ptr_type::HostMirror::const_type  ptrIn( ptrRaw , numRows+1 );
      typename ind_type::HostMirror::const_type  indIn( indRaw , nnz );
      typename val_type::HostMirror::const_type  valIn( valRaw , nnz );

      Kokkos::deep_copy (ptr, ptrIn);
      Kokkos::deep_copy (ind, indIn);
      Kokkos::deep_copy (val, valIn);
    }
    else { // whichMatrix != 0
      std::ostringstream os;
      os << "Invalid whichMatrix value " << whichMatrix
         << ".  Valid value(s) include " << 0 << ".";
      throw std::invalid_argument (os.str ());
    }
  }

  // Return the Kokkos::CrsMatrix corresponding to makeSparseMatrix().
  template<typename MemorySpace>
  Kokkos::CrsMatrix<double, int, MemorySpace>
  makeCrsMatrix ()
  {
    Kokkos::View<typename MemorySpace::size_type*, MemorySpace> ptr;
    Kokkos::View<int*, MemorySpace> ind;
    Kokkos::View<double*, MemorySpace> val;
    int numRows;
    int numCols;
    int nnz;

    const int whichMatrix = 0;
    makeSparseMatrix<MemorySpace> (ptr, ind, val, numRows, numCols, nnz, whichMatrix);
    typedef Kokkos::CrsMatrix<double, int, MemorySpace> crs_matrix_type;
    return crs_matrix_type ("A", numRows, numCols, nnz, val, ptr, ind);
  }

  // Create a Kokkos::CrsMatrix.  This mainly tests that the class
  // compiles.  However, it does need to initialize the MemorySpace's
  // default execution space, because it allocates Views and calls
  // deep_copy a few times.
  template<typename MemorySpace>
  void
  testCrsMatrix ()
  {
    Kokkos::initialize();

    typedef Kokkos::CrsMatrix<double, int, MemorySpace> crs_matrix_type;
    crs_matrix_type A = makeCrsMatrix<MemorySpace> ();
    // mfh 28 Sep 2013: Use A in some way, so the compiler can't
    // optimize it away completely.  This forces the compiler to
    // compile CrsMatrix, which is the whole point of this test.
    printf ("A is %d by %d\n", A.numRows (), A.numCols ());

    Kokkos::finalize ();
  }

  TEUCHOS_UNIT_TEST( CrsMatrix, Compile )
  {
    // For now, just test that CrsMatrix compiles.
    testCrsMatrix<Kokkos::DefaultExecutionSpace> ();
  }

} // namespace (anonymous)

