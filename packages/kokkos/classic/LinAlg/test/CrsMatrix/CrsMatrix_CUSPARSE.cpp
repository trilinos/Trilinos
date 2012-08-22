//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include "Kokkos_DefaultArithmetic.hpp"

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_ThrustGPUNode.hpp>
#include <Kokkos_CUSPARSEOps.hpp>

namespace {

  using std::endl;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Kokkos::ThrustGPUNode;
  using Teuchos::ParameterList;
  using Teuchos::null;

  TEUCHOS_UNIT_TEST( CrsMatrix, CUSPARSENative )
  {
    Kokkos::CUSPARSEdetails::Session::init();
    RCP<const cusparseHandle_t> sess = Kokkos::CUSPARSEdetails::Session::getHandle();
    
    ParameterList pl;
    RCP<ThrustGPUNode> gpunode = rcp(new ThrustGPUNode(pl));

    cusparseStatus_t status ; 

    // problem characteristics
    int n , nnz , nnz_vector ; 
    float dzero =0.0; 
    float dtwo =2.0; 
    float dthree =3.0; 
    float dfive =5.0;
    nnz_vector = 3; 
    n=4; nnz=9; 

    out << "Testing example" << endl;

    // device pointers
    ArrayRCP<int> cooRowIndex, cooColIndex, csrRowPtr;
    ArrayRCP<float> cooVal;
    ArrayRCP<int> xInd;
    ArrayRCP<float> xVal, y, z;
    cooRowIndex = gpunode->allocBuffer<int>(nnz);
    cooColIndex = gpunode->allocBuffer<int>(nnz);
    cooVal      = gpunode->allocBuffer<float>(nnz);
    y           = gpunode->allocBuffer<float>(2*n);
    xInd        = gpunode->allocBuffer<int>(nnz_vector);
    xVal        = gpunode->allocBuffer<float>(nnz_vector);
    csrRowPtr   = gpunode->allocBuffer<int>(n+1);
    z           = gpunode->allocBuffer<float>(2*(n+1));

    // init data on host
    {
      // host pointers
      ArrayRCP<int>   cooRowIndexHostPtr, cooColIndexHostPtr;
      ArrayRCP<float> cooValHostPtr;
      ArrayRCP<int> xIndHostPtr;
      ArrayRCP<float> xValHostPtr, yHostPtr;

      cooRowIndexHostPtr = gpunode->viewBufferNonConst(Kokkos::WriteOnly, cooRowIndex.size(), cooRowIndex);
      cooColIndexHostPtr = gpunode->viewBufferNonConst(Kokkos::WriteOnly, cooColIndex.size(), cooColIndex);
      cooValHostPtr      = gpunode->viewBufferNonConst(Kokkos::WriteOnly, cooVal.size(),      cooVal);
      yHostPtr           = gpunode->viewBufferNonConst(Kokkos::WriteOnly, y.size(),           y);
      xIndHostPtr        = gpunode->viewBufferNonConst(Kokkos::WriteOnly, xInd.size(),        xInd);
      xValHostPtr        = gpunode->viewBufferNonConst(Kokkos::WriteOnly, xVal.size(),        xVal);

      cooRowIndexHostPtr[0] = 0; cooColIndexHostPtr[0] = 0; cooValHostPtr[0] = 1.0; 
      cooRowIndexHostPtr[1] = 0; cooColIndexHostPtr[1] = 2; cooValHostPtr[1] = 2.0; 
      cooRowIndexHostPtr[2] = 0; cooColIndexHostPtr[2] = 3; cooValHostPtr[2] = 3.0; 
      cooRowIndexHostPtr[3] = 1; cooColIndexHostPtr[3] = 1; cooValHostPtr[3] = 4.0; 
      cooRowIndexHostPtr[4] = 2; cooColIndexHostPtr[4] = 0; cooValHostPtr[4] = 5.0; 
      cooRowIndexHostPtr[5] = 2; cooColIndexHostPtr[5] = 2; cooValHostPtr[5] = 6.0; 
      cooRowIndexHostPtr[6] = 2; cooColIndexHostPtr[6] = 3; cooValHostPtr[6] = 7.0; 
      cooRowIndexHostPtr[7] = 3; cooColIndexHostPtr[7] = 1; cooValHostPtr[7] = 8.0; 
      cooRowIndexHostPtr[8] = 3; cooColIndexHostPtr[8] = 3; cooValHostPtr[8] = 9.0; 
      yHostPtr[0] = 10.0; xIndHostPtr[0]=0; xValHostPtr[0]=100.0; 
      yHostPtr[1] = 20.0; xIndHostPtr[1]=1; xValHostPtr[1]=200.0; 
      yHostPtr[2] = 30.0; 
      yHostPtr[3] = 40.0; xIndHostPtr[2]=3; xValHostPtr[2]=400.0; 
      yHostPtr[4] = 50.0; 
      yHostPtr[5] = 60.0; 
      yHostPtr[6] = 70.0; 
      yHostPtr[7] = 80.0;

      cooRowIndexHostPtr = null;
      cooColIndexHostPtr = null;
      cooValHostPtr = null;
      yHostPtr = null;
      xIndHostPtr = null;
      xValHostPtr = null;
    }

    /* create and setup matrix descriptor */
    Teuchos::RCP<cusparseMatDescr_t> descr = Kokkos::CUSPARSEdetails::createMatDescr();
    cusparseSetMatType(      *descr , CUSPARSE_MATRIX_TYPE_GENERAL ) ; 
    cusparseSetMatIndexBase( *descr , CUSPARSE_INDEX_BASE_ZERO ) ;
    status = cusparseXcoo2csr(*sess,cooRowIndex.getRawPtr(),nnz,n, csrRowPtr.getRawPtr() , CUSPARSE_INDEX_BASE_ZERO ) ; 
    if ( status != CUSPARSE_STATUS_SUCCESS ) {
       success = false; return;
    }
    // scatter test
    status = cusparseSsctr(*sess , nnz_vector , xVal.getRawPtr() , xInd.getRawPtr() , y.getRawPtr()+n , CUSPARSE_INDEX_BASE_ZERO ) ;
    if ( status != CUSPARSE_STATUS_SUCCESS ) { 
      success = false; return;
    }
    // sparse matvec
    status = cusparseScsrmv ( *sess , CUSPARSE_OPERATION_NON_TRANSPOSE , n , n , nnz , &dtwo , *descr , cooVal.getRawPtr() , csrRowPtr.getRawPtr() , cooColIndex.getRawPtr() , y.getRawPtr() , &dthree , y.getRawPtr()+n) ; 
    if ( status != CUSPARSE_STATUS_SUCCESS ) {
       success = false; return;
    }
    // matvec result
    ArrayRCP<const float> yHostPtr, zHostPtr;
    yHostPtr = gpunode->viewBuffer<float>(y.size(), y);
    // mat-multivec
    cudaError_t cudaStat1 = cudaMemset((void *)z.getRawPtr(),0, 2*(n+1)*sizeof(z[0])); 
    if ( cudaStat1 != cudaSuccess ) {
       success = false; return;
    } 
    status = cusparseScsrmm( *sess , CUSPARSE_OPERATION_NON_TRANSPOSE , n , 2 , n , nnz , &dfive , *descr , cooVal.getRawPtr() , csrRowPtr.getRawPtr() , cooColIndex.getRawPtr() , y.getRawPtr(), n, &dzero, z.getRawPtr(), n+1); 
    if ( status != CUSPARSE_STATUS_SUCCESS ) {
       success = false; return;
    }
    zHostPtr = gpunode->viewBuffer<float>( z.size(), z );
    if ( status != CUSPARSE_STATUS_SUCCESS ) {
       success = false; return;
    }
    if ( status != CUSPARSE_STATUS_SUCCESS ) { 
       success = false; return;
    }
    if (    (zHostPtr[0] != 950.0)    || (zHostPtr[1] != 400.0)   || (zHostPtr[2] != 2550.0)  || (zHostPtr[3] != 2600.0) 
         || (zHostPtr[4] != 0.0)      || (zHostPtr[5] != 49300.0) || (zHostPtr[6] != 15200.0) || (zHostPtr[7] != 132300.0) 
         || (zHostPtr[8] != 131200.0) || (zHostPtr[9] != 0.0)     || (yHostPtr[0] != 10.0)    || (yHostPtr[1] != 20.0) 
         || (yHostPtr[2] != 30.0)     || (yHostPtr[3] != 40.0)    || (yHostPtr[4] != 680.0)   || (yHostPtr[5] != 760.0) 
         || (yHostPtr[6] != 1230.0)   || (yHostPtr[7] != 2240.0))
    {
       success = false;
    } else {
       success = true;
    }
  }

}
