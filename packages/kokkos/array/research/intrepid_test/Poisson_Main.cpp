/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
*/

/********************************************************************************/
/******** EXAMPLE: Stiffness matrix and right hand side with Intrepid ***********/
/********************************************************************************/
/*
      This example builds element stiffness matrices and RHS 
      contributions for a batch of hexahedral elements using
      first order HGRAD basis functions and 2nd order quadrature.

       Poisson system in weak form:
 
            (grad u, grad v) = (f,v)

       Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

*/
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <fstream>
/**************************************************************/
/*                          Includes                          */
/**************************************************************/

// Intrepid includes
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_CubatureDirectLineGauss.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_CubatureTensor.hpp>


// Teuchos includes
#include "Teuchos_RCP.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

#include "Original.hpp"

// KokkosArray includes
#include <KokkosArray_DeviceHost.hpp>
#include <KokkosArray_DeviceHost_MDArrayView.hpp>
#include <KokkosArray_DeviceHost_MultiVectorView.hpp>
#include <KokkosArray_DeviceHost_ValueView.hpp>
#include <KokkosArray_DeviceHost_ParallelFor.hpp>
#include <KokkosArray_DeviceHost_ParallelReduce.hpp>


#include <KokkosArray_DeviceHost_macros.hpp>
#include <Invert.hpp>
#include <Multiply.hpp>
#include <computeCellMeasure.hpp>
#include <Integrate.hpp>
#include <Determinant.hpp>
#include <Jacobian.hpp>
#include <Transform.hpp>
#include <simpleFill.hpp>
#include <TransformValue.hpp>
#include <Poisson_Driver.hpp>
#include <KokkosArray_DeviceClear_macros.hpp>

/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/


namespace Test {


	void poisson_host(int beg , int end)
	{	poisson_run< KokkosArray::DeviceHost>(beg , end);  };	

	void poisson_cuda(int,int);
	void poisson_tpi(int ,int,int);
	void poisson_tbb(int , int);
}


int main(int argc, char *argv[]) {

/*********************************************************/
/*               Define cell topology                    */
/*********************************************************/
	int end = 20;
	Test::poisson_host(10 , end);
	Test::poisson_cuda(10 , end);
	
	if (argc == 1) {
		Test::poisson_tpi(10,end,1);
	}
	else {
		Test::poisson_tpi(10,end,atoi(argv[1]));
	}
	Test::poisson_tbb(10,end);

	Original(10,end);

	return 0;
}

/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/
