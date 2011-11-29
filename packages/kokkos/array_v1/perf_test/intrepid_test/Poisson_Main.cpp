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

// Kokkos includes
#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>


#include <Kokkos_DeviceHost_macros.hpp>
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
#include <Kokkos_DeviceClear_macros.hpp>

/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/


namespace Test {


	void poisson_host(int beg , int end)
	{	poisson_run< Kokkos::DeviceHost>(beg , end);  };	

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
