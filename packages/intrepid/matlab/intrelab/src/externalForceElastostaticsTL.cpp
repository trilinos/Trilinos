#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n internalForceElastostaticsTL ..... MEX interface for function internalForceElastostaticsTL.\n\n"
                    "\t internalForceElastostaticsTL(internalForces,stress,inputFields)\n\n"
                    "\t<1-in/out> internalForces = Output fields array (3D array of size [#spaceDim x #numFields x #cells])\n"
                    "\t<2-in> 	   force = Input data array (3D array of size [#spaceDim x #cubPoints x #cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (4D array of size [#spaceDim x #cubPoints x #fields x #cells])\n");

    // Check the number of input arguments
    if(nInput != 3)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the output data values array
    Teuchos::Array<int> oData_dims;
    m2iGetArrayDims(oData_dims, pInput[0]);
    // Get the dimensions of the input data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> externalForce(oData_dims, oData_raw);
    FieldContainer<double> force(iData_dims, iData_raw);
    FieldContainer<double> wBF(iFields_dims, iFields_raw);

    // get sizes
    int numCells = wBF.dimension(0);
    int numNodes = wBF.dimension(1);
    int numQPs = wBF.dimension(2);
    int numDims = force.dimension(2);

    // compute external forces
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int node = 0; node < numNodes; ++node)
        {
            for(int qp = 0; qp < numQPs; ++qp)
            {
                for(int i = 0; i < numDims; i++)
                {
                    externalForce(cell, node, i) += force(cell, qp, i) * wBF(cell, node, qp);
                }
            }
        }
    }

}
