#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n tangentStiffnessElastostaticsTL ..... MEX interface for function tangentStiffnessElastostaticsTL.\n\n"
                    "\t tangentStiffnessElastostaticsTL(internalForces,stress,inputFields)\n\n"
                    "\t<1-in/out> tangentStiffness = Output data array (3D array of size [#numFields x #numFields x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (4D array of size [#rows x #numFields x #cubPoints x #cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (4D array of size [#rows x #numFields x #cubPoints x #cells])\n");

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
    Teuchos::Array<int> iDataBF_dims;
    m2iGetArrayDims(iDataBF_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iDataBF_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> tangentStiffness(oData_dims, oData_raw);
    FieldContainer<double> inputDataBF(iDataBF_dims, iDataBF_raw);
    FieldContainer<double> wGradBF(iFields_dims, iFields_raw);

    // get sizes
    int numCells = wGradBF.dimension(0);
    int numQPs = wGradBF.dimension(1);
    int numDofs = wGradBF.dimension(2);
    int matDim = wGradBF.dimension(3);

    // compute tangent stiffness matrix
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int qp = 0; qp < numQPs; ++qp)
        {
            for(int i = 0; i < numDofs; ++i)
            {
                for(int j = 0; j < numDofs; j++)
                {
                    for(int k = 0; k < matDim; k++)
                    {
                        tangentStiffness(cell, i, j) += wGradBF(cell, qp, i, k) * inputDataBF(cell, qp, j, k);
                    }
                }
            }
        }
    }

}
