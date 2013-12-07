#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n defGrad ..... MEX interface for function defGrad.\n\n"
                    "\t defGrad(outputFields,inputData)\n\n"
                    "\t<1-in/out> outputData = Output fields array (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (4D array of size [#spaceDim x #spaceDim x #numFields x #cells])\n");

    // Check the number of input arguments
    if(nInput != 2)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input vector data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);

    FieldContainer<double> outputData(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);

    // get sizes
    int invalRank = inputData.rank();
    int outvalRank = outputData.rank();
    int numCells = inputData.dimension(0);
    int numQPs = inputData.dimension(1);
    int numDims = inputData.dimension(2);

    // Compute DefGrad tensor from displacement gradient
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int qp = 0; qp < numQPs; ++qp)
        {
            for(int i = 0; i < numDims; ++i)
            {
                for(int j = 0; j < numDims; ++j)
                {
                    outputData(cell, qp, i, j) = inputData(cell, qp, i, j);
                }
                outputData(cell, qp, i, i) += 1.0;
            }
        }
    }
}
