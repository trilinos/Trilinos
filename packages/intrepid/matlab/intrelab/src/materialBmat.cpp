#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n materialBmat ..... MEX interface for function materialBmat.\n\n"
                    "\t materialBmat(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [2*#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in> 	   inputFields = Input fields array (4D array of size [#spaceDim x #cubPoints x #numFields x #cells])\n\n");

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

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    // get sizes
    int numCells = inputFields.dimension(0);
    int numNodes = inputFields.dimension(1);
    int numQPs = inputFields.dimension(2);
    int numDims = inputData.dimension(3);

    switch(numDims)
    {
        case 2:
            for(int cell = 0; cell < numCells; ++cell)
            {
                for(int node = 0; node < numNodes; ++node)
                {
                    for(int qp = 0; qp < numQPs; ++qp)
                    {
                        // row 1
                        outputFields(cell, numDims * node, qp, 0) = inputData(cell, qp, 0, 0)
                                * inputFields(cell, node, qp, 0);
                        outputFields(cell, numDims * node + 1, qp, 0) = inputData(cell, qp, 1, 0)
                                * inputFields(cell, node, qp, 0);
                        // row 2
                        outputFields(cell, numDims * node, qp, 1) = inputData(cell, qp, 0, 1)
                                * inputFields(cell, node, qp, 1);
                        outputFields(cell, numDims * node + 1, qp, 1) = inputData(cell, qp, 1, 1)
                                * inputFields(cell, node, qp, 1);
                        // row 3
                        outputFields(cell, numDims * node, qp, 2) =
                                (inputData(cell, qp, 0, 0) * inputFields(cell, node, qp, 1))
                                + (inputData(cell, qp, 0, 1) * inputFields(cell, node, qp, 0));
                        outputFields(cell, numDims * node + 1, qp, 2) =
                                (inputData(cell, qp, 1, 1) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 1, 0) * inputFields(cell, node, qp, 1));
                    }
                }
            }
            break;
        case 3:
            for(int cell = 0; cell < numCells; ++cell)
            {
                for(int node = 0; node < numNodes; ++node)
                {
                    for(int qp = 0; qp < numQPs; ++qp)
                    {
                        // row 1
                        outputFields(cell, numDims * node, qp, 0) = inputData(cell, qp, 0, 0)
                                * inputFields(cell, node, qp, 0);
                        outputFields(cell, numDims * node + 1, qp, 0) = inputData(cell, qp, 1, 0)
                                * inputFields(cell, node, qp, 0);
                        outputFields(cell, numDims * node + 2, qp, 0) = inputData(cell, qp, 2, 0)
                                * inputFields(cell, node, qp, 0);
                        // row 2
                        outputFields(cell, numDims * node, qp, 1) = inputData(cell, qp, 0, 1)
                                * inputFields(cell, node, qp, 1);
                        outputFields(cell, numDims * node + 1, qp, 1) = inputData(cell, qp, 1, 1)
                                * inputFields(cell, node, qp, 1);
                        outputFields(cell, numDims * node + 2, qp, 1) = inputData(cell, qp, 2, 1)
                                * inputFields(cell, node, qp, 1);
                        // row 3
                        outputFields(cell, numDims * node, qp, 2) = inputData(cell, qp, 0, 2)
                                * inputFields(cell, node, qp, 2);
                        outputFields(cell, numDims * node + 1, qp, 2) = inputData(cell, qp, 1, 2)
                                * inputFields(cell, node, qp, 2);
                        outputFields(cell, numDims * node + 2, qp, 2) = inputData(cell, qp, 2, 2)
                                * inputFields(cell, node, qp, 2);
                        // row 4
                        outputFields(cell, numDims * node, qp, 3) =
                                (inputData(cell, qp, 0, 2) * inputFields(cell, node, qp, 1))
                                + (inputData(cell, qp, 0, 1) * inputFields(cell, node, qp, 2));
                        outputFields(cell, numDims * node + 1, qp, 3) =
                                (inputData(cell, qp, 1, 2) * inputFields(cell, node, qp, 1))
                                + (inputData(cell, qp, 1, 1) * inputFields(cell, node, qp, 2));
                        outputFields(cell, numDims * node + 2, qp, 3) =
                                (inputData(cell, qp, 2, 2) * inputFields(cell, node, qp, 1))
                                + (inputData(cell, qp, 2, 1) * inputFields(cell, node, qp, 2));
                        // row 5
                        outputFields(cell, numDims * node, qp, 4) =
                                (inputData(cell, qp, 0, 2) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 0, 0) * inputFields(cell, node, qp, 2));
                        outputFields(cell, numDims * node + 1, qp, 4) =
                                (inputData(cell, qp, 1, 2) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 1, 0) * inputFields(cell, node, qp, 2));
                        outputFields(cell, numDims * node + 2, qp, 4) =
                                (inputData(cell, qp, 2, 2) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 2, 0) * inputFields(cell, node, qp, 2));
                        // row 6
                        outputFields(cell, numDims * node, qp, 5) =
                                (inputData(cell, qp, 0, 1) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 0, 0) * inputFields(cell, node, qp, 1));
                        outputFields(cell, numDims * node + 1, qp, 5) =
                                (inputData(cell, qp, 1, 1) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 1, 0) * inputFields(cell, node, qp, 1));
                        outputFields(cell, numDims * node + 2, qp, 5) =
                                (inputData(cell, qp, 2, 1) * inputFields(cell, node, qp, 0))
                                + (inputData(cell, qp, 2, 0) * inputFields(cell, node, qp, 1));
                    }
                }
            }
            break;
    }
}
