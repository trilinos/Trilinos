#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n evaluateVectorGradField ..... MEX interface for function evaluateVectorGradField.\n\n"
                    "\t evaluateVectorGradField(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (3D array of size [#spaceDim x #numFields x #cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (3D array of size [#cubPoints x #fields x #cells])\n");

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
    // Get the dimensions of the input vector data values array
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
    int invalRank = inputFields.rank();
    int outvalRank = outputFields.rank();
    int numCells = inputFields.dimension(0);
    int numNodes = inputFields.dimension(1);
    int numQPs = inputFields.dimension(2);
    int numDims = inputFields.dimension(3);

    // This is needed, since evaluate currently sums into
    for(int cell = 0; cell < numCells; ++cell)
    {
        for(int node = 0; node < numNodes; ++node)
        {
            for(int qp = 0; qp < numQPs; ++qp)
            {
                for(int i = 0; i < numDims; i++)
                {
                    for(int dim = 0; dim < numDims; dim++)
                    {
                        outputFields(cell, qp, i, dim) += inputData(cell, node, i) * inputFields(cell, node, qp, dim);
                    }
                }
            }
        }
    }

    // Intrepid::FunctionSpaceTools::evaluate<ScalarT>(grad_val_qp, val_node, GradBF);

}
