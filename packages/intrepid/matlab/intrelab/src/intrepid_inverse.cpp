#include "mex.h"

#include <Intrepid_FieldContainer.hpp>
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"

#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n intrepid_inverse ..... MEX interface for function intrepid_inverse.\n\n"
                    "\t intrepid_inverse(outputData,inputData)\n\n"
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

    // Get the dimensions of the output data values array
    Teuchos::Array<int> oData_dims;
    m2iGetArrayDims(oData_dims, pInput[0]);
    // Get the dimensions of the input vector data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);

    FieldContainer<double> outputData(oData_dims, oData_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);

    try
    {
        Intrepid::RealSpaceTools<double>::inverse(outputData, inputData);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getCubature.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }

}
