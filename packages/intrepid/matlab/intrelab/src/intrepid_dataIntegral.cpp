#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_dataIntegral ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::dataIntegral.\n\n"
                    "\tintrepid_dataIntegral(outValues,leftValues,rightValues,compEngine,sumInto)\n\n"
                    "\t<1-in/out> outputData = Output data array (1D array of size [#numCells])\n"
                    "\t<2-in>     inputDataLeft = Left data input array (4D array of size [#cubPoints x #cells x spaceDim x spaceDim] or 3D array of size [#cubPoints x #cells x spaceDim] or 2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in>     inputDataRight = Right data input array (4D array of size [#cubPoints x #cells x spaceDim x spaceDim] or 3D array of size [#cubPoints x #cells x spaceDim] or 2D array of size [#cubPoints x #cells])\n"
                    "\t<4-in>     compEngine = 'COMP_CPP' (for-loop computation) | 'COMP_BLAS' (BLAS computation) (string)\n"
                    "\t<5-in>     sumInto = 'true' | 'false' (default) determines if outValues will be summed into or overwritten (string)\n\n");

    // Check the number of input arguments
    if((nInput != 4) && (nInput != 5))
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the output values array
    Teuchos::Array<int> outData_dims;
    m2iGetArrayDims(outData_dims, pInput[0]);
    // Get the dimensions of the left values array
    Teuchos::Array<int> leftInputData_dims;
    m2iGetArrayDims(leftInputData_dims, pInput[1]);
    // Get the dimensions of the right values array
    Teuchos::Array<int> rightInputData_dims;
    m2iGetArrayDims(rightInputData_dims, pInput[2]);

    // Get the (pointers to) data
    double* outData_raw = mxGetPr(pInput[0]);
    double* leftInputData_raw = mxGetPr(pInput[1]);
    double* rightInputData_raw = mxGetPr(pInput[2]);

    const std::string ce_str(mxArrayToString(pInput[3]));
    ECompEngine compEngine = COMP_ENGINE_MAX;
    if(ce_str == "COMP_CPP")
    {
        compEngine = COMP_CPP;
    }
    else if(ce_str == "COMP_BLAS")
    {
        compEngine = COMP_BLAS;
    }

    bool sumInto = false;
    if(nInput == 4)
    {
        sumInto = false;
    }
    else if(nInput == 5)
    {
        const std::string si_str(mxArrayToString(pInput[4]));
        if(si_str == "true")
        {
            sumInto = true;
        }
    }

    FieldContainer<double> outData(outData_dims, outData_raw);
    FieldContainer<double> leftInputData(leftInputData_dims, leftInputData_raw);
    FieldContainer<double> rightInputData(rightInputData_dims, rightInputData_raw);

    try
    {
        FunctionSpaceTools::dataIntegral<double>(outData, leftInputData, rightInputData, compEngine, sumInto);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_dataIntegral.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
