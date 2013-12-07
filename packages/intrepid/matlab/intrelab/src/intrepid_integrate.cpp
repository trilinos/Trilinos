#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_integrate ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::integrate.\n\n"
                    "\tintrepid_integrate(outValues,leftValues,rightValues,compEngine,sumInto)\n\n"
                    "\t<1-in/out> outValues = Output  (variable-size array, see Intrepid documentation)\n"
                    "\t<2-in>     leftValues = Left inputs (variable-size array, see Intrepid documentation)\n"
                    "\t<3-in>     rightValues = Right inputs (variable-size array, see Intrepid documentation)\n"
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
    Teuchos::Array<int> outVals_dims;
    m2iGetArrayDims(outVals_dims, pInput[0], true);
    // Get the dimensions of the left values array
    Teuchos::Array<int> leftVals_dims;
    m2iGetArrayDims(leftVals_dims, pInput[1]);
    // Get the dimensions of the right values array
    Teuchos::Array<int> rightVals_dims;
    m2iGetArrayDims(rightVals_dims, pInput[2]);

    // Get the (pointers to) data
    double* outVals_raw = mxGetPr(pInput[0]);
    double* leftVals_raw = mxGetPr(pInput[1]);
    double* rightVals_raw = mxGetPr(pInput[2]);

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

    FieldContainer<double> outVals(outVals_dims, outVals_raw);
    FieldContainer<double> leftVals(leftVals_dims, leftVals_raw);
    FieldContainer<double> rightVals(rightVals_dims, rightVals_raw);

    try
    {
        FunctionSpaceTools::integrate<double>(outVals, leftVals, rightVals, compEngine, sumInto);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_integrate.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
