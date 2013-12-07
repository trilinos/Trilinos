#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_HGRADtransformVALUE ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::HGRADtransformVALUE.\n\n"
                    "\tintrepid_HGRADtransformVALUE(physValues,refValues)\n\n"
                    "\t<1-in/out> physValues = Values in physical space (3D array of size [#cubPoints x #fields x #cells])\n"
                    "\t<2-in>     refValues = Values on the reference cell (2D array of size [#cubPoints x #fields])\n");

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

    // Get the dimensions of the physical values array
    Teuchos::Array<int> pVals_dims;
    m2iGetArrayDims(pVals_dims, pInput[0]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> rVals_dims;
    m2iGetArrayDims(rVals_dims, pInput[1]);

    // Get the (pointers to) data
    double* pVals_raw = mxGetPr(pInput[0]);
    double* rVals_raw = mxGetPr(pInput[1]);

    FieldContainer<double> physVals(pVals_dims, pVals_raw);
    FieldContainer<double> refVals(rVals_dims, rVals_raw);

    try
    {
        FunctionSpaceTools::HGRADtransformVALUE<double>(physVals, refVals);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_HGRADtransformVALUE.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
