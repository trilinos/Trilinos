#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_evaluate ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::evaluate.\n\n"
                    "\tintrepid_evaluate(physValues,inCoeffs,inFields)\n\n"
                    "\t<1-in/out> physValues = Values in physical space (4D array of size [spaceDim x spaceDim x #evalPoints x #cells] or 3D array of size [spaceDim x #evalPoints x #cells] or 2D array size [#evalPoints x #cells])\n"
                    "\t<2-in>     basisCoeffs = Coefficients associated with the fields (basis) array (2D array of size [#fields x #cells])\n"
                    "\t<3-in>     basisFields = Field (basis) values (5D array of size [spaceDim x spaceDim x #evalPoints x #fields x #cells] or 4D array of size [spaceDim x #evalPoints x #fields x #cells] or 3D array size [#evalPoints x #fields x #cells])\n\n");

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

    // Get the dimensions of the point values array
    Teuchos::Array<int> phys_vals_dims;
    m2iGetArrayDims(phys_vals_dims, pInput[0]);
    // Get the dimensions of the coefficients associated with the fields array
    Teuchos::Array<int> coeffs_dims;
    m2iGetArrayDims(coeffs_dims, pInput[1]);
    // Get the dimensions of the field (basis) values array
    Teuchos::Array<int> fields_dims;
    m2iGetArrayDims(fields_dims, pInput[2], true);

    // Get the (pointers to) data
    double* phys_vals_raw = mxGetPr(pInput[0]);
    double* coeffs_raw = mxGetPr(pInput[1]);
    double* fields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> outPointVals(phys_vals_dims, phys_vals_raw);
    FieldContainer<double> inCoeffs(coeffs_dims, coeffs_raw);
    FieldContainer<double> inFields(fields_dims, fields_raw);

    try
    {
        FunctionSpaceTools::evaluate<double>(outPointVals, inCoeffs, inFields);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_setJacobian.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
