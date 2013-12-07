#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_scalarMultiplyDataField ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::scalarMultiplyDataField.\n\n"
                    "\tintrepid_scalarMultiplyDataField(outputFields,inputData,inputFields,reciprocal)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (5D array of size [#cubPoints x #fields x #cells x spaceDim x spaceDim] or 4D array of size [#cubPoints x #fields x #cells x spaceDim] or 3D array of size [#cubPoints x #fields x #cells])\n"
                    "\t<2-in> 	  inputData = Input data array (2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in> 	  inputFields = Input fields array ((5D array of size [#cubPoints x #fields x #cells x spaceDim x spaceDim] or 4D array of size [#cubPoints x #fields x #cells x spaceDim] or 3D array of size [#cubPoints x #fields x #cells])\n"
                    "\t<4-in>     reciprocal = If TRUE, divides input fields by the data (instead of multiplying). Default: FALSE");

    // Check the number of input arguments
    if((nInput != 3) && (nInput != 4))
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
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    bool reciprocal = false;
    if(nInput == 4)
    {
        reciprocal = false;
    }
    else if(nInput == 5)
    {
        const std::string si_str(mxArrayToString(pInput[4]));
        if(si_str == "true")
        {
            reciprocal = true;
        }
    }

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    try
    {
        FunctionSpaceTools::scalarMultiplyDataField<double>(outputFields, inputData, inputFields, reciprocal);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_scalarMultiplyDataField.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
