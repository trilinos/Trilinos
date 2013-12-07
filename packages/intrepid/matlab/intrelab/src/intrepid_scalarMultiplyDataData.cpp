#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_scalarMultiplyDataData ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::scalarMultiplyDataData.\n\n"
                    "\tintrepid_scalarMultiplyDataData(outputData,inputDataLeft,inputDataRight,reciprocal)\n\n"
                    "\t<1-in/out> outputData = Output data array (5D array of size [#cubPoints x #fields x #cells x spaceDim x spaceDim] or 4D array of size [#cubPoints x #fields x #cells x spaceDim] or 3D array of size [#cubPoints x #fields x #cells])\n"
                    "\t<2-in> 	  inputDataLeft = Input data array (2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in> 	  inputDataRight = Input data array (4D array of size [#cubPoints x #cells x spaceDim x spaceDim] or 3D array of size [#cubPoints x #cells x spaceDim] or 2D array of size [#cubPoints x #cells])\n"
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
    Teuchos::Array<int> oData_dims;
    m2iGetArrayDims(oData_dims, pInput[0]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> iDataLeft_dims;
    m2iGetArrayDims(iDataLeft_dims, pInput[1]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> iDataRight_dims;
    m2iGetArrayDims(iDataRight_dims, pInput[2]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iDataLeft_raw = mxGetPr(pInput[1]);
    double* iDataRight_raw = mxGetPr(pInput[2]);

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

    FieldContainer<double> outputData(oData_dims, oData_raw);
    FieldContainer<double> inputDataLeft(iDataLeft_dims, iDataLeft_raw);
    FieldContainer<double> inputDataRight(iDataRight_dims, iDataRight_raw);

    try
    {
        FunctionSpaceTools::scalarMultiplyDataData<double>(outputData, inputDataLeft, inputDataRight, reciprocal);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_scalarMultiplyDataData.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
