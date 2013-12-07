#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_tensorMultiplyDataData ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::tensorMultiplyDataData.\n\n"
                    "\tintrepid_tensorMultiplyDataData(outputFields,inputDataLeft,inputDataRight,transpose)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (3D array of size [spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in>     inputDataLeft = Input data array (4D array of size [spaceDim x spaceDim x #cubPoints x #cells] or 3D array of size [spaceDim x #cubPoints x #cells] or 2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in>     inputDataRight = Input data array (4D array of size [spaceDim x spaceDim x #cubPoints x #cells] or 3D array of size [spaceDim x #cubPoints x #cells])\n"
                    "\t<4-in>     transpose = If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'\n\n");

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
    Teuchos::Array<int> iDataLeft_dims;
    m2iGetArrayDims(iDataLeft_dims, pInput[1]);
    // Get the dimensions of the reference values array
    Teuchos::Array<int> iDataRight_dims;
    m2iGetArrayDims(iDataRight_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iDataLeft_raw = mxGetPr(pInput[1]);
    double* iDataRight_raw = mxGetPr(pInput[2]);

    char transpose = 'N';
    if(nInput == 3)
    {
        transpose = 'N';
    }
    else if(nInput == 4)
    {
        const std::string si_str(mxArrayToString(pInput[4]));
        if(si_str == "T")
        {
            transpose = 'T';
        }
    }

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputDataLeft(iDataLeft_dims, iDataLeft_raw);
    FieldContainer<double> inputDataRight(iDataRight_dims, iDataRight_raw);

    try
    {
        FunctionSpaceTools::tensorMultiplyDataData<double>(outputFields, inputDataLeft, inputDataRight, transpose);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_tensorMultiplyDataData.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
