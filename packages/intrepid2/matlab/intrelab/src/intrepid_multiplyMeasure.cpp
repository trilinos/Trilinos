// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "mex.h"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid2;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_multiplyMeasure ..... MEX interface for the Intrepid (Trilinos) function Intrepid2::FunctionSpaceTools::multiplyMeasure.\n\n"
                    "\tintrepid_multiplyMeasure(outFields,inMeasures,inFields)\n\n"
                    "\t<1-in/out> outFields = Output fields (variable-size array [dim1(optional) x dim2(optional) x #points x #basisFields #cells])\n"
                    "\t<2-in>     inMeasures = Cell measures (2D array of size [#points x #cells])\n"
                    "\t<3-in/out> inFields = Input fields (variable-size array [dim1(optional) x dim2(optional) x #points x #basisFields #cells(optional)])\n\n");

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

    // Get the dimensions of the output fields array
    Teuchos::Array<int> outFields_dims;
    m2iGetArrayDims(outFields_dims, pInput[0]);
    // Get the dimensions of the input measures array
    Teuchos::Array<int> inMeasures_dims;
    m2iGetArrayDims(inMeasures_dims, pInput[1]);
    // Get the dimensions of the input fields array
    Teuchos::Array<int> inFields_dims;
    m2iGetArrayDims(inFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* outFields_raw = mxGetPr(pInput[0]);
    double* inMeasures_raw = mxGetPr(pInput[1]);
    double* inFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> outFields(outFields_dims, outFields_raw);
    FieldContainer<double> inMeasures(inMeasures_dims, inMeasures_raw);
    FieldContainer<double> inFields(inFields_dims, inFields_raw);

    try
    {
        FunctionSpaceTools::multiplyMeasure<double>(outFields, inMeasures, inFields);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_multiplyMeasure.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
