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
            ("\nintrepid_HGRADtransformVALUE ..... MEX interface for the Intrepid (Trilinos) function Intrepid2::FunctionSpaceTools::HGRADtransformVALUE.\n\n"
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
