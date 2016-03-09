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
#include <Intrepid2_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid2;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n geometricBmat ..... MEX interface for function geometricBmat.\n\n"
                    "\t geometricBmat(outputFields,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [#spaceDim^2 x #cubPoints x numFields x #cells])\n"
                    "\t<2-in> 	   inputFields = Input fields array (4D array of size [#spaceDim x #cubPoints x #numFields x #cells])\n\n");

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

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[1]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iFields_raw = mxGetPr(pInput[1]);

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    // get sizes
    const int numCells = inputFields.dimension(0);
    const int numNodes = inputFields.dimension(1);
    const int numQPs = inputFields.dimension(2);
    const int spaceDim = inputFields.dimension(3);

    for(unsigned cell = 0; cell < numCells; ++cell)
    {
        for(unsigned node = 0; node < numNodes; ++node)
        {
            for(unsigned qp = 0; qp < numQPs; ++qp)
            {
                switch(spaceDim)
                {
                    case 2:
                    {
                        outputFields(cell, spaceDim * node, qp, 0) = inputFields(cell, node, qp, 0);
                        outputFields(cell, spaceDim * node, qp, 1) = inputFields(cell, node, qp, 1);
                        outputFields(cell, spaceDim * node + 1, qp, 2) = inputFields(cell, node, qp, 0);
                        outputFields(cell, spaceDim * node + 1, qp, 3) = inputFields(cell, node, qp, 1);
                        break;
                    }
                    case 3:
                    {
                        outputFields(cell, spaceDim * node, qp, 0) = inputFields(cell, node, qp, 0);
                        outputFields(cell, spaceDim * node, qp, 1) = inputFields(cell, node, qp, 1);
                        outputFields(cell, spaceDim * node, qp, 2) = inputFields(cell, node, qp, 2);
                        outputFields(cell, spaceDim * node + 1, qp, 3) = inputFields(cell, node, qp, 0);
                        outputFields(cell, spaceDim * node + 1, qp, 4) = inputFields(cell, node, qp, 1);
                        outputFields(cell, spaceDim * node + 1, qp, 5) = inputFields(cell, node, qp, 2);
                        outputFields(cell, spaceDim * node + 2, qp, 6) = inputFields(cell, node, qp, 0);
                        outputFields(cell, spaceDim * node + 2, qp, 7) = inputFields(cell, node, qp, 1);
                        outputFields(cell, spaceDim * node + 2, qp, 8) = inputFields(cell, node, qp, 2);
                        break;
                    }
                } // switch statement
            } // P-loop
        } // F-loop
    } // C-loop
}
