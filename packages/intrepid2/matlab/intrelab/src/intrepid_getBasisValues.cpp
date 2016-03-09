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
#include "m2i_helpers.hpp"

using namespace Intrepid2;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getBasisValues ..... MEX interface for the Intrepid (Trilinos) function Intrepid2::Basis::getValues.\n\n"
                    "\tintrepid_getBasisValues(outputValues,inputPoints,diffOp,cellType,degree)\n\n"
                    "\t<1-in/out> outputValues = Values of diff. op. applied to basis (variable-rank array)\n"
                    "\t<2-in>     inputPoints = Evaluation points (2D array of size [spaceDim x #inputPoints])\n"
                    "\t<3-in>     diffOp = 'OPERATOR_VALUE' | 'OPERATOR_GRAD' (string)\n"
                    "\t<4-in>     cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<5-in>     degree = Degree (order) of the basis (integer)\n\n");

    // Check the number of input arguments
    if(nInput != 5)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the values array
    Teuchos::Array<int> values_dims;
    m2iGetArrayDims(values_dims, pInput[0]);
    // Get the dimensions of the evaluation point array
    Teuchos::Array<int> points_dims;
    m2iGetArrayDims(points_dims, pInput[1]);

    // Get the (pointers to) data
    const std::string diff_op_str(mxArrayToString(pInput[2]));
    const std::string cell_type_str(mxArrayToString(pInput[3]));
    const int basisDegree = (int) mxGetScalar(pInput[4]);
    double* values_raw = mxGetPr(pInput[0]);
    double* points_raw = mxGetPr(pInput[1]);

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type_str, descriptor);

    Teuchos::RCP<Basis<double, FieldContainer<double> > > basisPtr;
    m2iGetBasis(basisPtr, cell_type_str, basisDegree, descriptor);
    EOperator diff_op = m2iGetDiffOperator(diff_op_str, descriptor);

    FieldContainer<double> basis_values(values_dims, values_raw);
    FieldContainer<double> eval_points(points_dims, points_raw);

    try
    {
        basisPtr->getValues(basis_values, eval_points, diff_op);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getBasisValues.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
