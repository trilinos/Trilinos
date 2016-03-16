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
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid2;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getNumCubaturePoints ..... MEX interface for the Intrepid (Trilinos) function Intrepid2::Cubature::getNumPoints.\n\n"
                    "\t[numPts] = intrepid_getNumCubaturePoints(cellType,degree)\n\n"
                    "\t<return> numPts = Number of cubature points (integer)\n"
                    "\t<1-in>   cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<2-in>   degree = Degree of polynomials to be integrated exactly (integer)\n\n");

    // Check the number of input arguments
    if(nInput != 2)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput > 1)
    {
        std::string ioError = descriptor + "Incorrect number of output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the (pointers to) data
    const std::string cell_type(mxArrayToString(pInput[0]));
    const int cubDegree = (int) mxGetScalar(pInput[1]);

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    // Create the cubature factory
    DefaultCubatureFactory<double> cubFactory;
    Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(*(cellTopo.get()), cubDegree);

    pOutput[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    int* numPoints = (int*) mxGetPr(pOutput[0]);

    try
    {
        numPoints[0] = cellCub->getNumPoints();
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getNumCubaturePoints.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
