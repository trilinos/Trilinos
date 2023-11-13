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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file   StructuredIntegrationTests_TagDefs.hpp
    \brief  Definitions of tags and enums for use in StructuredIntegrationTests.
    \author Nathan V. Roberts
*/

#include "Intrepid2_Types.hpp"

using namespace Intrepid2;

enum FormulationChoice
{
  Poisson, // (grad, grad)
  Hgrad,   // (grad, grad) + (value, value)
  Hdiv,    // (div, div)   + (value, value)
  Hcurl,   // (curl, curl) + (value, value)
  L2       // (value, value)
};

enum AlgorithmChoice
{
  Standard,
  AffineNonTensor,
  NonAffineTensor,
  AffineTensor,
  DiagonalJacobian,
  Uniform
};

enum BasisFamilyChoice
{
  Nodal,
  Hierarchical,
  Serendipity
};

// tags to allow us to use templated Teuchos tests
class PoissonFormulation {
public:
  static const FormulationChoice formulation = Poisson;
};
class HgradFormulation {
public:
  static const FormulationChoice formulation = Hgrad;
};
class HdivFormulation {
public:
  static const FormulationChoice formulation = Hdiv;
};
class HcurlFormulation {
public:
  static const FormulationChoice formulation = Hcurl;
};
class L2Formulation {
public:
  static const FormulationChoice formulation = L2;
};
class StandardAlgorithm
{
public:
  static const AlgorithmChoice algorithm = Standard;
};
class AffineNonTensorAlgorithm
{
public:
  static const AlgorithmChoice algorithm = AffineNonTensor;
};
class NonAffineTensorAlgorithm
{
public:
  static const AlgorithmChoice algorithm = NonAffineTensor;
};
class AffineTensorAlgorithm
{
public:
  static const AlgorithmChoice algorithm = AffineTensor;
};
class UniformAlgorithm
{
public:
  static const AlgorithmChoice algorithm = Uniform;
};
class DiagonalJacobianAlgorithm // note that DiagonalJacobian is not yet supported by getMesh()
{
public:
  static const AlgorithmChoice algorithm = DiagonalJacobian;
};
class D1
{
public:
  static const int spaceDim = 1;
};
class D2
{
public:
  static const int spaceDim = 2;
};
class D3
{
public:
  static const int spaceDim = 3;
};
class P1
{
public:
  static const int polyOrder = 1;
};
class P2
{
public:
  static const int polyOrder = 2;
};
class P3
{
public:
  static const int polyOrder = 3;
};
class P4
{
public:
  static const int polyOrder = 4;
};

class HGRAD
{
public:
  static const Intrepid2::EFunctionSpace functionSpace = Intrepid2::FUNCTION_SPACE_HGRAD;
};
class HDIV
{
public:
  static const Intrepid2::EFunctionSpace functionSpace = Intrepid2::FUNCTION_SPACE_HDIV;
};
class HCURL
{
public:
  static const Intrepid2::EFunctionSpace functionSpace = Intrepid2::FUNCTION_SPACE_HCURL;
};
class HVOL
{
public:
  static const Intrepid2::EFunctionSpace functionSpace = Intrepid2::FUNCTION_SPACE_HVOL;
};

class GRAD
{
public:
  static const Intrepid2::EOperator op = Intrepid2::OPERATOR_GRAD;
};
class DIV
{
public:
  static const Intrepid2::EOperator op = Intrepid2::OPERATOR_DIV;
};
class CURL
{
public:
  static const Intrepid2::EOperator op = Intrepid2::OPERATOR_CURL;
};
class VALUE
{
public:
  static const Intrepid2::EOperator op = Intrepid2::OPERATOR_VALUE;
};
