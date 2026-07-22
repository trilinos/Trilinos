// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  L2,      // (value, value)
  VectorWeightedPoisson // (a dot grad, b dot grad)
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
class VectorWeightedPoissonFormulation {
public:
  static const FormulationChoice formulation = VectorWeightedPoisson;
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
