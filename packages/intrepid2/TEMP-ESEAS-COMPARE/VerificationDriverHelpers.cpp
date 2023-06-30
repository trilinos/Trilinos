//
//  VerificationDriverHelpers.cpp
//  VerificationDriver
//
//  Created by Roberts, Nathan V on 6/17/22.
//

#include "Intrepid2_ESEAS_Interface.hpp"
#include "VerificationDriverHelpers.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_TestUtils.hpp"

using DeviceType = Kokkos::DefaultExecutionSpace::device_type;

static const double REL_TOL = 4e-14;
static const double ABS_TOL = 1e-14;

using namespace Intrepid2;

bool valuesMatchToTolerance(double actual, double expected)
{
  bool valuesAreBothSmall = valuesAreSmall(actual, expected, ABS_TOL);
  bool absDiffIsSmall = std::abs(actual-expected) < ABS_TOL;
  bool relDiffIsSmall = essentiallyEqual(actual, expected, REL_TOL);
  return valuesAreBothSmall || absDiffIsSmall || relDiffIsSmall;
}

std::vector< std::vector<int> > getTestCasesUpToDegree(int spaceDim, int minDegree, int polyOrder_x, int polyOrder_y, int polyOrder_z)
{
  std::vector<int> degrees(spaceDim);
  degrees[0] = polyOrder_x;
  if (spaceDim > 1) degrees[1] = polyOrder_y;
  if (spaceDim > 2) degrees[2] = polyOrder_z;
  
  int numCases = degrees[0];
  for (int d=1; d<int(degrees.size()); d++)
  {
    numCases = numCases * (degrees[d] + 1 - minDegree);
  }
  std::vector< std::vector<int> > subBasisDegreeTestCases(numCases);
  for (int caseOrdinal=0; caseOrdinal<numCases; caseOrdinal++)
  {
    std::vector<int> subBasisDegrees(degrees.size());
    int caseRemainder = caseOrdinal;
    for (int d=degrees.size()-1; d>=0; d--)
    {
      int subBasisDegree = caseRemainder % (degrees[d] + 1 - minDegree);
      caseRemainder = caseRemainder / (degrees[d] + 1 - minDegree);
      subBasisDegrees[d] = subBasisDegree + minDegree;
    }
    subBasisDegreeTestCases[caseOrdinal] = subBasisDegrees;
//    std::cout << "Test case " << caseOrdinal << ": {";
//    for (int d=0; d<degrees.size(); d++)
//    {
//      std::cout << subBasisDegrees[d];
//      if (d<degrees.size()-1) std::cout << ",";
//    }
//    std::cout << "}\n";
  }
  return subBasisDegreeTestCases;
}

bool testBasis(shards::CellTopology &cellTopo, Intrepid2::EFunctionSpace fs, Teuchos::FancyOStream &out, bool &success)
{
  using Scalar = double;
  using namespace Intrepid2;
  int numPoints_1D = 15;
  const bool defineVertexFunctions = true; // ESEAS doesn't have a notion of the false case; we don't test that case here
  Teuchos::RCP< Basis<DeviceType,Scalar,Scalar> > basis = getHierarchicalBasis<defineVertexFunctions>(cellTopo, fs, N, N, N);
  
  auto inputPoints = getInputPointsView<double,DeviceType>(cellTopo, numPoints_1D);
  int numPoints = inputPoints.extent_int(0);
  const int spaceDim = inputPoints.extent_int(1);
  
  const bool isHypercube = (cellTopo.getKey() == shards::Line<>::key) || (cellTopo.getKey() == shards::Quadrilateral<>::key) || (cellTopo.getKey() == shards::Hexahedron<>::key);
  const bool isWedge     = (cellTopo.getKey() == shards::Wedge<>::key);
  const bool isPyramid   = (cellTopo.getKey() == shards::Pyramid<>::key);
  
  std::map<int,int> dofMap;
  if (isHypercube && (spaceDim > 1))
  {
    if (spaceDim == 2)
    {
      dofMap = getESEASOrdinalMap(cellTopo, fs, N, N);
    }
    else if (spaceDim == 3)
    {
      dofMap = getESEASOrdinalMap(cellTopo, fs, N, N, N);
    }
  }
  else if (isWedge)
  {
    dofMap = getESEASOrdinalMap(cellTopo, fs, N, N);
  }
  else if (isPyramid)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Pyramid map not yet established");
  }
  else
  {
    int basisCardinality = basis->getCardinality();
    for (int i=0; i<basisCardinality; i++)
    {
      // for triangles and tets, our dof ordering should match that of ESEAS.
      dofMap[i] = i;
    }
  }
  
  ViewType<double,DeviceType> outputValues, outputDerivatives;
  ViewType<double,DeviceType> outputScalarValues, outputVectorValues;
  ViewType<double,DeviceType> outputVectorValues2; // in cases where both value and derivative are vector-valued -- H(curl) -- this will contain the derivative values
  
  Intrepid2::EOperator derivative_op = OPERATOR_VALUE;
  
  // multipliers to account for different reference elements
  // these are non-unitary for derivatives of hypercube topologies
  double scalarMultiplier = 1.0;
  
  std::vector<double> derivativeMultipliers(spaceDim);
  for (int d=0; d<spaceDim; d++)
  {
    derivativeMultipliers[d] = isHypercube ? fromZeroOne_dx(1.0) : 1.0;
  }
  if (cellTopo.getKey() == shards::Wedge<>::key)
  {
    derivativeMultipliers[2] = fromZeroOne_dx(1.0);
  }
  std::vector<double> vectorMultipliers(spaceDim,1.0);
  
  // vectorMultiplier2 is only used in a case where the vector values are always derivative values
  std::vector<double> vectorMultipliers2 = derivativeMultipliers;
  
  int basisCardinality = basis->getCardinality();
  
  if (int(dofMap.size()) != basisCardinality)
  {
    out << "internal test error: dofMap.size() != basisCardinality.\n";
    out << dofMap.size() << " != " << basisCardinality << std::endl;
    success = false;
    return success;
  }
  
  if (fs == FUNCTION_SPACE_HGRAD)
  {
    outputValues      = ViewType<double,DeviceType>("H^1 value output",basisCardinality,numPoints);
    outputDerivatives = ViewType<double,DeviceType>("H^1 derivative output",basisCardinality,numPoints,spaceDim);
    outputScalarValues = outputValues;
    outputVectorValues = outputDerivatives;
    // vector values are derivatives, so they get a multiplier
    vectorMultipliers = derivativeMultipliers;
    derivative_op = OPERATOR_GRAD;
  }
  else if (fs == FUNCTION_SPACE_HCURL)
  {
    if (spaceDim == 2)
    {
      outputValues      = ViewType<double,DeviceType>("H(curl) value output",basisCardinality,numPoints,spaceDim);
      outputDerivatives = ViewType<double,DeviceType>("H(curl) derivative output",basisCardinality,numPoints);
      outputScalarValues = outputDerivatives;
      outputVectorValues = outputValues;
      // scalar values are derivatives, so they get a multiplier
      scalarMultiplier = derivativeMultipliers[0]; // 2D curl derivative multipliers are always isotropic
    }
    else if (spaceDim == 3)
    {
      // CURL and VALUE are both vector-valued in 3D, so we need to treat as a special case
      outputValues      = ViewType<double,DeviceType>("H(curl) value output",basisCardinality,numPoints,spaceDim);
      outputDerivatives = ViewType<double,DeviceType>("H(curl) derivative output",basisCardinality,numPoints,spaceDim);
      outputVectorValues  = outputValues;
      outputVectorValues2 = outputDerivatives;
      if (cellTopo.getBaseKey() == shards::Wedge<>::key)
      {
        // CURL for wedge is a special case where we have to treat the two families differently
        // so for vectorMultipliers2, which corresponds to OPERATOR_CURL, just do unit weighting
        vectorMultipliers2 = std::vector<double>(spaceDim,1.0);
      }
    }
    derivative_op = OPERATOR_CURL;
  }
  else if (fs == FUNCTION_SPACE_HDIV)
  {
    outputValues      = ViewType<double,DeviceType>("H(div) value output",basisCardinality,numPoints,spaceDim);
    outputDerivatives = ViewType<double,DeviceType>("H(div) derivative output",basisCardinality,numPoints);
    outputScalarValues = outputDerivatives;
    outputVectorValues = outputValues;
    // scalar values are derivatives, so they get a multiplier
    scalarMultiplier = derivativeMultipliers[0]; // we special-case the z-derivatives for H(div) wedges, below, in the context of the actual comparison.
    derivative_op = OPERATOR_DIV;
  }
  else if (fs == FUNCTION_SPACE_HVOL)
  {
    outputValues      = ViewType<double,DeviceType>("H(VOL) value output",basisCardinality,numPoints);
    outputScalarValues = outputValues;
    outputVectorValues = outputDerivatives;
    derivative_op = OPERATOR_VALUE; // no derivative op
  }
  
  basis->getValues(outputValues,      inputPoints, OPERATOR_VALUE);
  if (derivative_op != OPERATOR_VALUE)
  {
    basis->getValues(outputDerivatives, inputPoints, derivative_op);
  }
  
  int dofCount = -1;
  int expectedDofCount = -1;
  
  const int size[2] = {N,basisCardinality};
  const int sizeNPlus1[2] = {N+1,basisCardinality}; // used for HVOL
  
  std::cout << "numPoints: " << numPoints << std::endl;
  
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    double scalar_values[basisCardinality];
    double vector_values [basisCardinality*spaceDim];
    double vector_values2[basisCardinality*spaceDim]; // used for derivative values in H(curl)
    const double x_intrepid2 = inputPoints(pointOrdinal,0);
    const double y_intrepid2 = (spaceDim > 1) ? inputPoints(pointOrdinal,1) : -2.0; // -2.0 value just to indicate a value that should not be used
    const double z_intrepid2 = (spaceDim > 2) ? inputPoints(pointOrdinal,2) : -2.0;
    if (cellTopo.getBaseKey() == shards::Quadrilateral<4>::key)
    {
      const double xi[2] = {toZeroOne(x_intrepid2), toZeroOne(y_intrepid2)};
      int polyOrder[5], polyOrderPlusOne[5];
      getQuadPolyOrderESEAS(polyOrder, N);  // used for everything but HVOL
      getQuadPolyOrderESEAS(polyOrderPlusOne, N+1); // used for HVOL
      const int orientations[4] = {0,0,0,0}; // corresponds to all unpermuted edges.  We'll need to add support for permuting edges, but we'll do it in a somewhat different way (not true embedding of the orientations).
      
      if (fs == FUNCTION_SPACE_HGRAD)
      {
        expectedDofCount = (N+1)*(N+1);
        shape2dhquad_(xi, polyOrder, orientations, size, dofCount, scalar_values, vector_values);
      }
      else if (fs == FUNCTION_SPACE_HCURL)
      {
        expectedDofCount = 2*N*(N+1);
        shape2dequad_(xi, polyOrder, orientations, size, dofCount, vector_values, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HDIV)
      {
        expectedDofCount = 2*N*(N+1);
        shape2dvquad_(xi, polyOrder, orientations, size, dofCount, vector_values, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HVOL)
      {
        expectedDofCount = (N+1)*(N+1);
        const int sizeNPlus1[2] = {N+1,basisCardinality};
        shape2dqquad_(xi, polyOrderPlusOne, sizeNPlus1, dofCount, scalar_values);
      }
      else
      {
        success = false;
        out << "FAILURE: function space unsupported by test.  Exiting.\n";
        return success;
      }
      if (expectedDofCount != dofCount)
      {
        success = false;
        out << "FAILURE: expected dof count returned by ESEAS to be " << expectedDofCount << " but it was " << dofCount << std::endl;
        return success;
      }
    }
    else if (cellTopo.getBaseKey() == shards::Hexahedron<8>::key)
    {
      const double xi[3] = {toZeroOne(x_intrepid2), toZeroOne(y_intrepid2), toZeroOne(z_intrepid2)};
      int polyOrder[19], polyOrderPlusOne[19];
      getHexPolyOrderESEAS(polyOrder, N);  // used for everything but HVOL
      getHexPolyOrderESEAS(polyOrderPlusOne, N+1); // used for HVOL
      
      const int edgeOrientations[12] = {0,0,0,0,0,0,0,0,0,0,0,0}; // corresponds to all unpermuted edges.  We don't yet have a mechanism for permuting edges in Intrepid2
      const int faceOrientations[6]  = {0,0,0,0,0,0}; // unpermuted faces
      
      if (fs == FUNCTION_SPACE_HVOL)
      {
        expectedDofCount = (N+1)*(N+1)*(N+1);
        shape3dqhexa_(xi, polyOrderPlusOne, sizeNPlus1, dofCount, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HGRAD)
      {
        //        for (int i=0; i<19; i++)
        //        {
        //          out << "Nord[" << i+1 << "] = " << polyOrder[i] << std::endl;
        //        }
        expectedDofCount = (N+1)*(N+1)*(N+1);
        shape3dhhexa_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, scalar_values, vector_values);
      }
      else if (fs == FUNCTION_SPACE_HCURL)
      {
        expectedDofCount = 3*(N+1)*(N+1)*N;
        shape3dehexa_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, vector_values, vector_values2);
      }
      else if (fs == FUNCTION_SPACE_HDIV)
      {
        expectedDofCount = 3*(N+1)*N*N;
        shape3dvhexa_(xi, polyOrder, faceOrientations, size, dofCount, vector_values, scalar_values);
      }
      else
      {
        success = false;
        out << "FAILURE: function space unsupported by test.  Exiting.\n";
        return success;
      }
    }
    else if (cellTopo.getBaseKey() == shards::Triangle<>::key)
    {
      const double xi[2] {x_intrepid2, y_intrepid2};
      int polyOrder[4], polyOrderPlusOne[4];
      getTriPolyOrderESEAS(polyOrder, N);          // used for everything but HVOL
      getTriPolyOrderESEAS(polyOrderPlusOne, N+1); // used for HVOL
      const int edgeOrientations[3] {0,0,0}; // corresponds to all unpermuted edges.
      
      if (fs == FUNCTION_SPACE_HVOL)
      {
        expectedDofCount = (N+1)*(N+2)/2;
        shape2dqtri_(xi, polyOrderPlusOne, sizeNPlus1, dofCount, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HGRAD)
      {
        expectedDofCount = (N+1)*(N+2)/2;
        shape2dhtri_(xi, polyOrder, edgeOrientations, size, dofCount, scalar_values, vector_values);
      }
      else if (fs == FUNCTION_SPACE_HCURL)
      {
        expectedDofCount = 3 * N + N * (N-1);
        shape2detri_(xi, polyOrder, edgeOrientations, size, dofCount, vector_values, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HDIV)
      {
        expectedDofCount = 3 * N + N * (N-1);
        shape2dvtri_(xi, polyOrder, edgeOrientations, size, dofCount, vector_values, scalar_values);
      }
      else
      {
        success = false;
        out << "FAILURE: function space unsupported by test.  Exiting.\n";
        return success;
      }
    }
    else if (cellTopo.getBaseKey() == shards::Tetrahedron<>::key)
    {
      const double xi[3] {x_intrepid2, y_intrepid2, z_intrepid2};
      int polyOrder[11], polyOrderPlusOne[11];
      getTetPolyOrderESEAS(polyOrder, N);  // used for everything but HVOL
      getTetPolyOrderESEAS(polyOrderPlusOne, N+1); // used for HVOL
      const int edgeOrientations[6] {0,0,0,0,0,0}; // corresponds to all unpermuted edges.  We don't yet have a mechanism for permuting edges in Intrepid2
      const int faceOrientations[4] {0,0,0,0}; // unpermuted faces
      
      if (fs == FUNCTION_SPACE_HVOL)
      {
        expectedDofCount = -1; // TODO: set this appropriately
        shape3dqtet_(xi, polyOrderPlusOne, sizeNPlus1, dofCount, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HGRAD)
      {
        expectedDofCount = (N+1)*(N+2)*(N+3)/6;
        shape3dhtet_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, scalar_values, vector_values);
      }
      else if (fs == FUNCTION_SPACE_HCURL)
      {
        expectedDofCount = -1; // TODO: set this appropriately
        shape3detet_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, vector_values, vector_values2);
      }
      else if (fs == FUNCTION_SPACE_HDIV)
      {
        expectedDofCount = -1; // TODO: set this appropriately
        shape3dvtet_(xi, polyOrder, faceOrientations, size, dofCount, vector_values, scalar_values);
      }
      else
      {
        success = false;
        out << "FAILURE: function space unsupported by test.  Exiting.\n";
        return success;
      }
    }
    else if (cellTopo.getBaseKey() == shards::Wedge<>::key)
    {
      const double xi[3] {x_intrepid2, y_intrepid2, toZeroOne(z_intrepid2)};
      int polyOrder[15], polyOrderPlusOne[15];
      getWedgePolyOrderESEAS(polyOrder, N);  // used for everything but HVOL
      getWedgePolyOrderESEAS(polyOrderPlusOne, N+1); // used for HVOL
      const int edgeOrientations[9] {0,0,0,0,0,0,0,0,0}; // corresponds to all unpermuted edges.  We don't have a mechanism for permuting edges in Intrepid2
      const int faceOrientations[5] {0,0,0,0,0}; // unpermuted faces
      
      if (fs == FUNCTION_SPACE_HVOL)
      {
        expectedDofCount = (N+1)*(N+1)*(N+2)/2; // card of line x card of tri
        shape3dqpris_(xi, polyOrderPlusOne, sizeNPlus1, dofCount, scalar_values);
      }
      else if (fs == FUNCTION_SPACE_HGRAD)
      {
        expectedDofCount = (N+1)*(N+1)*(N+2)/2; // card of line x card of tri
        shape3dhpris_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, scalar_values, vector_values);
      }
      else if (fs == FUNCTION_SPACE_HCURL)
      {
        expectedDofCount = -1; // TODO: set this appropriately
        shape3depris_(xi, polyOrder, edgeOrientations, faceOrientations, size, dofCount, vector_values, vector_values2);
      }
      else if (fs == FUNCTION_SPACE_HDIV)
      {
        expectedDofCount = -1; // TODO: set this appropriately
        shape3dvpris_(xi, polyOrder, faceOrientations, size, dofCount, vector_values, scalar_values);
      }
      else
      {
        success = false;
        out << "FAILURE: function space unsupported by test.  Exiting.\n";
        return success;
      }
    }
    else
    {
      success = false;
      out << "FAILURE: ESEAS interface for " << cellTopo.getName() << " is not yet implemented\n";
      return success;
    }
    
    if ((basisCardinality != expectedDofCount) && (expectedDofCount > 0))
    {
      success = false;
      out << "FAILURE: expected dof count does not match Intrepid2 basis cardinality.  Exiting test." << std::endl;
      return success;
    }
    
    if (dofCount != basisCardinality)
    {
      success = false;
      out << "FAILURE: Intrepid2 basis cardinality does not match ESEAS basis cardinality.  Exiting test." << std::endl;
      return success;
    }
    
    bool enableScalarTests = true;
    bool haveScalarValues = outputScalarValues.size() > 0;
    if (enableScalarTests && haveScalarValues)
    {
      bool pointPassed = true;
      for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
      {
        auto fieldOrdinal_ESEAS = dofMap[fieldOrdinal];
        auto actual_scalar_value   = outputScalarValues(fieldOrdinal,pointOrdinal);
        auto expected_scalar_value = scalarMultiplier * scalar_values[fieldOrdinal_ESEAS];
        
        {
          // ACCOUNTING FOR A DIFFERENCE between ESEAS and our implemention
          // Face fields with non-zero y components in ESEAS are negated relative to our implementation.
          // We tried just adding a -1 factor to the y components in our H(div) implementation, but then we get a
          // failure on the interior degrees of freedom, thanks to the tensor structure.
          // It seems that ESEAS is inconsistent between face weighting and interior weighting;
          // this may be a bug in ESEAS relative to the paper documenting ESEAS -- it
          // looks like the signs should be consistent between face and interior dofs.
          if ((fs == FUNCTION_SPACE_HDIV) && (cellTopo.getBaseKey() == shards::Hexahedron<8>::key)) // HDIV HEX, y component
          {
            auto subcellDim = basis->getDofTag(fieldOrdinal)(0);
            if (subcellDim == 2) // face
            {
              auto subcellOrdinal = basis->getDofTag(fieldOrdinal)(1);
              if ((subcellOrdinal == 0) || (subcellOrdinal == 2)) // top or bottom faces -- these are the ones nonzero in y component
              {
                expected_scalar_value *= -1.0;
              }
            }
          }
          // H(div) on wedges are a special case, because we need to scale all z derivatives
          // Family I does not have z derivatives, and Family II does.  So we need to scale Family II
          const bool opIsDiv = (fs == FUNCTION_SPACE_HDIV);
          if (opIsDiv && (cellTopo.getBaseKey() == shards::Wedge<>::key)) // H(curl) wedge, x or y component
          {
            if (fieldOrdinal >= VALUE_LENGTH_WEDGE_HDIV_FAMILY_I_BASIS) // this fieldOrdinal is in Family II, so its divergenc has a factor involving a z derivative
            {
              expected_scalar_value *= derivativeMultipliers[2]; // z derivative weight
            }
          }
        }
        
        bool values_match = valuesMatchToTolerance(actual_scalar_value, expected_scalar_value);
        
        if (!values_match)
        {
          pointPassed = false;
          if (fs == FUNCTION_SPACE_HCURL)
          {
            // scalar values are the curls
            out << "curl ";
          }
          else if (fs == FUNCTION_SPACE_HDIV)
          {
            // scalar values are the div values
            out << "div ";
          }
          if (spaceDim == 2)
            out << "values for ("  << x_intrepid2 << "," << y_intrepid2 << ") differ for field ordinal " << fieldOrdinal;
          else
            out << "values for ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2 << ") differ for field ordinal " << fieldOrdinal;
          out << ": expected " << expected_scalar_value << "; actual " << actual_scalar_value;
          out << " (diff: " << expected_scalar_value-actual_scalar_value << ")" << std::endl;
          success = false;
        }
      }
      if (!pointPassed)
      {
        if (spaceDim == 2)
          out << "All actual values for the point ("  << x_intrepid2 << "," << y_intrepid2 << ") are as follows:\n";
        else
          out << "All actual values for the point ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2 << ") are as follows:\n";
        for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
        {
          out << fieldOrdinal << ": " << outputScalarValues(fieldOrdinal,pointOrdinal) << std::endl;
        }
        
        if (spaceDim == 2)
          out << "All expected values for the point ("  << x_intrepid2 << "," << y_intrepid2 << ") are as follows:\n";
        else
          out << "All expected values for the point ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2 << ") are as follows:\n";
        for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
        {
          auto fieldOrdinal_ESEAS = dofMap[fieldOrdinal];
          out << fieldOrdinal << " (ESEAS ordinal " << fieldOrdinal_ESEAS << "): " << scalarMultiplier * scalar_values[fieldOrdinal_ESEAS] << std::endl;
        }
      }
    }
    else if (!enableScalarTests)
    {
      static bool haveWarned = false;
      if (!haveWarned)
      {
        out << "Skipping scalar tests (for debugging vector issues).  These should be re-enabled when vector tests pass!\n";
        haveWarned = true;
      }
    }
    bool enableVectorTests = true;
    bool haveVectorValue = (derivative_op != OPERATOR_VALUE); // if the "derivative_op" is "VALUE", this is H(VOL) and no vector values are defined.
    if (enableVectorTests && haveVectorValue)
    {
      std::vector<ViewType<double,DeviceType>> outputVectorValuesIntrepid2 = {outputVectorValues};
      std::vector<std::vector<double>> valueMultipliers = {vectorMultipliers};
      if (outputVectorValues2.size() > 0)
      {
        outputVectorValuesIntrepid2.push_back(outputVectorValues2);
        valueMultipliers.push_back(vectorMultipliers2);
      }
      std::vector<std::vector<std::vector<double>>> vectorValuesESEAS(outputVectorValuesIntrepid2.size());
      for (int valuesContainerOrdinal=0; valuesContainerOrdinal<int(vectorValuesESEAS.size()); valuesContainerOrdinal++)
      {
        vectorValuesESEAS[valuesContainerOrdinal] = std::vector<std::vector<double>>(basisCardinality);
        for (int fieldOrdinal_ESEAS=0; fieldOrdinal_ESEAS<basisCardinality; fieldOrdinal_ESEAS++)
        {
          vectorValuesESEAS[valuesContainerOrdinal][fieldOrdinal_ESEAS] = std::vector<double>(spaceDim);
          for (int d=0; d<spaceDim; d++)
          {
            auto index = flatVectorIndex(spaceDim, fieldOrdinal_ESEAS, d);
            double value = (valuesContainerOrdinal == 0) ? vector_values[index] : vector_values2[index];
            vectorValuesESEAS[valuesContainerOrdinal][fieldOrdinal_ESEAS][d] = value;
          }
        }
      }
      
      for (int valuesContainerOrdinal=0; valuesContainerOrdinal<int(vectorValuesESEAS.size()); valuesContainerOrdinal++)
      {
        auto valuesIntrepid2 = outputVectorValuesIntrepid2[valuesContainerOrdinal];
        auto valuesESEAS     = vectorValuesESEAS[valuesContainerOrdinal];
        bool pointPassed = true;
        auto vectorMultipliersForContainer = valueMultipliers[valuesContainerOrdinal];
        for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
        {
          auto fieldOrdinal_ESEAS = dofMap[fieldOrdinal];
          if (derivative_op != OPERATOR_VALUE) // VALUE is the only one for which no vector-valued container is defined
          {
            for (int d=0; d<spaceDim; d++)
            {
              auto actual_vector_value   = valuesIntrepid2(fieldOrdinal,pointOrdinal,d);
              auto expected_vector_value = vectorMultipliersForContainer[d] * valuesESEAS[fieldOrdinal_ESEAS][d];
              
              {
                // ACCOUNTING FOR A DIFFERENCE between ESEAS and our implemention
                // Face fields with non-zero y components in ESEAS are negated relative to our implementation.
                // We tried just adding a -1 factor to the y components in our H(div) implementation, but then we get a
                // failure on the interior degrees of freedom, thanks to the tensor structure.
                // It seems that ESEAS is inconsistent between face weighting and interior weighting;
                // this may be a bug in ESEAS relative to the paper documenting ESEAS -- it
                // looks like the signs should be consistent between face and interior dofs.
                if ((fs == FUNCTION_SPACE_HDIV) && (cellTopo.getBaseKey() == shards::Hexahedron<8>::key) && (d==1)) // HDIV HEX, y component
                {
                  auto subcellDim = basis->getDofTag(fieldOrdinal)(0);
                  if (subcellDim == 2) // face
                  {
                    expected_vector_value *= -1.0;
                  }
                }
                // H(curl) on wedges are a special case, because we need to scale all z derivatives
                // and z derivatives enter the curl in the (x,y) slots for our Family I but not our Family II.
                const bool opIsCurl = (fs == FUNCTION_SPACE_HCURL) && (valuesContainerOrdinal == 1);
                if (opIsCurl && (cellTopo.getBaseKey() == shards::Wedge<>::key) && (d != 2)) // H(curl) wedge, x or y component
                {
                  if (fieldOrdinal < VALUE_LENGTH_WEDGE_HCURL_FAMILY_I_BASIS) // this fieldOrdinal is in Family I, so its x,y CURL has a factor involving a z derivative
                  {
                    expected_vector_value *= derivativeMultipliers[2]; // z derivative weight
                  }
                }
              }
              
              bool components_match = valuesMatchToTolerance(actual_vector_value, expected_vector_value);
              
              if (!components_match)
              {
                pointPassed = false;
                if (fs == FUNCTION_SPACE_HGRAD)
                {
                  out << "gradient ";
                }
                else if ((fs == FUNCTION_SPACE_HCURL) && (valuesContainerOrdinal==1))
                {
                  out << "curl ";
                }
                if (spaceDim == 2)
                {
                  out << "values for ("  << x_intrepid2 << "," << y_intrepid2;
                }
                else
                {
                  out << "values for ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2;
                }
                
                out << ") differ in component " << d << " for field ordinal " << fieldOrdinal;
                out << ": expected " << expected_vector_value << "; actual " << actual_vector_value;
                out << " (diff: " << expected_vector_value-actual_vector_value << ")" << std::endl;
                success = false;
              }
            }
          }
        }
        if (!pointPassed)
        {
          if (spaceDim == 2)
            out << "All actual values for the point ("  << x_intrepid2 << "," << y_intrepid2 << ") are as follows:\n";
          else
            out << "All actual values for the point ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2 << ") are as follows:\n";
          for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
          {
            out << fieldOrdinal << ": (";
            for (int d=0; d<spaceDim; d++)
            {
              out << valuesIntrepid2(fieldOrdinal,pointOrdinal,d);
              if (d<spaceDim-1) out << ",";
            }
            out << ")\n";
          }
          
          if (spaceDim == 2)
            out << "All expected values for the point ("  << x_intrepid2 << "," << y_intrepid2 << ") are as follows:\n";
          else
            out << "All expected values for the point ("  << x_intrepid2 << "," << y_intrepid2 << "," << z_intrepid2 << ") are as follows:\n";
          for (int fieldOrdinal=0; fieldOrdinal<basisCardinality; fieldOrdinal++)
          {
            auto fieldOrdinal_ESEAS = dofMap[fieldOrdinal];
            out << fieldOrdinal << " (ESEAS ordinal " << fieldOrdinal_ESEAS << "): (";
            for (int d=0; d<spaceDim; d++)
            {
              out << vectorMultipliersForContainer[d] * valuesESEAS[fieldOrdinal_ESEAS][d];
              if (d<spaceDim-1) out << ",";
            }
            out << ")\n";
          }
        }
      }
    }
    else if (!enableVectorTests)
    {
      static bool haveWarned = false;
      if (!haveWarned)
      {
        out << "Skipping vector tests (for debugging scalar issues).  These should be re-enabled when scalar tests pass!\n";
        haveWarned = true;
      }
    }
  }
  if (success == false)
  {
    //then to assist in assessment/debugging, print a report of which dofs go where
    out << "testBasis() failed for basis named " << basis->getName() << std::endl;
    out << "******    Dofs in this basis:    ******\n";
    
    
    /** \brief  Retrieves all DoF tags.
     
     \return reference to a vector of vectors with dimensions (basisCardinality_, 4) such that \n
     \li     element [DofOrd][0] = tag field 0 for the DoF with the specified ordinal
     \li     element [DofOrd][1] = tag field 1 for the DoF with the specified ordinal
     \li     element [DofOrd][2] = tag field 2 for the DoF with the specified ordinal
     \li     element [DofOrd][3] = tag field 3 for the DoF with the specified ordinal
     */
    
    using std::setw;
    using std::endl;
    auto tagToOrdinal = basis->getAllDofTags();
    const int wCol0 = 12;
    const int wCol1 = 12;
    const int wCol2 = 12;
    const int wCol3 = 20;
    const int wCol4 = 25;
    out << setw(wCol0);
    out << "Dof Ordinal"  << setw(wCol1);
    out << "Subcell Dim"  << setw(wCol2);
    out << "Subcell Ord"  << setw(wCol3);
    out << "Dof Ord in Subcell" << setw(wCol4);
    out << "Total Dofs in Subcell";
    out << endl;
    for (int dofOrdinal=0; dofOrdinal<tagToOrdinal.extent_int(0); dofOrdinal++)
    {
      out << setw(wCol0);
      out << dofOrdinal << setw(wCol1);
      out << tagToOrdinal(dofOrdinal,0) << setw(wCol2);
      out << tagToOrdinal(dofOrdinal,1) << setw(wCol3);
      out << tagToOrdinal(dofOrdinal,2) << setw(wCol4);
      out << tagToOrdinal(dofOrdinal,3);
      out << endl;
    }
  }
  
  return success;
}
