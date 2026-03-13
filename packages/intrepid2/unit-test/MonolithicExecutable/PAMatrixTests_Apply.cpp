// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Comparing PAMatrix apply versus standard matrix assembly and apply().
    \author Nathan V. Roberts
*/

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellGeometryTestUtils.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>
#include <Intrepid2_IntegrationTools.hpp>
#include <Intrepid2_Kernels.hpp>
#include <Intrepid2_NodalBasisFamily.hpp>
#include <Intrepid2_TensorArgumentIterator.hpp>
#include <Intrepid2_TestUtils.hpp>

#include <Intrepid2_CellGeometry.hpp>
#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

#include "PAMatrixAssembly.hpp"

#include "StructuredIntegrationTests_TagDefs.hpp"
#include "StructuredIntegrationTests_Utils.hpp"

namespace
{
using namespace Intrepid2;
template<typename DeviceType, typename Scalar>
using PAMatrix = ::Intrepid2::PAMatrix<DeviceType,Scalar>;

template<class Scalar, typename DeviceType>
void testPAMatrixApply(PAMatrix<DeviceType,Scalar> &paMatrix, const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
{
  // compare columns of the fully assembled matrix with the action of PAMatrix::apply() on unit vectors with a "1" in corresponding column position.
  
  const bool printDebugOutput = false;
  if (printDebugOutput)
  {
    // DEBUGGING
    using namespace std;
    
    {
      // print transformed basis values
      auto basisValuesLeft = paMatrix._basisValuesLeft;
      const int numCells = basisValuesLeft.extent_int(0);
      const int numFieldsLeft = basisValuesLeft.extent_int(1);
      const int numPoints = basisValuesLeft.extent_int(2);
      const int spaceDim  = basisValuesLeft.extent_int(3);
      int rank = basisValuesLeft.rank();
      cout << "physical basisValuesLeft: { ";
      for (int c=0; c<numCells; c++)
      {
        cout << "{ ";
        for (int f=0; f<numFieldsLeft; f++)
        {
          cout << "{ ";
          for (int p=0; p<numPoints; p++)
          {
            if (rank == 3)
            {
              cout << basisValuesLeft(c,f,p);
            }
            else
            {
              cout << "(";
              for (int d=0; d<spaceDim; d++)
              {
                cout << basisValuesLeft(c,f,p,d);
                if (d<spaceDim-1) cout << ",";
              }
              cout << ")";
            }
            if (p < numPoints-1) cout << ", ";
          }
          cout << "}";
          if (f<numFieldsLeft-1) cout << ", ";
        }
        cout << "}";
        if (c<numCells-1) cout << ", ";
      }
      cout << "}\n";
      
      rank = basisValuesLeft.basisValues().rank();
      cout << "reference basisValuesLeft: { ";
      for (int f=0; f<numFieldsLeft; f++)
      {
        cout << "{ ";
        for (int p=0; p<numPoints; p++)
        {
          if (rank == 2)
          {
            cout << basisValuesLeft.basisValues()(f,p);
          }
          else
          {
            cout << "(";
            for (int d=0; d<spaceDim; d++)
            {
              cout << basisValuesLeft.basisValues()(f,p,d);
              if (d<spaceDim-1) cout << ",";
            }
            cout << ")";
          }
          if (p < numPoints-1) cout << ", ";
        }
        cout << "}";
        if (f<numFieldsLeft-1) cout << ", ";
      }
      
      auto transform = basisValuesLeft.transform();
      cout << "basisValuesLeft transform: { ";
      if (transform.rank() == 2)
      {
        // (C,P)
        for (int c=0; c<numCells; c++)
        {
          cout << "{ ";
          for (int p=0; p<numPoints; p++)
          {
            cout << transform(c,p) << ", ";
          }
          cout << "}";
          if (c<numCells-1) cout << ", ";
        }
        cout << "}\n";
      }
      else if (transform.rank() == 4)
      {
        // (C,P,D,D)
        const int spaceDim = transform.extent_int(2);
        cout << "transform: { ";
        for (int c=0; c<numCells; c++)
        {
          cout << "{ ";
          for (int p=0; p<numPoints; p++)
          {
            cout << "[";
            for (int d1=0; d1<spaceDim; d1++)
            {
              cout << "[";
              for (int d2=0; d2<spaceDim; d2++)
              {
                cout << transform(c,p,d1,d2) << " ";
              }
              cout << "]";
            }
            cout << "]";
            if (p<numPoints-1) cout << ", ";
          }
          cout << "}";
          if (c<numCells-1) cout << ", ";
        }
        cout << "}\n";
      }

      auto basisValuesRight = paMatrix._basisValuesRight;
      const int numFieldsRight = basisValuesRight.extent_int(1);
      rank = basisValuesRight.rank();
      cout << "physical basisValuesRight: { ";
      for (int c=0; c<numCells; c++)
      {
        cout << "{ ";
        for (int f=0; f<numFieldsRight; f++)
        {
          cout << "{ ";
          for (int p=0; p<numPoints; p++)
          {
            if (rank == 3)
            {
              cout << basisValuesRight(c,f,p);
            }
            else
            {
              cout << "(";
              for (int d=0; d<spaceDim; d++)
              {
                cout << basisValuesRight(c,f,p,d);
                if (d<spaceDim-1) cout << ",";
              }
              cout << ")";
            }
            if (p < numPoints-1) cout << ", ";
          }
          cout << "}";
          if (f<numFieldsRight-1) cout << ", ";
        }
        cout << "}";
        if (c<numCells-1) cout << ", ";
      }
      cout << "}\n";
      
      rank = int(basisValuesRight.basisValues().rank());
      cout << "reference basisValuesRight: { ";
      for (int f=0; f<numFieldsRight; f++)
      {
        cout << "{ ";
        for (int p=0; p<numPoints; p++)
        {
          if (rank == 2)
          {
            cout << basisValuesRight.basisValues()(f,p);
          }
          else
          {
            cout << "(";
            for (int d=0; d<spaceDim; d++)
            {
              cout << basisValuesRight.basisValues()(f,p,d);
              if (d<spaceDim-1) cout << ",";
            }
            cout << ")";
          }
          if (p < numPoints-1) cout << ", ";
        }
        cout << "}";
        if (f<numFieldsRight-1) cout << ", ";
      }
      cout << "}\n";
    }
    
    auto basisValuesLeft      = paMatrix._basisValuesLeft.basisValues().vectorData(); // ref space values
    const int numComponents   = basisValuesLeft.numComponents();
    const int numFamiliesLeft = basisValuesLeft.numFamilies();
    for (int family=0; family<numFamiliesLeft; family++)
    {
      for (int d=0; d<numComponents; d++)
      {
        cout << "basisValuesLeft family " << family << "(" << d << "):\n";
        
        auto tensorData = basisValuesLeft.getComponent(family,d);
        
        int numTensorComponents = tensorData.numTensorComponents();
        for (int r=0; r<numTensorComponents; r++)
        {
          auto data = tensorData.getTensorComponent(r);
          
          const int numFields = data.extent_int(0);
          const int numPoints = data.extent_int(1);
          
          cout << "{ ";
          for (int f=0; f<numFields; f++)
          {
            cout << "{ ";
            for (int p=0; p<numPoints; p++)
            {
              if (data.rank() == 2)
              {
                cout << data(f,p);
              }
              else
              {
                const int spaceDim = data.extent_int(2);
                if (spaceDim == 1)
                {
                  cout << data(f,p,0);
                }
                else
                {
                  cout << "(";
                  for (int d1=0; d1<spaceDim; d1++)
                  {
                    cout << tensorData(f,p,d1);
                    if (d1 < spaceDim-1) cout << ",";
                  }
                  cout << ")";
                }
              }
              if (p<numPoints-1) cout << ", ";
            }
            cout << "}";
            if (f < numFields-1) cout << ", ";
          }
          cout << " }";
          if (r<numTensorComponents-1) cout << " x ";
          else cout << "\n";
        }
      }
    }
    
//    auto basisValuesRight = paMatrix._basisValuesRight;
//    const int numFieldsRight = basisValuesRight.extent_int(1);
//    cout << "basisValuesRight: { ";
//    for (int c=0; c<numCells; c++)
//    {
//      cout << "{ ";
//      for (int f=0; f<numFieldsRight; f++)
//      {
//        cout << "{ ";
//        for (int p=0; p<numPoints; p++)
//        {
//          cout << basisValuesRight(c,f,p) << ", ";
//        }
//        cout << "}";
//        if (f<numFieldsRight-1) cout << ", ";
//      }
//      cout << "}";
//      if (c<numCells-1) cout << ", ";
//    }
//    cout << "}\n";
//    
    auto transform = paMatrix._composedWeightedTransform;
    const int numCells  = transform.extent_int(0);
    const int numPoints = transform.extent_int(1);
    if (transform.rank() == 2)
    {
      // assuming (C,P)
      cout << "composed weighted transform: { ";
      for (int c=0; c<numCells; c++)
      {
        cout << "{ ";
        for (int p=0; p<numPoints; p++)
        {
          cout << transform(c,p) << ", ";
        }
        cout << "}";
        if (c<numCells-1) cout << ", ";
      }
      cout << "}\n";
    }
    else if (transform.rank() == 4)
    {
      const int spaceDim = transform.extent_int(2);
      cout << "transform: { ";
      for (int c=0; c<numCells; c++)
      {
        cout << "{ ";
        for (int p=0; p<numPoints; p++)
        {
          cout << "[";
          for (int d1=0; d1<spaceDim; d1++)
          {
            cout << "[";
            for (int d2=0; d2<spaceDim; d2++)
            {
              cout << transform(c,p,d1,d2) << " ";
            }
            cout << "]";
          }
          cout << "]";
          if (p<numPoints-1) cout << ", ";
        }
        cout << "}";
        if (c<numCells-1) cout << ", ";
      }
      cout << "}\n";
    }
  }
  
  auto fullMatrix = paMatrix.allocateMatrixStorage();
  paMatrix.assemble(fullMatrix);
  
  auto fullMatrixView = fullMatrix.getUnderlyingView(); // I believe this will always have shape (C,F1,F2) for our present inputs
  TEUCHOS_TEST_FOR_EXCEPTION(fullMatrixView.rank() != 3, std::invalid_argument, "underlying view does not have the expected rank");
  
  // construct a unit vector to extract a corresponding column
  const int numCells = fullMatrix.extent_int(0);
  const int numRows  = fullMatrix.extent_int(1); // (C,F1,F2)
  const int numCols  = fullMatrix.extent_int(2); // (C,F1,F2)
  
  ScalarView<Scalar, DeviceType>  inputVector( "inputVector",numCells,numCols); // C,F2
  ScalarView<Scalar, DeviceType> outputVector("outputVector",numCells,numRows); // C,F1
  
  Kokkos::deep_copy(inputVector,0.0);
  auto workspace = paMatrix.allocateWorkspace(numCells);
  
  auto divideOutputVectorValuesByTwo = [&]()
  {
    using ExecutionSpace = typename DeviceType::execution_space;
    auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{numCells,numRows});
    Kokkos::parallel_for("testPAMatrixApply(): subtractOutputVectorValues", policy,
                         KOKKOS_LAMBDA(const ordinal_type &c, const ordinal_type &row)
                         {
      outputVector(c,row) *= 0.5;
    });
    ExecutionSpace().fence();
  };
  
  for (int colOrdinal=0; colOrdinal<numCols; colOrdinal++)
  {
    auto colSubView = Kokkos::subview(inputVector, Kokkos::ALL(), colOrdinal);
    Kokkos::deep_copy(colSubView, 1.0);
    
    paMatrix.apply(outputVector, inputVector, workspace);
  
    auto fullMatrixColumn = Kokkos::subview(fullMatrixView, Kokkos::ALL(), Kokkos::ALL(), colOrdinal);
    
    out << "checking col: " << colOrdinal << std::endl;
    testFloatingEquality2(fullMatrixColumn, outputVector, relTol, absTol, out, success, "fullMatrixColumn", "outputVector");
    
    // now do it again but with sumInto = true
    const bool sumInto = true;
    paMatrix.apply(outputVector, inputVector, workspace, sumInto);
    // now we have double the correct value; divide by two
    divideOutputVectorValuesByTwo();
    
    out << "checking col: " << colOrdinal << " (sumInto = true) " << std::endl;
    testFloatingEquality2(fullMatrixColumn, outputVector, relTol, absTol, out, success, "fullMatrixColumn", "outputVector");
    
    // zero out the column entry we changed for next pass through the loop
    Kokkos::deep_copy(colSubView, 0.0);
  }
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
void testPAMatrixApply(const int &meshWidth, const int &worksetSize,
                       const EFunctionSpace &fs1, const EOperator &op1, const int &p1, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight1,
                       const EFunctionSpace &fs2, const EOperator &op2, const int &p2, Teuchos::RCP< Kokkos::Array<PointScalar,spaceDim> > vectorWeight2,
                       const double &relTol, const double &absTol, Teuchos::FancyOStream &out, bool &success)
{
  // compare columns of the fully assembled matrix with the action of PAMatrix::apply() on unit vectors with a "1" in corresponding column position.
  
  using namespace std;
  
  Kokkos::Array<int,spaceDim> gridCellCounts;
  for (int d=0; d<spaceDim; d++)
  {
    gridCellCounts[d] = meshWidth;
  }
  
  auto geometry = getMesh<PointScalar, spaceDim, DeviceType>(Standard, gridCellCounts);
  shards::CellTopology cellTopo = geometry.cellTopology();
  
  double flopCountIntegration = 0, flopCountJacobian = 0;
  auto paMatrix = constructPAMatrix<Scalar,BasisFamily>(geometry,
                                                        p1, fs1, op1, vectorWeight1,
                                                        p2, fs2, op2, vectorWeight2);
  
  testPAMatrixApply<Scalar, DeviceType>(paMatrix, relTol, absTol, out, success);
}

template<class Scalar, class BasisFamily, class PointScalar, int spaceDim, typename DeviceType>
void testPAMatrixApply(const int &meshWidth, const int &worksetSize,
                                             const EFunctionSpace &fs1, const EOperator &op1, const int &p1,
                                             const EFunctionSpace &fs2, const EOperator &op2, const int &p2,
                                             const double &relTol, const double &absTol,
                                             Teuchos::FancyOStream &out, bool &success)
{
  testPAMatrixApply<Scalar, BasisFamily, PointScalar, spaceDim, DeviceType>(meshWidth, worksetSize,
                                                                            fs1, op1, p1, Teuchos::null,
                                                                            fs2, op2, p2, Teuchos::null,
                                                                            relTol, absTol, out, success);
}

// MARK: Apply_D1_P1_P1
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

// MARK: Apply_D1_P1_P2
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D1_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D1_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 1;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D2_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D2_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 2;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 2;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D3_P1_P2, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 2;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, Apply_D3_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, fs2, op2, p2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorWeighted_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  
  weight = 0.5;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight *= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorWeighted_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  
  weight = 0.5;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight *= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorWeighted_D2_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  
  weight = 0.5;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight *= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyScalarAgainstVectorDotVector_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyScalarAgainstVectorDotVector_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyScalarAgainstVectorDotVector_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight1; // no vector weight on scalar term
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight /= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorDotVectorAgainstScalar_D1_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 1;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorDotVectorAgainstScalar_D2_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 2;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorDotVectorAgainstScalar_D3_P1_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 1;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  Teuchos::RCP<Kokkos::Array<double,spaceDim> > vectorWeight2; // no vector weight on scalar term
  
  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PAMatrix, ApplyVectorWeighted_D3_P2_P1, FS1Tag, Op1Tag, FS2Tag, Op2Tag)
{
  using DataScalar  = double;
  using PointScalar = double;
  const int meshWidth = 1;
  const int spaceDim = 3;
  const int p1 = 2;
  const int p2 = 1;
  const int worksetSize = meshWidth;
  
  auto vectorWeight1 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  auto vectorWeight2 = Teuchos::rcp(new Kokkos::Array<double,spaceDim>);
  
  double weight = 1.0;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight1)[d] = weight;
    weight /= 2.0;
  }
  
  weight = 0.5;
  for (int d=0; d<spaceDim; d++)
  {
    (*vectorWeight2)[d] = weight;
    weight *= 2.0;
  }

  using DeviceType = DefaultTestDeviceType;
  using BasisFamily = DerivedNodalBasisFamily<DeviceType>;
  
  const EFunctionSpace fs1 = FS1Tag::functionSpace;
  const EFunctionSpace fs2 = FS2Tag::functionSpace;
  const EOperator op1 = Op1Tag::op;
  const EOperator op2 = Op2Tag::op;
  
  double relTol = 1e-12;
  double absTol = 1e-12;
  
  testPAMatrixApply<DataScalar, BasisFamily, PointScalar, spaceDim, DeviceType>
    (meshWidth, worksetSize, fs1, op1, p1, vectorWeight1, fs2, op2, p2, vectorWeight2, relTol, absTol, out, success);
}

// asymmetric tests (mostly -- a couple symmetric ones tossed in as sanity checks on the test itself)

// 1D tests: H(grad) and H(vol) bases defined
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P2, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P1_P2, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P2_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D1_P2_P1, HVOL,  VALUE, HGRAD, VALUE)

// 1D vector-weighted test
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D1_P1_P1, HGRAD, GRAD, HGRAD, GRAD)

// 1D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D1_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D1_P1_P1, HGRAD, VALUE, HGRAD, GRAD)

// 1D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D1_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D1_P1_P1, HGRAD, GRAD, HGRAD, VALUE)

// 2D tests: curls of H(curl) are scalars.
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, VALUE, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HGRAD, VALUE, HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HCURL, CURL,  HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D2_P1_P2, HVOL,  VALUE, HGRAD, VALUE)

// 2D vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D2_P1_P1, HGRAD, GRAD, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D2_P2_P1, HGRAD, GRAD, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D2_P1_P1, HCURL, VALUE, HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D2_P2_P1, HCURL, VALUE, HDIV,  VALUE)

// 2D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D2_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D2_P1_P1, HGRAD, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D2_P1_P1, HGRAD, VALUE, HDIV, VALUE)

// 2D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D2_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D2_P1_P1, HGRAD, GRAD, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D2_P1_P1, HDIV, VALUE, HGRAD, VALUE)

// 3D tests: curls of H(curl) are vectors
// p1, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P1, HVOL,  VALUE, HGRAD, VALUE)
// p2, p1:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P2_P1, HVOL,  VALUE, HGRAD, VALUE)
// p1, p2:
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HGRAD, GRAD,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HGRAD, GRAD,  HCURL, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HGRAD, GRAD,  HCURL, CURL)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HGRAD, VALUE, HDIV,  DIV)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HDIV,  DIV,   HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HCURL, CURL,  HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, Apply_D3_P1_P2, HVOL,  VALUE, HGRAD, VALUE)

// 3D vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D3_P2_P1, HGRAD, GRAD,  HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D3_P2_P1, HCURL, VALUE, HDIV,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorWeighted_D3_P2_P1, HCURL, CURL,  HGRAD, GRAD)

// 3D scalar against vector-weighted tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D3_P1_P1, HVOL, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D3_P1_P1, HGRAD, VALUE, HGRAD, GRAD)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyScalarAgainstVectorDotVector_D3_P1_P1, HGRAD, VALUE, HDIV, VALUE)

// 3D vector-weighted against scalar tests
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D3_P1_P1, HGRAD, GRAD, HVOL,  VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D3_P1_P1, HGRAD, GRAD, HGRAD, VALUE)
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PAMatrix, ApplyVectorDotVectorAgainstScalar_D3_P1_P1, HDIV, VALUE, HGRAD, VALUE)

} // anonymous namespace
