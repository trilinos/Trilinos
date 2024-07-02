// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_UnitTestRepository.hpp"

#include "Kokkos_Core.hpp"

#include <fstream>

const int numFields = 8;
const int numNodes = 8;
const int numDims  = 3;
const int numIntg  = 8;

using DeviceType = Kokkos::DefaultExecutionSpace::device_type;

int worksetSize = 512;

Kokkos::DynRankView<double, DeviceType> points;   // (P,D) or (C,P,D), allocated/initialized in main below
Kokkos::DynRankView<double, DeviceType> gradient; // (C,F,P,D), allocated/initialized in main below

Kokkos::DynRankView<double, DeviceType> jacobian;
Kokkos::DynRankView<double, DeviceType> jacobian_inv;
Kokkos::DynRankView<double, DeviceType> tempGrad;

shards::CellTopology cell_topology = shards::getCellTopologyData<shards::Hexahedron<8> >();

using BasisFamily = Intrepid2::DerivedNodalBasisFamily<DeviceType>;
BasisFamily::BasisPtr basis;

/** \brief Create a uniform Cartesian mesh, with origin at 0, and domain extents and mesh widths that can be different in different coordinate dimensions.
   \param [in] domainExtents - array specifying the extent of the domain in each coordinate dimension.
   \param [in] gridCellCounts - array specifying the number of cells in each coordinate dimension.
   \return a uniform Cartesion mesh, with origin at 0, with the specified domain extents and grid cell counts.
*/
template<class PointScalar, int spaceDim, typename DeviceType>
inline Intrepid2::CellGeometry<PointScalar,spaceDim,DeviceType> uniformCartesianMesh(const Kokkos::Array<PointScalar,spaceDim> &domainExtents,
                                                                                     const Kokkos::Array<int,spaceDim> &gridCellCounts)
{
  Kokkos::Array<PointScalar,spaceDim> origin;
  for (int d=0; d<spaceDim; d++)
  {
    origin[d] = 0.0;
  }
  
  using CellGeometry = ::Intrepid2::CellGeometry<PointScalar, spaceDim, DeviceType>;
  const auto NO_SUBDIVISION = CellGeometry::NO_SUBDIVISION;
  const auto HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS = CellGeometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
  
  return CellGeometry(origin, domainExtents, gridCellCounts, NO_SUBDIVISION, HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS);
}

void intrepid2_gradient_operator(const int n, const double* coords, double* grad, double* det)
{
  using FunctionSpaceTools = Intrepid2::FunctionSpaceTools<DeviceType>;
  using CellToolsReal2 = Intrepid2::CellTools<DeviceType>;
 
  using FieldContainer = Kokkos::DynRankView<double, DeviceType>;
  
  const FieldContainer coordVec(const_cast<double*>(coords), n, numNodes, numDims);
 
  FieldContainer jacobian_det(det, n, numIntg);
 
  // compute the jacobian, its inverse, and its determinant
  CellToolsReal2::setJacobian(jacobian, points, coordVec, basis);
  CellToolsReal2::setJacobianInv(jacobian_inv, jacobian);
  CellToolsReal2::setJacobianDet(jacobian_det, jacobian);
 
  // compute the Basis Function Gradients
  FunctionSpaceTools::HGRADtransformGRAD(tempGrad, jacobian_inv, gradient);
 
  //
  //  Transform gradient operator from the intrepid order: ( Element, Node Number, Intg Point,  Component)
  //  to index order consistent with the rest of Sierra:   ( Element, Intg Point,  Node Number, Component)
  //
  FieldContainer gradOpContainer(grad, n, numIntg, numNodes, numDims);
  FieldContainer tempGradLocal = tempGrad; // local View to avoid CUDA complaints.
  Kokkos::parallel_for("reorder gradients", n,
  KOKKOS_LAMBDA(const int &ielem)
  {
    for (int ip = 0; ip < numIntg; ++ip) {
      for (int inode = 0; inode < numNodes; ++inode) {
        for (int d = 0; d<numDims; ++d) {
          gradOpContainer(ielem, ip, inode, d) = tempGradLocal(ielem, inode, ip, d);
        }
      }
    }
  });
}

int main(int argc, char *argv[])
{
  int numInvocations = 10000;
  
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  
  std::string test_name = "Intrepid2 Sierra Test";
  
  {
    // Note that the dtor for GlobalMPISession will call Kokkos::finalize_all() but ctor does not call Kokkos::initialize()...
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize(argc,argv);
    Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
    
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
    auto comm = Teuchos::DefaultComm<int>::getComm();
    
    {
      Teuchos::CommandLineProcessor cmdp(false,true); // false: don't throw exceptions; true: do return errors for unrecognized options

      cmdp.setOption("numInvocations", &numInvocations, "How many times to run the main loop");
      cmdp.setOption("worksetSize", &worksetSize, "The size of the workset - how many cells to process at a time");
      
      if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      {
        return -1;
      }
    }
    
    {
      stacked_timer = rcp(new Teuchos::StackedTimer(test_name.c_str()));
      Teuchos::RCP<Teuchos::FancyOStream> verbose_out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
      verbose_out->setShowProcRank(true);
      stacked_timer->setVerboseOstream(verbose_out);
      Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
    }
    
    Kokkos::DynRankView<double, DeviceType> coords, det;
    {
      Teuchos::TimeMonitor setupTimer(*Teuchos::TimeMonitor::getNewTimer(std::string("allocations and other one-time setup")));
      
      gradient = Kokkos::DynRankView<double, DeviceType>("gradient", worksetSize, numFields, numIntg, numDims);
      points   = Kokkos::DynRankView<double, DeviceType>("points", numIntg, numDims); // (P,D)
      
      jacobian     = Kokkos::DynRankView<double, DeviceType>("jacobian", worksetSize, numIntg, numDims, numDims);
      jacobian_inv = Kokkos::DynRankView<double, DeviceType>("jacobian_inv", worksetSize, numIntg, numDims, numDims);
      tempGrad     = Kokkos::DynRankView<double, DeviceType>("tempGrad", worksetSize, numNodes, numIntg, numDims);
      
      basis = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_C1_FEM<DeviceType>() );
      
      const int cellCount = 25000;
      
      Kokkos::Array<double,numDims> domainExtents{1.,1.,1.};
      Kokkos::Array<int,numDims> gridCellCounts{50,50,10}; // 25000 total cells
      
      auto geometry = uniformCartesianMesh<double, numDims, DeviceType>(domainExtents, gridCellCounts);
      
      det    = Kokkos::DynRankView<double, DeviceType>("det", worksetSize, numIntg);
      coords = Kokkos::DynRankView<double, DeviceType>("coords", cellCount, numNodes, numDims);
      Kokkos::parallel_for("initialize coords", cellCount,
      KOKKOS_LAMBDA(const int &cellOrdinal)
      {
        for (int nodeOrdinal=0; nodeOrdinal<numNodes; nodeOrdinal++)
        {
          for (int d=0; d<numDims; d++)
          {
            coords(cellOrdinal,nodeOrdinal,d) = geometry(cellOrdinal,nodeOrdinal,d);
          }
        }
      });
      Kokkos::fence();
    }
    
    {
      Teuchos::TimeMonitor intrepid2Timer(*Teuchos::TimeMonitor::getNewTimer(std::string("intrepid2_gradient_operator")));
      for (int i=0; i<numInvocations; i++)
      {
        // this is only approximately what will be happening in Sierra; there, the pointers here will move from one call to the next
        // corresponding to the cell ordinal.  But we start with this; probably that difference will not be significant.
        intrepid2_gradient_operator(worksetSize, coords.data(), gradient.data(), det.data());
      }
    }
    
    {
      Teuchos::TimeMonitor dtorTimer(*Teuchos::TimeMonitor::getNewTimer(std::string("post-run cleanup")));
      gradient = Kokkos::DynRankView<double, DeviceType>();
      points   = Kokkos::DynRankView<double, DeviceType>();
      basis = Teuchos::null;
      jacobian     = Kokkos::DynRankView<double, DeviceType>();
      jacobian_inv = Kokkos::DynRankView<double, DeviceType>();
      tempGrad     = Kokkos::DynRankView<double, DeviceType>();
    }
    
    stacked_timer->stop(test_name);
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stacked_timer->report(*out, comm, options);
    
    auto xmlOut = stacked_timer->reportWatchrXML(test_name + ' ' + std::to_string(comm->getSize()) + " ranks", comm);
    if(xmlOut.length())
      std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
  }
  
  return 0;
}

