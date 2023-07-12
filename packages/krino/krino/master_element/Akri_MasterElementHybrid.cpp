// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementHybrid.hpp>
#include <Akri_MasterElementBasis.hpp>
#include <Akri_MasterElement.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/MetaData.hpp> // for get_cell_topology
#include <Akri_MasterElementCalc.hpp>

#ifdef __INTEL_COMPILER
#include <Intrepid_FieldContainer.hpp>
#else
//FieldContainer has shadowed variables
//this disables the checking on GCC only
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Intrepid_FieldContainer.hpp>
#pragma GCC diagnostic pop
#endif
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Cubature.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_HGRAD_LINE_Cn_FEM.hpp>

namespace krino {

MasterElementHybrid::MasterElementHybrid(
    stk::topology topology,
    std::unique_ptr<Basis> basis)
: m_topology(topology),
  m_Basis(std::move(basis))
{
  // set the cubature
  Intrepid::DefaultCubatureFactory<double> cubatureFactory;
  Teuchos::RCP<Intrepid::Cubature<double>> intrepidCubature = cubatureFactory.create(stk::mesh::get_cell_topology(topology), 2*m_Basis->degree());
  m_numIntgPts = intrepidCubature->getNumPoints();

  m_numNodes = topology.num_nodes();
  m_numElemDims  = intrepidCubature->getDimension();

  // Allocate reference data
  m_shapeFuncs.resize(m_numIntgPts*m_numNodes);
  m_pointGrads.resize(m_numIntgPts*m_numNodes*m_numElemDims );
  m_refPoints.resize(m_numIntgPts*m_numElemDims);
  m_refWeights.resize(m_numIntgPts);
  m_refCoords.resize(m_numNodes*m_numElemDims);
  m_refVolume = m_Basis->parametric_volume();

  // retrieve the cubature points and weights
  std::vector<int> refPointsDims = {m_numIntgPts, m_numElemDims};
  Intrepid::FieldContainer<double> refPointsFC( refPointsDims, m_refPoints.data() );
  std::vector<int> refWeightsDims = {m_numIntgPts};
  Intrepid::FieldContainer<double> refWeightsFC( refWeightsDims, m_refWeights.data() );
  intrepidCubature->getCubature(refPointsFC, refWeightsFC);

  // compute the reference values and gradients at the integration points
  m_Basis->shape_fcn(m_numIntgPts, m_refPoints.data(), m_shapeFuncs.data());
  m_Basis->shape_fcn_deriv(m_numIntgPts, m_refPoints.data(), m_pointGrads.data());
  m_Basis->nodal_parametric_coordinates(m_refCoords.data());

  m_centroidParCoords.resize(m_numElemDims, 0.);
  for(int n=0; n < m_numNodes; ++n)
  {
    for(int d=0; d < m_numElemDims; ++d)
    {
      m_centroidParCoords[d] += m_refCoords[n*m_numElemDims + d];
    }
  }
}

void
MasterElementHybrid::determinant(
    const int numCoordDims,
    const int nelem,
    const double* coords,  // (numCoordDims,npe,nelem)
    double* det_J,         // (nelem,nint)
    double* error ) const  // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, numCoordDims, m_numIntgPts, m_numNodes, m_pointGrads.data(), nelem, coords, det_J, error);
}

void
MasterElementHybrid::shape_fcn(
    const int nint,  // returns array(npe,nint)
    const double* p_coords,
    double* result) const
{
  m_Basis->shape_fcn(nint, p_coords, result);
}

void
MasterElementHybrid::shape_fcn_deriv(
    const int nint,
    const double* p_coords,
    double* result ) const
{
  m_Basis->shape_fcn_deriv(nint, p_coords, result);
}


void
MasterElementHybrid::interpolate_point(
    const int  npar_coord,
    const double * par_coord,      // (npar_coord)
    const int  ncomp_field,
    const double * field,          // (ncomp_field,num_nodes)
    double * result ) const        // (ncomp_field)
{
  std::vector<double> shape(m_numNodes);
  shape_fcn(1, par_coord, shape.data());
  for ( int comp(0); comp < ncomp_field; ++comp ) {
    result[ comp ] = 0.0;
    for ( int node(0); node < m_numNodes; ++node ) {
      result[ comp ] += shape[node] * field[ ncomp_field * node + comp ];
    }
  }
}

void
MasterElementHybrid::scalar_gradient(
    const int   nelem,      //:  number of elements to process
    const double* gradop,     //: (nvec,npe,nelem,nint)
    const double* det_J,      //: (nelem,nint)
    const double* sfield,     //: (npe,nelem)
    double* vector ) const    //: (nvec,nelem,nint)
{
  MasterElementCalc::scalar_gradient(m_numIntgPts, nelem, m_numElemDims, m_numNodes, gradop, det_J, sfield, vector);
}

void
MasterElementHybrid::scalar_gradient(
        const int   nint,       //:  number of intg points
        const int   nelem,      //:  number of elements to process
        const double* gradop,     //: (nvec,npe,nelem,nint)
        const double* det_J,      //: (nelem,nint)
        const double* sfield,     //: (npe,nelem)
        double* vector ) const    //: (nvec,nelem,nint)
{
  MasterElementCalc::scalar_gradient(nint, nelem, m_numElemDims, m_numNodes, gradop, det_J, sfield, vector);
}

void
MasterElementHybrid::determinant(
      const int numCoordDims,
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (m_numElemDims,npe_g,nint)
      const int   nelem,
      const double* coords,   // (numCoordDims,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error ) const   // (nelem)
{
  MasterElementCalc::determinant(m_numElemDims, numCoordDims, nint, npe_g, deriv_g, nelem, coords, det_J, error);
}

void
MasterElementHybrid::gradient_operator(
        const int numCoordDims,
        const int nint,
        const int npe_g,
        const double* deriv_g,   // (nvec,npe_g,nint)
        const int npe_f,
        const double* deriv_f,   // (nvec,npe_f,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* gradop,          // (nvec,npe,nelem,nint)
        double* det_J,           // (nelem,nint)
        double* error) const
{
  ThrowRequireMsg(m_numElemDims == numCoordDims, "MasterElementHybrid::gradient_operator does not support lower rank elements in higher dimensions (e.g. BAR,QUAD,TRI in 3D).");
  MasterElementCalc::gradient_operator(m_numElemDims, nint, npe_g, deriv_g, npe_f, deriv_f, nelem, coords, gradop, det_J, error);
}

} // namespace krino
