// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MasterElementIntrepid_h
#define Akri_MasterElementIntrepid_h

#include <memory>
#include <stk_topology/topology.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Intrepid2_Basis.hpp>
#include <Intrepid2_Cubature.hpp>

namespace krino {

class MasterElementIntrepid
{
public:
  using ExecutionSpace = Kokkos::DefaultHostExecutionSpace;
  using OutputType = double;
  using PointType = double;
  using IntrepidBasis = Intrepid2::Basis<ExecutionSpace, OutputType, PointType>;
  using IntrepidCubature = Intrepid2::Cubature<ExecutionSpace, OutputType, PointType>;

  MasterElementIntrepid(stk::topology topology);

  // Copy and assignment are not allowed
  MasterElementIntrepid( const MasterElementIntrepid & ) = delete;
  MasterElementIntrepid & operator=( const MasterElementIntrepid & ) = delete;

  stk::topology get_topology() const { return m_topology; }
  unsigned topology_dimension() const { return m_topology.dimension(); }

  // returns the number of integration points
  unsigned num_intg_pts() const { return m_numIntgPts; }

  // returns the number of nodes
  unsigned num_nodes() const { return m_numNodes; }

  //: Query the integration weights
  const double* intg_weights() const { return m_refWeights.data(); }

  //: Query the integration points/stations
  const double* intg_pt_locations() const { return m_refPoints.data(); }

  const double * nodal_parametric_coordinates() const { return m_refCoords.data(); }

  double parametric_volume() const { return m_refVolume; }

  void determinant(
     const int numCoordDims,
     const int nelem,
     const double* coords,  // (numCoordDims,npe,nelem)
     double* det_J,         // (nelem,nint)
     double* error ) const; // (nelem)
  void determinant(
      const int numCoordDims,
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (m_numElemDims,npe_g,nint)
      const int   nelem,
      const double* coords,   // (numCoordDims,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error ) const ; // (nelem)

  //: Returns the values of the nodal interpolation shape functions
  //: at the integration stations.
  const double* shape_fcn() const { return m_shapeFuncs.data(); }
  void shape_fcn(
      const int nint,  // returns array(npe,nint)
      const double* p_coords,
      double* result) const;

  //: Returns the derivatives of the nodal interpolation shape functions
  //: with respect to the local parametric coordinates at the integration
  //: stations.
  const double* shape_fcn_deriv() const { return m_pointGrads.data(); }
  void shape_fcn_deriv(
      const int nint,
      const double* p_coords,
      double* result ) const;

  void interpolate_point(
      const int  npar_coord,
      const double * par_coord,      // (npar_coord)
      const int  ncomp_field,
      const double * field,          // (num_nodes,ncomp_field)
      double * result ) const;       // (ncomp_field)

  void gradient_operator(
      const int numCoordDims,
      const int nint,
      const int npe_g,
      const double* deriv_g,   // (numElemDims,npe_g,nint)
      const int npe_f,
      const double* deriv_f,   // (numElemDims,npe_f,nint)
      const int   nelem,
      const double* coords,    // (numElemDims,npe,nelem)
      double* gradop,          // (numElemDims,npe,nelem,nint)
      double* det_J,           // (nelem,nint)
      double* error) const;

  void scalar_gradient(
      const int   nelem,      //:  number of elements to process
      const double* gradop,     //: (nvec,npe,nelem,nint)
      const double* det_J,      //: (nelem,nint)
      const double* sfield,     //: (npe,nelem)
      double* vector ) const ;  //: (nvec,nelem,nint)
  void scalar_gradient(
      const int   nint,       //:  number of intg points
      const int   nelem,      //:  number of elements to process
      const double* gradop,     //: (nvec,npe,nelem,nint)
      const double* det_J,      //: (nelem,nint)
      const double* sfield,     //: (npe,nelem)
      double* vector ) const ;  //: (nvec,nelem,nint)

private:
  stk::topology m_topology;
  int m_numNodes;
  int m_numElemDims;
  int m_numIntgPts;

  std::unique_ptr<IntrepidBasis> m_intrepidBasis;

  double m_refVolume;

  // Local FieldContainers
  std::vector<double> m_shapeFuncs;
  std::vector<double> m_pointGrads;
  std::vector<double> m_refPoints;
  std::vector<double> m_refWeights;
  std::vector<double> m_refCoords;
};

}  // end namespace krino

#endif // Akri_MasterElementIntrepid_h
