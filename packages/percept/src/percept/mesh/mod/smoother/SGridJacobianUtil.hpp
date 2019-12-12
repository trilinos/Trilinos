// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SGridJacobianUtil_hpp
#define SGridJacobianUtil_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)


#include <percept/PerceptMesh.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/MeshType.hpp>
#include <percept/structured/StructuredCellIndex.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>

namespace percept {

inline bool jacobian_matrix_3D(double &detJ,
		DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
{
  A(0,0) = (x1[0] - x0[0]);
  A(0,1) = (x2[0] - x0[0]);
  A(0,2) = (x3[0] - x0[0]);

  A(1,0) = (x1[1] - x0[1]);
  A(1,1) = (x2[1] - x0[1]);
  A(1,2) = (x3[1] - x0[1]);

  A(2,0) = (x1[2] - x0[2]);
  A(2,1) = (x2[2] - x0[2]);
  A(2,2) = (x3[2] - x0[2]);

  detJ = det(A);

  return detJ < 0.0;
}

// only for 3D structured grids (Hex mesh)
  template<>
  class JacobianUtilImpl<StructuredGrid> : public  JacobianUtilBase<StructuredGrid>
  {

  public:
    using Array4D = typename MTSGridField::Array4D;

    using Base =  JacobianUtilBase<StructuredGrid>;

    using  Base::Vec3D;
    using  Base::NNODES_MAX;
    using  Base::m_detJ;
    using  Base::m_J;
    using  Base::m_dMetric_dA;
    using  Base::m_grad;
    using  Base::m_num_nodes;
    using  Base::m_scale_to_unit;
    using  Base::m_use_approximate_quadratic_jacobian;

    using MeshType = StructuredGrid;

    enum { NELEM_TYPES = 1 };
    //using  Base::NELEM_TYPES;


    JacobianUtilImpl(bool use_approximate_quadratic_jacobian=true) : Base(use_approximate_quadratic_jacobian)
    {
    }

    // evaluate jacobian at each node of the cell given by cell_bijk (block, ijk)
    bool operator()(double& averageJ, PerceptMesh& eMesh, typename StructuredGrid::MTElement element, typename StructuredGrid::MTField *coord_field,
                    const typename StructuredGrid::MTCellTopology * topology_data_in = 0 );


    bool grad_metric_util( PerceptMesh& eMesh, typename MeshType::MTElement element, typename MeshType::MTField *coord_field,
                           const typename MeshType::MTCellTopology * topology_data );



  private:

    inline bool jacobian_matrix_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
    {
      A(0,0) = (x1[0] - x0[0]);
      A(0,1) = (x2[0] - x0[0]);
      A(0,2) = (x3[0] - x0[0]);

      A(1,0) = (x1[1] - x0[1]);
      A(1,1) = (x2[1] - x0[1]);
      A(1,2) = (x3[1] - x0[1]);

      A(2,0) = (x1[2] - x0[2]);
      A(2,1) = (x2[2] - x0[2]);
      A(2,2) = (x3[2] - x0[2]);

      detJ = det(A);

      return detJ < 0.0;
    }

    void grad_util(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);

  };

    // evaluate jacobian at each node of the cell given by cell_bijk (block, ijk)
    bool SGridJacobianUtil(double& averageJ, double detJ[8], typename StructuredGrid::MTField::Array4D coords,
    		typename StructuredGrid::MTElement element);
}

#endif
#endif
