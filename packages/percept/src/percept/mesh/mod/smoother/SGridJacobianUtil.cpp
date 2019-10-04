// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>

#if !defined(NO_GEOM_SUPPORT)

#include "SGridJacobianUtil.hpp"
#include <percept/structured/BlockStructuredGrid.hpp>

namespace percept {

  template<size_t N>
  std::ostream& operator<<(std::ostream& os, const std::array<double, N>& arr)
  {
    for (unsigned i=0; i < N; ++i)
      os << arr[i] << (i == N-1 ? "" : " ");
    return os;
  }

  /**
   *   Structured grid nodes numbering
   *
   *         6                      7
   *           o------------------o
   *          /|                 /|
   *         / |                / |
   *        /  |               /  |
   *       /   |              /   |
   *      /    |             /    |
   *     /     |            /     |
   *  4 /      |         5 /      |
   *   o------------------o       |
   *   |       |          |       |
   *   |   2   o----------|-------o  3
   *   |      /           |      /
   *   |     /            |     /
   *   |    /             |    /
   *   |   /              |   /
   *   |  /               |  /
   *   | /                | /
   *   |/                 |/
   *   o------------------o
   *  0                    1
   *
   */

  static const int locs_hex[8][4] = {{0, 1, 2, 4}, {1, 3, 0, 5},
                                     {2, 0, 3, 6}, {3, 2, 1, 7},
                                     {4, 6, 5, 0}, {5, 4, 7, 1},
                                     {6, 7, 4, 2}, {7, 5, 6, 3}};

  /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate()

  bool JacobianUtilImpl<StructuredGrid>::
  operator()(double& averageJ, PerceptMesh& eMesh, typename StructuredGrid::MTElement cell_ijkb, typename StructuredGrid::MTField *coord_field,
             const typename StructuredGrid::MTCellTopology * topology_data_in )
  {
    EXCEPTWATCH;

    DenseMatrix<3,3> J;
    int i=0;

    bool metric_invalid = false;

    using Array4D = typename StructuredGrid::MTField::Array4D;
    std::shared_ptr<BlockStructuredGrid> bgrid = eMesh.get_block_structured_grid();
    unsigned iblock = cell_ijkb[3];
    VERIFY_OP_ON(iblock, <, bgrid->m_sblocks.size(), "bad iblock");
    std::shared_ptr<StructuredBlock> sgrid = bgrid->m_sblocks[iblock];
    Array4D& a4d = (*coord_field->m_block_fields[iblock]);

    const int A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];

    m_num_nodes = 8; //v_i.size();
    std::vector<std::array<double,3> > v_i(8);
    std::array<unsigned,3> indx{{0,0,0}}, II{{0,0,0}};

    unsigned cnt = 0;
    for ( indx[2]=0; indx[2] < 2; ++indx[2])
      {
        II[2] = indx[2] + cell_ijkb[2];
        for ( indx[1]=0; indx[1] < 2; ++indx[1])
          {
            II[1] = indx[1] + cell_ijkb[1];
            for ( indx[0]=0; indx[0] < 2; ++indx[0])
              {
                II[0] = indx[0] + cell_ijkb[0];
                for (unsigned ic=0; ic < 3; ++ic)
                  {
                    v_i[cnt][ic] = a4d(II[A0], II[A1], II[A2], ic);
                  }
                ++cnt;
              }
          }
      }

#define VERTEX(vi)  (vi.data())

    for (i = 0; i < 8; ++i) {
      bool mi = jacobian_matrix_3D(m_detJ[i], m_J[i],
                                   VERTEX(v_i[locs_hex[i][0]]),
                                   VERTEX(v_i[locs_hex[i][1]]),
                                   VERTEX(v_i[locs_hex[i][2]]),
                                   VERTEX(v_i[locs_hex[i][3]]));

      metric_invalid = metric_invalid || mi;
    }

    averageJ = average_metrics(m_detJ, 8);

    return metric_invalid;
  }

  /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate(), and TargetMetricUtil
  /// fills the mGrad member variable given the array of (member variable) m_dMetric_dA terms

  bool JacobianUtilImpl<StructuredGrid>::grad_metric_util(  PerceptMesh& eMesh, typename StructuredGrid::MTElement element, typename StructuredGrid::MTField *coord_field,
                                                            const typename StructuredGrid::MTCellTopology * topology_data )
  {
    int i=0;
    bool metric_valid = true;

    for (i = 0; i < 8; ++i) {
      const int *indices_hex = locs_hex[i];
      grad_util(m_dMetric_dA[i], m_grad[i], 8, 3, indices_hex, 4);
    }

    return metric_valid;
  }


  /// calculates dMetric_dA (@param dMdA) times dA/dx_n_i to get the gradient of the metric dMetric_dx_n_i, where
  /// x_n_i is the i'th coordinate of the n'th node in the element.
  /// @param indices are the 4 indices associated with this corner of the element, passed in from the grad_metric function
  /// Note: this is the dMetric_dx_n_i term associated with a particular corner, so there are up to nnode of these passed
  /// in from the grad_metric function
  /// @see jacobian_matrix_3D


  void JacobianUtilImpl<StructuredGrid>::grad_util(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind)
  {
    for (int i=0; i < nnode; i++)
      for (int j=0; j < 3; j++)
        grad[i][j]=0.0;

    grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
    grad[indices[2]][0] += dMdA(0,1)*(+1); grad[indices[0]][0] += dMdA(0,1)*(-1);
    grad[indices[3]][0] += dMdA(0,2)*(+1); grad[indices[0]][0] += dMdA(0,2)*(-1);

    grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
    grad[indices[2]][1] += dMdA(1,1)*(+1); grad[indices[0]][1] += dMdA(1,1)*(-1);
    grad[indices[3]][1] += dMdA(1,2)*(+1); grad[indices[0]][1] += dMdA(1,2)*(-1);

    grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
    grad[indices[2]][2] += dMdA(2,1)*(+1); grad[indices[0]][2] += dMdA(2,1)*(-1);
    grad[indices[3]][2] += dMdA(2,2)*(+1); grad[indices[0]][2] += dMdA(2,2)*(-1);

  }


  bool SGridJacobianUtil(double& averageJ, double detJ[8], typename StructuredGrid::MTField::Array4D coords,
		  typename StructuredGrid::MTElement cell_ijkb)
  {
	  DenseMatrix<3,3> J;

	  bool metric_invalid = false;

	  //const int A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];
	  const int A0 = 0, A1 = 1, A2 = 2;

	  std::vector<std::array<double,3> > v_i(8);
	  std::array<unsigned,3> indx{{0,0,0}}, II{{0,0,0}};

	  unsigned cnt = 0;
	  for ( indx[2]=0; indx[2] < 2; ++indx[2])
	  {
		  II[2] = indx[2] + cell_ijkb[2];
		  for ( indx[1]=0; indx[1] < 2; ++indx[1])
		  {
			  II[1] = indx[1] + cell_ijkb[1];
			  for ( indx[0]=0; indx[0] < 2; ++indx[0])
			  {
				  II[0] = indx[0] + cell_ijkb[0];
				  for (unsigned ic=0; ic < 3; ++ic)
				  {
					  v_i[cnt][ic] = coords(II[A0], II[A1], II[A2], ic);
				  }
				  ++cnt;
			  }
		  }
	  }

#define VERTEX(vi)  (vi.data())

	  //double x0[3], x1[3], x2[3], x3[3];

		  for (int i = 0; i < 8; ++i) {
			  bool mi = jacobian_matrix_3D(detJ[i], J,
					  VERTEX(v_i[locs_hex[i][0]]),
					  VERTEX(v_i[locs_hex[i][1]]),
					  VERTEX(v_i[locs_hex[i][2]]),
					  VERTEX(v_i[locs_hex[i][3]]));

			  metric_invalid = metric_invalid || mi;
		  }

		  //averageJ = average_metrics(m_detJ, 8);
		  //double average_metrics(const double detJ[NNODES_MAX], const unsigned n)
		  {
			  averageJ=0.0;
			  for (unsigned i=0; i < 8; i++) averageJ += detJ[i];
			  averageJ/=double(8);
		  }

		  return metric_invalid;
  }
} // namespace percept

#endif
