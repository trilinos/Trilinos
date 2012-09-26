#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include "JacobianUtil.hpp"

#include <MeanRatioFunctions.hpp>

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;

    // Prism and Hex element descriptions
    static const int locs_prism[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
                                         {2, 0, 1, 5}, {3, 5, 4, 0},
                                         {4, 3, 5, 1}, {5, 4, 3, 2}};
    static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
                                       {2, 3, 1, 6}, {3, 0, 2, 7},
                                       {4, 7, 5, 0}, {5, 4, 6, 1},
                                       {6, 5, 7, 2}, {7, 6, 4, 3}};

    static inline void set2d(JacobianUtil::Vec3D v, const double *x)
    {
      v[0] = x[0];
      v[1] = x[1];
      v[2] = 0.0;
    }
    static inline void set3d(JacobianUtil::Vec3D v, const double *x)
    {
      v[0] = x[0];
      v[1] = x[1];
      v[2] = x[2];
    }

    static inline Vector3D vector_3D(double *x)
    {
      return Vector3D(x[0], x[1], x[2]);
    }
    static inline Vector3D vector_2D(double *x)
    {
      return Vector3D(x[0], x[1], 0.0);
    }

    void scale_to_unit(MsqMatrix<3,3>& A)
    {
      for (int jvert=0; jvert < 3; jvert++)
        {
          double sum=0.0;
          for (int ixyz=0; ixyz < 3; ixyz++)
            {
              sum += A(ixyz, jvert)*A(ixyz, jvert);
            }
          sum = std::max(1.e-10, std::sqrt(sum));
          for (int ixyz=0; ixyz < 3; ixyz++)
            {
              A(ixyz, jvert) /= sum;
            }
        }
    }

    bool JacobianUtil::jacobian_matrix_2D(double &detJ, MsqMatrix<3,3>& A, const JacobianUtil::Vec3D *x)
    {
      /* Calculate A */
      // x_xi, x_eta, x_zeta => A(ixyz, ixietazeta) = dx_i/dxi_j
      A(0,0) = (x[1][0] - x[0][0]);  
      A(0,1) = (x[2][0] - x[0][0]);
      A(0,2) = 0;

      A(1,0) = (x[1][1] - x[0][1]);
      A(1,1) = (x[2][1] - x[0][1]);
      A(1,2) = 0;

      A(2,0) = (x[1][2] - x[0][2]);
      A(2,1) = (x[2][2] - x[0][2]);
      A(2,2) = 1.0;
      if (m_scale_to_unit) scale_to_unit(A);

      detJ = det(A);
      return detJ < MSQ_MIN;
    }

    bool JacobianUtil::jacobian_matrix_3D(double &detJ, MsqMatrix<3,3>& A, const JacobianUtil::Vec3D *x)
    {
      A(0,0) = (x[1][0] - x[0][0]);  
      A(0,1) = (x[2][0] - x[0][0]);
      A(0,2) = (x[3][0] - x[0][0]);

      A(1,0) = (x[1][1] - x[0][1]);
      A(1,1) = (x[2][1] - x[0][1]);
      A(1,2) = (x[3][1] - x[0][1]);

      A(2,0) = (x[1][2] - x[0][2]);
      A(2,1) = (x[2][2] - x[0][2]);
      A(2,2) = (x[3][2] - x[0][2]);
      if (m_scale_to_unit) scale_to_unit(A);

      detJ = det(A);
      return detJ < MSQ_MIN;
    }

    /// calculates dMetric_dA (@param dMdA) times dA/dx_n_i to get the gradient of the metric dMetric_dx_n_i, where
    /// x_n_i is the i'th coordinate of the n'th node in the element.
    /// @param indices are the 4 indices associated with this corner of the element, passed in from the grad_metric function
    /// Note: this is the dMetric_dx_n_i term associated with a particular corner, so there are up to nnode of these passed
    /// in from the grad_metric function
    /// @see jacobian_matrix_3D
    void JacobianUtil::grad_util_2d(const MsqMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < spd; j++)
          grad[i][j]=0.0;

      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      grad[indices[2]][0] += dMdA(0,1)*(+1); grad[indices[0]][0] += dMdA(0,1)*(-1);
      //grad[indices[3]][0] += dMdA(0,2)*(+1); grad[indices[0]][0] += dMdA(0,2)*(-1);

      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      grad[indices[2]][1] += dMdA(1,1)*(+1); grad[indices[0]][1] += dMdA(1,1)*(-1);
      //grad[indices[3]][1] += dMdA(1,2)*(+1); grad[indices[0]][1] += dMdA(1,2)*(-1);

      //grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      //grad[indices[2]][2] += dMdA(2,1)*(+1); grad[indices[0]][2] += dMdA(2,1)*(-1);
      //grad[indices[3]][2] += dMdA(2,2)*(+1); grad[indices[0]][2] += dMdA(2,2)*(-1);

    }

    void JacobianUtil::grad_util(const MsqMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < spd; j++)
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

    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate()
    bool JacobianUtil::operator()(double& m,  PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field,
                                  const CellTopologyData * topology_data )
    {
      MsqError err;
      static MsqMatrix<3,3> J;

      //static Vector3D n;			// Surface normal for 2D objects

      int i=0;

      //static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = stk::percept::PerceptMesh::get_cell_topology(element);

      stk::mesh::PairIterRelation v_i = element.relations(eMesh.node_rank());
      m_num_nodes = v_i.size();

      //#define VERTEX_2D(vi) vector_2D( eMesh.field_data(coord_field, *vi.entity() ) )
      //#define VERTEX_3D(vi) vector_3D( eMesh.field_data(coord_field, *vi.entity() ) )

#define VERTEX(vi)  stk::mesh::field_data( *static_cast<const VectorFieldType *>(coord_field) , *vi.entity() )

      switch(topology_data->key) 
        {
        case shards::Triangle<3>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          set2d(mCoords[0], VERTEX(v_i[0]));
          set2d(mCoords[1], VERTEX(v_i[1]));
          set2d(mCoords[2], VERTEX(v_i[2]));
          metric_valid = jacobian_matrix_2D(m, J, mCoords);
          for (i = 0; i < 3; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;
    
        case shards::Quadrilateral<4>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            set2d(mCoords[0], VERTEX(v_i[locs_hex[i][0]]));
            set2d(mCoords[1], VERTEX(v_i[locs_hex[i][1]]));
            set2d(mCoords[2], VERTEX(v_i[locs_hex[i][2]]));
            metric_valid = jacobian_matrix_2D(m_detJ[i], m_J[i], mCoords);
          }
          m = average_metrics(m_detJ, 4, err); MSQ_ERRZERO(err);
          break;

        case shards::Tetrahedron<4>::key:
          set3d(mCoords[0], VERTEX(v_i[0]));
          set3d(mCoords[1], VERTEX(v_i[1]));
          set3d(mCoords[2], VERTEX(v_i[2]));
          set3d(mCoords[3], VERTEX(v_i[3]));
          metric_valid = jacobian_matrix_3D(m, J, mCoords);
          for (i = 0; i < 4; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;

        case shards::Pyramid<5>::key:
          for (i = 0; i < 4; ++i) {
            set3d(mCoords[0], VERTEX(v_i[ i     ]));
            set3d(mCoords[1], VERTEX(v_i[(i+1)%4]));
            set3d(mCoords[2], VERTEX(v_i[(i+3)%4]));
            set3d(mCoords[3], VERTEX(v_i[ 4     ]));
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords);
          }
          // FIXME
          m_J[4] = (m_J[0]+m_J[1]+m_J[2]+m_J[3]);
          m_J[4] *= 0.25;
          m_detJ[4] = det(m_J[4]);
          m = average_metrics(m_detJ, 5, err); MSQ_ERRZERO(err);
          break;

        case shards::Wedge<6>::key:
          for (i = 0; i < 6; ++i) {
            set3d(mCoords[0], VERTEX(v_i[locs_prism[i][0]]));
            set3d(mCoords[1], VERTEX(v_i[locs_prism[i][1]]));
            set3d(mCoords[2], VERTEX(v_i[locs_prism[i][2]]));
            set3d(mCoords[3], VERTEX(v_i[locs_prism[i][3]]));
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords);
          }
          m = average_metrics(m_detJ, 6, err); MSQ_ERRZERO(err);
          break;

        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            set3d(mCoords[0], VERTEX(v_i[locs_hex[i][0]]));
            set3d(mCoords[1], VERTEX(v_i[locs_hex[i][1]]));
            set3d(mCoords[2], VERTEX(v_i[locs_hex[i][2]]));
            set3d(mCoords[3], VERTEX(v_i[locs_hex[i][3]]));
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords);
          }
          m = average_metrics(m_detJ, 8, err); MSQ_ERRZERO(err);
          break;


          // unimplemented
        case shards::Node::key:
        case shards::Particle::key:
        case shards::Line<2>::key:
        case shards::Line<3>::key:
        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:
      
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:
      
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
      
        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<10>::key:
        case shards::Tetrahedron<11>::key:
      
        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
      
        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
      
        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
      
        case shards::Pentagon<5>::key:
        case shards::Hexagon<6>::key:

        default:
          shards::CellTopology topology(topology_data);
          std::cout << "topology = " << topology.getName() << std::endl;
          throw std::runtime_error("unknown/unhandled topology in JacobianUtil");
          break;

        } // end switch over element type

      return metric_valid;
    }

    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate(), and TargetMetricUtil
    /// fills the mGrad member variable given the array of (member variable) m_dMetric_dA terms
    bool JacobianUtil::grad_metric_util( PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field,
                                        const CellTopologyData * topology_data )
    {
      MsqError err;
      static MsqMatrix<3,3> J;

      //static Vector3D n;			// Surface normal for 2D objects

      int i=0;

      //static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = stk::percept::PerceptMesh::get_cell_topology(element);

      stk::mesh::PairIterRelation v_i = element.relations(eMesh.node_rank());
      m_num_nodes = v_i.size();

      //#define VERTEX_2D(vi) vector_2D( eMesh.field_data(coord_field, *vi.entity() ) )
      //#define VERTEX_3D(vi) vector_3D( eMesh.field_data(coord_field, *vi.entity() ) )

#define VERTEX(vi)  stk::mesh::field_data( *static_cast<const VectorFieldType *>(coord_field) , *vi.entity() )

      const int indices_tri[3] = {0,1,2};
      const int indices_tet[4] = {0,1,2,3};
      switch(topology_data->key) 
        {
        case shards::Triangle<3>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 3; i++) { 
            //void JacobianUtil::grad_util(const MsqMatrix<3,3>& dMetric_dA, double grad[NNODES_MAX][3], int nnode, int spd, int *indices, int nind)
            grad_util_2d(m_dMetric_dA[i], m_grad[i], 3, 2, indices_tri, 3);
          }
          break;
    
        case shards::Quadrilateral<4>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            const int *indices_quad = locs_hex[i];
            grad_util_2d(m_dMetric_dA[i], m_grad[i], 4, 2, indices_quad, 3);
          }
          break;

        case shards::Tetrahedron<4>::key:
          for (i = 0; i < 4; ++i) {
            grad_util(m_dMetric_dA[i], m_grad[i], 4, 3, indices_tet, 4);
          }
          break;

        case shards::Pyramid<5>::key:
          for (i = 0; i < 4; ++i) {
            const int indices_pyr[4] = {i, (i+1)%4, (i+3)%4, 4};
            grad_util(m_dMetric_dA[i], m_grad[i], 5, 3, indices_pyr, 4);
          }
          // FIXME
          for ( i=0; i < 5; i++) 
            {
              for (int j=0; j < 3; j++)
                m_grad[4][i][j] = 0.0;
            }
          for (int k=0; k < 4; k++) 
            {
              for (i=0; i < 5; i++) 
                {
                  for (int j=0; j < 3; j++)
                    m_grad[4][i][j] += m_grad[k][i][j]*0.25;
                }
            }
          break;

        case shards::Wedge<6>::key:
          for (i = 0; i < 6; ++i) {
            const int *indices_prism = locs_prism[i];
            grad_util(m_dMetric_dA[i], m_grad[i], 6, 3, indices_prism, 4);
          }
          break;

        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            const int *indices_hex = locs_hex[i];
            grad_util(m_dMetric_dA[i], m_grad[i], 8, 3, indices_hex, 4);
          }
          break;


          // unimplemented
        case shards::Node::key:
        case shards::Particle::key:
        case shards::Line<2>::key:
        case shards::Line<3>::key:
        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:
      
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:
      
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
      
        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<10>::key:
        case shards::Tetrahedron<11>::key:
      
        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
      
        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
      
        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
      
        case shards::Pentagon<5>::key:
        case shards::Hexagon<6>::key:

        default:
          shards::CellTopology topology(topology_data);
          std::cout << "topology = " << topology.getName() << std::endl;
          throw std::runtime_error("unknown/unhandled topology in JacobianUtil");
          break;

        } // end switch over element type

      return metric_valid;
    }

  }
}

#endif
