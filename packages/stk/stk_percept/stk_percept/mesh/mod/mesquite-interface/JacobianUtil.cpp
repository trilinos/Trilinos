#if !defined(__IBMCPP__) && defined(STK_BUILT_IN_SIERRA)

#include "JacobianUtil.hpp"

#include <MeanRatioFunctions.hpp>

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;

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

    bool JacobianUtil::jacobian_matrix_2D(double &detJ, MsqMatrix<3,3>& A, const Vector3D *x, const Vector3D &n, const Vector3D &d)
    {
      /* Calculate A */
      // x_xi, x_eta, x_zeta => A(ixyz, ixietazeta) = dx_i/dxi_j
      A(0,0) = d[0]*(x[1][0] - x[0][0]);  
      A(0,1) = d[1]*(x[2][0] - x[0][0]);
      A(0,2) = n[0];

      A(1,0) = d[0]*(x[1][1] - x[0][1]);
      A(1,1) = d[1]*(x[2][1] - x[0][1]);
      A(1,2) = n[1];

      A(2,0) = d[0]*(x[1][2] - x[0][2]);
      A(2,1) = d[1]*(x[2][2] - x[0][2]);
      A(2,2) = n[2];
      if (m_scale_to_unit) scale_to_unit(A);

      detJ = det(A);
      return detJ < MSQ_MIN;
    }

    bool JacobianUtil::jacobian_matrix_3D(double &detJ, MsqMatrix<3,3>& A, const Vector3D *x, const Vector3D &n, const Vector3D &d)
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

    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate()
    bool JacobianUtil::operator()(double& m,  PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field,
                                  const CellTopologyData * topology_data )
    {
      MsqError err;
      static MsqMatrix<3,3> J;

      static Vector3D n;			// Surface normal for 2D objects

      // Prism and Hex element descriptions
      static const int locs_prism[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
                                           {2, 0, 1, 5}, {3, 5, 4, 0},
                                           {4, 3, 5, 1}, {5, 4, 3, 2}};
      static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
                                         {2, 3, 1, 6}, {3, 0, 2, 7},
                                         {4, 7, 5, 0}, {5, 4, 6, 1},
                                         {6, 5, 7, 2}, {7, 6, 4, 3}};
      int i=0;

      static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = stk::percept::PerceptMesh::get_cell_topology(element);

      stk::mesh::PairIterRelation v_i = element.relations(eMesh.node_rank());
      m_num_nodes = v_i.size();

      //#define VERTEX_2D(vi) vector_2D( eMesh.field_data(coord_field, *vi.entity() ) )
      //#define VERTEX_3D(vi) vector_3D( eMesh.field_data(coord_field, *vi.entity() ) )

#define VERTEX_2D(vi) vector_2D( stk::mesh::field_data( *static_cast<const VectorFieldType *>(coord_field) , *vi.entity() ) )
#define VERTEX_3D(vi) vector_3D( stk::mesh::field_data( *static_cast<const VectorFieldType *>(coord_field) , *vi.entity() ) )


      switch(topology_data->key) 
        {
        case shards::Triangle<3>::key:
          n[0] = 0; n[1] = 0; n[2] = 1;
          mCoords[0] = VERTEX_2D(v_i[0]);
          mCoords[1] = VERTEX_2D(v_i[1]);
          mCoords[2] = VERTEX_2D(v_i[2]);
          metric_valid = jacobian_matrix_2D(m, J, mCoords, n, d_con);
          for (i = 0; i < 4; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;
    
        case shards::Quadrilateral<4>::key:
          n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            mCoords[0] = VERTEX_2D(v_i[locs_hex[i][0]]);
            mCoords[1] = VERTEX_2D(v_i[locs_hex[i][1]]);
            mCoords[2] = VERTEX_2D(v_i[locs_hex[i][2]]);
            metric_valid = jacobian_matrix_2D(m_detJ[i], m_J[i], mCoords, n, d_con);
          }
          m = average_metrics(m_detJ, 4, err); MSQ_ERRZERO(err);
          break;

        case shards::Tetrahedron<4>::key:
          mCoords[0] = VERTEX_3D(v_i[0]);
          mCoords[1] = VERTEX_3D(v_i[1]);
          mCoords[2] = VERTEX_3D(v_i[2]);
          mCoords[3] = VERTEX_3D(v_i[3]);
          metric_valid = jacobian_matrix_3D(m, J, mCoords, n, d_con);
          for (i = 0; i < 4; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;

        case shards::Pyramid<5>::key:
          for (i = 0; i < 4; ++i) {
            mCoords[0] = VERTEX_3D(v_i[ i     ]);
            mCoords[1] = VERTEX_3D(v_i[(i+1)%4]);
            mCoords[2] = VERTEX_3D(v_i[(i+3)%4]);
            mCoords[3] = VERTEX_3D(v_i[ 4     ]);
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords, n, d_con);
          }
          // FIXME
          m_J[4] = (m_J[0]+m_J[1]+m_J[2]+m_J[3]);
          m_J[4] *= 0.25;
          m_detJ[4] = det(m_J[4]);
          m = average_metrics(m_detJ, 5, err); MSQ_ERRZERO(err);
          break;

        case shards::Wedge<6>::key:
          for (i = 0; i < 6; ++i) {
            mCoords[0] = VERTEX_3D(v_i[locs_prism[i][0]]);
            mCoords[1] = VERTEX_3D(v_i[locs_prism[i][1]]);
            mCoords[2] = VERTEX_3D(v_i[locs_prism[i][2]]);
            mCoords[3] = VERTEX_3D(v_i[locs_prism[i][3]]);
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords, n, d_con);
          }
          m = average_metrics(m_detJ, 6, err); MSQ_ERRZERO(err);
          break;

        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            mCoords[0] = VERTEX_3D(v_i[locs_hex[i][0]]);
            mCoords[1] = VERTEX_3D(v_i[locs_hex[i][1]]);
            mCoords[2] = VERTEX_3D(v_i[locs_hex[i][2]]);
            mCoords[3] = VERTEX_3D(v_i[locs_hex[i][3]]);
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i], mCoords, n, d_con);
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
  }
}

#endif
