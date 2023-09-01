// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include "VolumeUtil.hpp"

#include "mpi.h"

  namespace percept {

#define MIN_JAC (0.0)

    // Prism and Hex element descriptions
    static const int locs_prism[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
                                         {2, 0, 1, 5}, {3, 5, 4, 0},
                                         {4, 3, 5, 1}, {5, 4, 3, 2}};
    static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
                                       {2, 3, 1, 6}, {3, 0, 2, 7},
                                       {4, 7, 5, 0}, {5, 4, 6, 1},
                                       {6, 5, 7, 2}, {7, 6, 4, 3}};

    bool VolumeUtil::volume_matrix_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
    {
      /* Calculate A */
      // x_xi, x_eta, x_zeta => A(ixyz, ixietazeta) = dx_i/dxi_j
      A(0,0) = (x[1][0] - x[0][0]);
      A(0,1) = (x[2][0] - x[0][0]);
      A(0,2) = 0;

      A(1,0) = (x[1][1] - x[0][1]);
      A(1,1) = (x[2][1] - x[0][1]);
      A(1,2) = 0;

      A(2,0) = 0; // (x[1][2] - x[0][2]);
      A(2,1) = 0; // (x[2][2] - x[0][2]);
      A(2,2) = 1.0;

      detJ = det(A);
      return detJ < 0.0;
    }

    bool VolumeUtil::volume_matrix_2D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
    {
      double x0[3] = {(x[1][0] - x[0][0]), (x[1][1] - x[0][1]), (x[1][2] - x[0][2])};
      double x1[3] = {(x[2][0] - x[0][0]), (x[2][1] - x[0][1]), (x[2][2] - x[0][2])};
      double v[3] = {0,0,0};
      Math::cross_3d(x0, x1, v);

      double vl = Math::norm_3d(v);
      if (vl == 0.0) vl = 1.0;
      A(0,0) = v[0]/vl;
      A(0,1) = v[1]/vl;
      A(0,2) = v[2]/vl;

      A(1,0) = x0[0];
      A(1,1) = x0[1];
      A(1,2) = x0[2];

      A(2,0) = x1[0];
      A(2,1) = x1[1];
      A(2,2) = x1[2];

      detJ =  det(A);
      double td = std::sqrt(Math::dot_3d(v,v));
      VERIFY_OP_ON(std::fabs(td - detJ), < , 1.e-6*(std::fabs(td)+std::fabs(detJ)), "bad det");
      //std::cout << "A= \n" << A << " detJ= " << detJ << std::endl;
      return detJ < 0.0;
    }

    bool VolumeUtil::volume_matrix_1D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
    {
      A(0,0) = 0;
      A(0,1) = 0;
      A(0,2) = 0;

      A(1,0) = 0;
      A(1,1) = 0;
      A(1,2) = 0;

      A(2,0) = 0;
      A(2,1) = 0;
      A(2,2) = 0;

      double x0[3] = {(x[1][0] - x[0][0]), (x[1][1] - x[0][1]), (x[1][2] - x[0][2])};
      detJ =  std::sqrt(Math::dot_3d(x0,x0));
      return false;
    }

    bool VolumeUtil::volume_matrix_1D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
    {
      A(0,0) = 0;
      A(0,1) = 0;
      A(0,2) = 0;

      A(1,0) = 0;
      A(1,1) = 0;
      A(1,2) = 0;

      A(2,0) = 0;
      A(2,1) = 0;
      A(2,2) = 0;

      double x0[3] = {(x[1][0] - x[0][0]), (x[1][1] - x[0][1]), 0.0};
      detJ =  std::sqrt(Math::dot_3d(x0,x0));
      return false;
    }

    void VolumeUtil::check_unhandled_topo(PerceptMesh& eMesh, const CellTopologyData * topology_data)
    {
      if (!m_use_approximate_quadratic_volume)
        {
          shards::CellTopology topology(topology_data);
          std::cout << "topology = " << topology.getName() << std::endl;
          throw std::runtime_error("unknown/unhandled topology in VolumeUtil");
        }
      else
        {
          static int count = 0;
          if (count < 3)
            {
              shards::CellTopology topology(topology_data);
              if (eMesh.get_rank() == 0)
                std::cout << "WARNING[percept::VolumeUtil]: topology [" << topology.getName() << "] is unhandled, and will use linear volume estimate" << std::endl;
              ++count;
            }
        }
    }
    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate()
    bool VolumeUtil::operator()(double& m,  PerceptMesh& eMesh, stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field,
                                  const CellTopologyData * topology_data )
    {
      EXCEPTWATCH;
      static DenseMatrix<3,3> J;
      int spatialDimension = eMesh.get_spatial_dim();

      //static Vector3D n;			// Surface normal for 2D objects

      int i=0;

      //static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = eMesh.get_cell_topology(element);

      const MyPairIterRelation v_i(eMesh, element, eMesh.node_rank() );
      m_num_nodes = v_i.size();

      const double *x2d[3] = {0,0,0};
      //const double *x3d[4] = {0,0,0,0};

#define VERTEX(vi)  static_cast<double*>(stk::mesh::field_data( *coord_field, vi.entity() ))

      switch(topology_data->key)
        {
        case shards::Particle::key:
          m = 1.0;
          return metric_valid;
          break;
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<6>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::ShellTriangle<3>::key:
        case shards::Triangle<3>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          x2d[0] = VERTEX(v_i[0]);
          x2d[1] = VERTEX(v_i[1]);
          x2d[2] = VERTEX(v_i[2]);
          if (spatialDimension == 3)
            metric_valid = volume_matrix_2D_in_3D(m, J, x2d);
          else
            metric_valid = volume_matrix_tri_2D(m, J, x2d);
          for (i = 0; i < 3; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;

        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::ShellQuadrilateral<4>::key:
        case shards::Quadrilateral<4>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            x2d[0] =  VERTEX(v_i[locs_hex[i][0]]);
            x2d[1] =  VERTEX(v_i[locs_hex[i][1]]);
            x2d[2] =  VERTEX(v_i[locs_hex[i][2]]);
            if (spatialDimension == 3)
              metric_valid = volume_matrix_2D_in_3D(m_detJ[i], m_J[i], x2d);
            else
              metric_valid = volume_matrix_2D(m_detJ[i], m_J[i], x2d);
          }
          m = average_metrics(m_detJ, 4);
          break;

        // FIXME - tmp
        case shards::Tetrahedron<10>::key:
          metric_valid = true;
          for (i = 0; i < 10; i++) { m_detJ[i] = 1.0; identity(m_J[i]); }
          break;

        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<11>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Tetrahedron<4>::key:
          metric_valid = volume_matrix_tet_3D(m, m_J[0],
                                                VERTEX(v_i[0]),
                                                VERTEX(v_i[1]),
                                                VERTEX(v_i[2]),
                                                VERTEX(v_i[3]) );
          m_detJ[0] = m;
          for (i = 1; i < 4; i++) { m_detJ[i] = m; m_J[i] = m_J[0]; }
          break;

        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Pyramid<5>::key:
          {
            bool err=false;
            bool useApprox = true;
            if (useApprox)
              {
                double face_centroid[3] = {0,0,0};
                for (unsigned jc=0; jc < 3; ++jc)
                  {
                    for (i = 0; i < 4; ++i)
                      {
                        face_centroid[jc] += 0.25*VERTEX(v_i[i])[jc];
                      }
                  }
                m_J[4].set(0.0);
                DenseMatrix<3,3> J4;
                double detJ4 = 0.0;
                m_detJ[4] = 0;
                for (i = 0; i < 4; ++i)
                  {
                    metric_valid = volume_matrix_tet_3D(m_detJ[i], m_J[i],
                                                        VERTEX(v_i[i]),
                                                        VERTEX(v_i[(i+1)%4]),
                                                        face_centroid,
                                                        VERTEX(v_i[4]) );
                    double fac = 4.0;
                    m_J[i] = m_J[i]*fac;
                    m_detJ[i] *= fac;

                    volume_matrix_tet_3D(detJ4, J4,
                                         VERTEX(v_i[4]),
                                         face_centroid,
                                         VERTEX(v_i[(i+1)%4]),
                                         VERTEX(v_i[i])
                                         );
                    m_J[4] += 0.25*J4*fac;
                    m_detJ[4] += 0.25*detJ4*fac;
                  }
                m = average_metrics(m_detJ, 4);
              }
            else
              {
                for (i = 0; i < 5; ++i) {
                  metric_valid = volume_matrix_pyramid_3D_new(i,
                                                              m_detJ[i], m_J[i],
                                                              VERTEX(v_i[0]),
                                                              VERTEX(v_i[1]),
                                                              VERTEX(v_i[2]),
                                                              VERTEX(v_i[3]),
                                                              VERTEX(v_i[4]));
                  if (m_detJ[i] < MIN_JAC) err=true;
                }

                // FIXME
                m = average_metrics(m_detJ, 5);
              }
            if (0 && (m < MIN_JAC || err))
              {
                std::cout << "pyramid detJ= " << m << std::endl;
                for (i = 0; i < 5; ++i) {
                  std::cout << " detJ[" << i << "]= " << m_detJ[i] << std::endl;
                }
                for (i = 0; i < 5; ++i) {
                  std::cout << " J[" << i << "]= " << m_J[i] << std::endl;
                }
                eMesh.print_entity(element);
                std::cout << "pyramid negative volume, stacktrace=\n" << eMesh.demangled_stacktrace() << std::endl;
                std::cout << "pyramid parts= " << eMesh.print_entity_parts_string(element, "\n");
              }
          }
          break;

        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Wedge<6>::key:
          {
            bool useApprox = true;

            for (i = 0; i < 6; ++i) {
              metric_valid = volume_matrix_wedge_3D(m_detJ[i], m_J[i],
                                                    VERTEX(v_i[locs_prism[i][0]]),
                                                    VERTEX(v_i[locs_prism[i][1]]),
                                                    VERTEX(v_i[locs_prism[i][2]]),
                                                    VERTEX(v_i[locs_prism[i][3]]));
            }
            m = average_metrics(m_detJ, 6);

            if (useApprox)
              {
                double centroid[3] = {0,0,0};
                for (i = 0; i < 6; ++i)
                  {
                    for (int j=0; j < 3; ++j)
                      {
                        centroid[j] += VERTEX(v_i[i])[j]/6.0;
                      }
                  }
                static const int quad_faces[][4] =  {
                  { 0, 1, 4, 3},
                  { 1, 2, 5, 4},
                  { 0, 3, 5, 2} };
                static const int tri_faces[][3] = {
                  { 0, 2, 1},
                  { 3, 4, 5} };

                double morig=m;
                m = 0;
                double vol=0.0;
                for (i = 0; i < 2; ++i)
                  {
                    DenseMatrix<3,3> AJ;
                    double Jac = 0.0;
                    metric_valid = volume_matrix_tet_3D(Jac, AJ,
                                                        centroid ,
                                                        VERTEX(v_i[tri_faces[i][0]]),
                                                        VERTEX(v_i[tri_faces[i][1]]),
                                                        VERTEX(v_i[tri_faces[i][2]])
                                                        );
                    //m_J[i] = m_J[i]*2;
                    vol += Jac/6.0;
                    if (0 && Jac < 0.0)
                      std::cout << "wedge neg Jac= " << Jac << std::endl;
                  }
                for (i = 0; i < 3; ++i)
                  {
                    double face_centroid[3] = {0,0,0};
                    for (int k = 0; k < 4; ++k)
                      {
                        for (int j=0; j < 3; ++j)
                          {
                            face_centroid[j] += VERTEX(v_i[quad_faces[i][k]])[j] / 4.0;
                          }
                      }

                    DenseMatrix<3,3> AJ;
                    double Jac = 0.0;
                    for (int k = 0; k < 4; k++)
                      {
                        metric_valid = volume_matrix_tet_3D(Jac, AJ,
                                                            VERTEX(v_i[quad_faces[i][(k+1)%4]]),
                                                            VERTEX(v_i[quad_faces[i][k]]),
                                                            face_centroid,
                                                            centroid );
                        //m_J[i] = m_J[i]*2;
                        vol += Jac/6.0;
                        if (0 && Jac < 0.0)
                          std::cout << "wedge neg Jac= " << Jac << std::endl;
                      }
                  }
                m = vol*2.0;
                if (0) std::cout << "m= " << m << " morig= " << morig << std::endl;
              }
          }
          break;

        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            metric_valid = volume_matrix_3D(m_detJ[i], m_J[i],
                                              VERTEX(v_i[locs_hex[i][0]]),
                                              VERTEX(v_i[locs_hex[i][1]]),
                                              VERTEX(v_i[locs_hex[i][2]]),
                                              VERTEX(v_i[locs_hex[i][3]]));
          }
          m = average_metrics(m_detJ, 8);
          if (0 && m < MIN_JAC)
            {
              std::cout << "hex detJ= " << m << std::endl;
              for (i = 0; i < 8; ++i) {
                std::cout << " detJ[" << i << "]= " << m_detJ[i] << std::endl;
              }
              for (i = 0; i < 8; ++i) {
                std::cout << " J[" << i << "]= " << m_J[i] << std::endl;
              }
              eMesh.print_entity(element);
              std::cout << "hex negative volume, stacktrace=\n" << eMesh.demangled_stacktrace() << std::endl;
              std::cout << "hex parts= " << eMesh.print_entity_parts_string(element, "\n");
            }

          break;

        case shards::Line<2>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          x2d[0] = VERTEX(v_i[0]);
          x2d[1] = VERTEX(v_i[1]);
          if (spatialDimension == 3)
            metric_valid = volume_matrix_1D_in_3D(m, J, x2d);
          else
            metric_valid = volume_matrix_1D(m, J, x2d);
          for (i = 0; i < 2; i++) { m_detJ[i] = m; m_J[i] = J; }

          break;

          // unimplemented
        case shards::Node::key:
        case shards::Line<3>::key:
        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:

        case shards::Pentagon<5>::key:
        case shards::Hexagon<6>::key:

        default:
          shards::CellTopology topology(topology_data);
          //double *x=0;
          //std::cout << "topology = " << *x;
          std::cout << "topology = " << topology.getName() << "\n" << eMesh.demangled_stacktrace() << std::endl;

          throw std::runtime_error("unknown/unhandled topology in VolumeUtil: "+std::string(topology.getName()));
          break;

        } // end switch over element type

      return metric_valid;
    }

    double VolumeUtil::getJacobianToVolumeScale(shards::CellTopology& cell_topo)
    {
      double volEqui = 1.0;
      switch(cell_topo.getKey() )
        {

          // Tet cells
        case shards::Tetrahedron<4>::key:
        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<10>::key:
          volEqui = 1./6.;
          break;

          // Hex cells
        case shards::Hexahedron<8>::key:
        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
          volEqui = 1.0;
          break;

          // Pyramid cells
        case shards::Pyramid<5>::key:
        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
          volEqui = 1./6.;
          break;

          // Wedge cells
        case shards::Wedge<6>::key:
        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
          volEqui = 1./2.;
          break;

        case shards::Triangle<3>::key:
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
          volEqui = 1./2.;
          break;

        case shards::Quadrilateral<4>::key:
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:
          volEqui = 1./2.;
          break;

        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:
          volEqui = 1.0;
          break;

        default:
          break;
        }//cell key
      return volEqui;
    }

  }


#endif
