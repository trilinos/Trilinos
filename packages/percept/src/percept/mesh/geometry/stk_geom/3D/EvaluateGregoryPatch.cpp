// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GregoryPatch.hpp"
#include "EvaluateGregoryPatch.hpp"
#include "FitGregoryPatches.hpp"

namespace percept {

  static bool s_enableAllDebug = true;
  static double s_maxstep = 10;

  void EvaluateGregoryPatch::project_to_face_boundary(stk::mesh::Entity face, double *xyz_in, double *uv) const
  {
    double xyz[3] = {0,0,0}, xyz_final[3] = {0,0,0};
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    double dmin = std::numeric_limits<double>::max();
    unsigned nn = face_nodes.size();
    VERIFY_OP_ON(nn, != , 2, "bad nn");
    static double uvp_quad[][2] = {{0,0},{1,0},{1,1},{0,1}};
    static double uvp_tri[][2] = {{0,0},{1,0},{0,1}};

    for (unsigned ii=0; ii < nn; ++ii)
      {
        unsigned i0 = ii;
        unsigned i1 = (ii+1) % nn;
        stk::mesh::Entity node0 = face_nodes[i0].entity();
        stk::mesh::Entity node1 = face_nodes[i1].entity();

        double *n0 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , node0 ));
        double *n1 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , node1 ));
        double d0 = Math::distance_3d(xyz_in, n0);
        double d1 = Math::distance_3d(xyz_in, n1);
        double ulocal = 0.0;
        Math::copy_3d(xyz, xyz_in);
        double d01 = Math::project_to_line_3d(n0, n1, xyz, ulocal);
        if (d0 < dmin)
          {
            dmin = d0;
            Math::copy_3d(xyz_final, n0);
            for (unsigned j=0; j < 2; ++j)
              uv[j] = (nn == 4 ? uvp_quad[i0][j] : uvp_tri[i0][j]);
          }
        if (d1 < dmin)
          {
            dmin = d1;
            Math::copy_3d(xyz_final, n1);
            for (unsigned j=0; j < 2; ++j)
              uv[j] = (nn == 4 ? uvp_quad[i1][j] : uvp_tri[i1][j]);
          }
        if (d01 < dmin)
          {
            dmin = d01;
            Math::copy_3d(xyz_final, xyz);
            for (unsigned j=0; j < 2; ++j)
              uv[j] = (nn == 4 ?
                       uvp_quad[i0][j] + ulocal*(uvp_quad[i1][j] - uvp_quad[i0][j])
                       : uvp_tri[i0][j] + ulocal*(uvp_tri[i1][j] - uvp_tri[i0][j])
                       );
          }
      }
    Math::copy_3d(xyz_in, xyz_final);
  }

  //static
  bool EvaluateGregoryPatch::
  evaluateGregoryPatch(const stk::mesh::BulkData& bulk, const double *uv, stk::mesh::Entity face, double *xyz,
           const GregoryControlPointsType *sideset_field, const GregoryControlPointsType *shell_field)
  {
    const GregoryControlPointsType *field = 0;
    if (bulk.entity_rank(face) == stk::topology::ELEMENT_RANK)
      {
        if (!shell_field)
          shell_field = bulk.mesh_meta_data().get_field<GregoryControlPointsType::value_type>(stk::topology::ELEMENT_RANK, "gregory_control_points_shell");
        field = shell_field;
      }
    else
      {
        if (!sideset_field)
          sideset_field = bulk.mesh_meta_data().get_field<GregoryControlPointsType::value_type>(bulk.mesh_meta_data().side_rank(), "gregory_control_points");
        field = sideset_field;
      }
    VERIFY_OP_ON(field, !=, 0, "Null GregoryControlPointsType field");
    double *Cp = stk::mesh::field_data( *field , face);
    if (Cp == 0)
      {
        //PerceptMesh pp(&bulk.mesh_meta_data(), const_cast<stk::mesh::BulkData*>(&bulk), true);
        //std::cout << "Cp = 0 face= " << bulk.identifier(face) << pp.print_part_vector_string(bulk.bucket(face).supersets()) << std::endl;
        return true;
      }
    VERIFY_OP_ON(Cp, !=, 0, "Null GregoryControlPointsType field_data");

    int dim = FitGregoryPatches::MaxControlPoints();
    MDArray Cpx(Cp, 1, dim),
      Cpy(Cp+dim, 1, dim),
      Cpz(Cp+2*dim, 1, dim);

    // std::cout << "Cpx= \n" << Cpx
    //           << "Cpy= \n" << Cpy
    //           << "Cpz= \n" << Cpz  << std::endl;

    switch(bulk.bucket(face).topology().value())
      {
      case stk::topology::SHELL_QUAD_4:
      case stk::topology::QUAD_4:
        xyz[0] = GregoryPatch::evalQuad(uv[0], uv[1], Cpx);
        xyz[1] = GregoryPatch::evalQuad(uv[0], uv[1], Cpy);
        xyz[2] = GregoryPatch::evalQuad(uv[0], uv[1], Cpz);
        break;
      case stk::topology::SHELL_TRI_3:
      case stk::topology::TRI_3:
        xyz[0] = GregoryPatch::evalTri(uv[0], uv[1], Cpx);
        xyz[1] = GregoryPatch::evalTri(uv[0], uv[1], Cpy);
        xyz[2] = GregoryPatch::evalTri(uv[0], uv[1], Cpz);
        break;
      default:
        break;
      }
    return false;
  }

  //static
  bool EvaluateGregoryPatch::
  normalGregoryPatch(const stk::mesh::BulkData& bulk, const double *uv, stk::mesh::Entity face, double *normal,
         const GregoryControlPointsType *sideset_field, const GregoryControlPointsType *shell_field)
  {
    const GregoryControlPointsType *field = 0;
    if (bulk.entity_rank(face) == stk::topology::ELEMENT_RANK)
      {
        if (!shell_field)
          shell_field = bulk.mesh_meta_data().get_field<GregoryControlPointsType::value_type>(stk::topology::ELEMENT_RANK, "gregory_control_points_shell");
        field = shell_field;
      }
    else
      {
        if (!sideset_field)
          sideset_field = bulk.mesh_meta_data().get_field<GregoryControlPointsType::value_type>(bulk.mesh_meta_data().side_rank(), "gregory_control_points");
        field = sideset_field;
      }
    VERIFY_OP_ON(field, !=, 0, "Null GregoryControlPointsType field");
    double *Cp = stk::mesh::field_data( *field , face);
    VERIFY_OP_ON(Cp, !=, 0, "Null GregoryControlPointsType field_data");

    int dim = FitGregoryPatches::MaxControlPoints();
    MDArray Cpx(Cp, 1, dim),
      Cpy(Cp+dim, 1, dim),
      Cpz(Cp+2*dim, 1, dim);

    // std::cout << "Cpx= \n" << Cpx
    //           << "Cpy= \n" << Cpy
    //           << "Cpz= \n" << Cpz  << std::endl;

    double Grad[3][2] = {{0}};
    switch(bulk.bucket(face).topology().value())
      {
      case stk::topology::SHELL_QUAD_4:
      case stk::topology::QUAD_4:
        GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpx, Grad[0]);
        GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpy, Grad[1]);
        GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpz, Grad[2]);
        break;
      case stk::topology::SHELL_TRI_3:
      case stk::topology::TRI_3:
        GregoryPatch::evalTriGrad(uv[0], uv[1], Cpx, Grad[0]);
        GregoryPatch::evalTriGrad(uv[0], uv[1], Cpy, Grad[1]);
        GregoryPatch::evalTriGrad(uv[0], uv[1], Cpz, Grad[2]);
        break;
      default:
        break;
      }
    double GradT[2][3] = {{0}};
    for (unsigned ii=0; ii < 2; ii++)
      for (unsigned jj=0; jj < 3; jj++)
        GradT[ii][jj] = Grad[jj][ii];

    Math::cross_3d(GradT[0], GradT[1], normal);
    Math::normalize_3d(normal);

    return false;
  }

  bool EvaluateGregoryPatch::
  evaluateGregoryPatch(const double *uv, stk::mesh::Entity face, double *xyz)
  {
    return evaluateGregoryPatch(*m_eMesh.get_bulk_data(), uv, face, xyz,
                    m_eMesh.m_gregory_control_points_field, m_eMesh.m_gregory_control_points_field_shell);
  }

  bool EvaluateGregoryPatch::
  normalGregoryPatch(const double *uv, stk::mesh::Entity face, double *norm)
  {
    return normalGregoryPatch(*m_eMesh.get_bulk_data(), uv, face, norm,
                    m_eMesh.m_gregory_control_points_field, m_eMesh.m_gregory_control_points_field_shell);
  }

  // finite-difference gradient of given function F
  template<class F>
  class FDGradient {
  public:
    std::string getName() { return "FDGradient"; }

    void operator()(F& f, const double u[F::Size], double g[F::Size]) const
    {
      double eps = f.eps();
      unsigned n = F::Size;
      double v[n];
      for (unsigned i = 0; i < n; ++i)
        {
          v[i] = u[i];
        }
      for (unsigned i = 0; i < n; ++i)
        {
          v[i] = u[i] + eps;
          double fp = f(v);
          v[i] = u[i] - eps;
          double fm = f(v);
          g[i] = (fp -fm)/(2.0*eps);
          v[i] = u[i];
        }
    }
  };

  // finite-difference Hessian of given function F
  template<class F>
  class FDHessian {
  public:
    void operator()(F& f, const double u[F::Size], double h[F::Size][F::Size]) const
    {
      double eps = f.eps();
      unsigned n = F::Size;
      double v[n];
      for (unsigned i = 0; i < n; ++i)
        {
          for (unsigned j = 0; j < n; ++j)
            {
              v[i] = u[i];
              v[j] = u[j];
              v[i] += eps;
              v[j] += eps;
              double fpp = f(v);

              v[i] = u[i];
              v[j] = u[j];
              v[i] -= eps;
              v[j] += eps;
              double fmp = f(v);

              v[i] = u[i];
              v[j] = u[j];
              v[i] += eps;
              v[j] -= eps;
              double fpm = f(v);

              v[i] = u[i];
              v[j] = u[j];
              v[i] -= eps;
              v[j] -= eps;
              double fmm = f(v);

              h[i][j] = (fpp - fmp - fpm + fmm)/(4.0*eps*eps);
            }
        }
    }
  };

  /// the cost function (distance-squared from evaluated point to input point)
  struct FGregory {
    PerceptMesh& m_eMesh;
    stk::mesh::Entity m_face;
    double xyz[3];
    const double *m_point;
    FGregory(PerceptMesh& eMesh, stk::mesh::Entity face, const double *point) : m_eMesh(eMesh), m_face(face),m_point(point) {}
    enum { Size = 2 };
    std::string getName() { return "FGregory"; }
    double eps() const { return 1.e-6; }
    double operator()(const double u[Size])
    {
      EvaluateGregoryPatch::evaluateGregoryPatch(*m_eMesh.get_bulk_data(), u, m_face, xyz,
                                     m_eMesh.m_gregory_control_points_field, m_eMesh.m_gregory_control_points_field_shell);

      double sum=0.0;
      for (unsigned ii=0; ii < 3; ++ii)
        {
          sum += (xyz[ii] - m_point[ii])*(xyz[ii] - m_point[ii]);
        }
      return sum;
    }
  };

  /// the cost function (distance-squared from evaluated point to input point)
  struct FLinear {
    PerceptMesh& m_eMesh;
    stk::mesh::Entity m_face;
    double xyz[3];
    const double *m_point;
    FLinear(PerceptMesh& eMesh, stk::mesh::Entity face, const double *point) : m_eMesh(eMesh), m_face(face),m_point(point) {}
    enum { Size = 2 };
    std::string getName() { return "FLinear"; }
    double eps() const { return 1.e-6; }
    double operator()(const double uu[Size])
    {
      double u=uu[0];
      double v=uu[1];

      const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), m_face, stk::topology::NODE_RANK );
      unsigned nn = face_nodes.size();
      VERIFY_OP_ON(nn, != , 2, "bad nn");

      double basesTri[] = { 1.0-u-v, u, v};
      double basesQuad[] = { (1.0-u)*(1.0-v), u*(1.0-v), u*v, (1.0-u)*v};
      double *bases = (nn == 3 ? basesTri : basesQuad);
      for (unsigned jc = 0; jc < 3; jc++)
        {
          xyz[jc] = 0.0;
        }
      for (unsigned ib = 0; ib < nn; ++ib)
        {
          double *nc = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , face_nodes[ib].entity() ));
          for (unsigned jc = 0; jc < 3; jc++)
            {
              xyz[jc] += bases[ib]*nc[jc];
            }
        }
      double sum=0.0;
      for (unsigned ii=0; ii < 3; ++ii)
        {
          sum += (xyz[ii] - m_point[ii])*(xyz[ii] - m_point[ii]);
        }
      return sum;
    }
  };

  // interfaces to analytic gradient
  void EvaluateGregoryPatch::evaluateLinearQuad(PerceptMesh& eMesh, const double *uv, stk::mesh::Entity face, double *xyz, double G[3][2], double H[3][2][2])
  {
    double u=uv[0], v=uv[1];
    double bases[4] = {(1-u)*(1-v), u*(1-v), u*v, (1-u)*v};
    double dbases[4][2] = {{-(1-v),-(1-u)}, {(1-v),-u}, {v, u}, {-v, (1-u)}};
    double ddbases[4][2][2] = {{{0,1},{1,0}}, {{0,-1},{-1,0}}, {{0,1},{1,0}}, {{0,-1},{-1,0}}};
    const MyPairIterRelation nodes(eMesh, face, stk::topology::NODE_RANK);

    for (unsigned jj=0; jj < 3; ++jj)
      {
        xyz[jj] = 0.0;
        for (unsigned kk=0; kk < 2; ++kk)
          {
            G[jj][kk] = 0.0;
            for (unsigned mm=0; mm < 2; ++mm)
              {
                H[jj][kk][mm] = 0.0;
              }
          }
      }
    for (unsigned ii=0; ii < nodes.size(); ++ii)
      {
        double *nc = static_cast<double*>(stk::mesh::field_data(*eMesh.get_coordinates_field(), nodes[ii].entity()));
        for (unsigned jj=0; jj < 3; ++jj)
          {
            xyz[jj] += bases[ii]*nc[jj];
            for (unsigned kk=0; kk < 2; ++kk)
              {
                G[jj][kk] += dbases[ii][kk]*nc[jj];
                for (unsigned mm=0; mm < 2; ++mm)
                  {
                    H[jj][kk][mm] += ddbases[ii][kk][mm]*nc[jj];
                  }
              }
          }
      }
  }

  // interfaces to analytic gradient
  void evaluateLinearTri(PerceptMesh& eMesh, const double *uv, stk::mesh::Entity face, double *xyz, double G[3][2], double H[3][2][2])
  {
    double u=uv[0], v=uv[1];
    double bases[3] = {1 - u - v, u, v};
    double dbases[3][2] = {{-1,-1},{1,0},{0,1}};
    const MyPairIterRelation nodes(eMesh, face, stk::topology::NODE_RANK);

    for (unsigned jj=0; jj < 3; ++jj)
      {
        xyz[jj] = 0.0;
        for (unsigned kk=0; kk < 2; ++kk)
          {
            G[jj][kk] = 0.0;
            for (unsigned mm=0; mm < 2; ++mm)
              {
                H[jj][kk][mm] = 0.0;
              }
          }
      }
    for (unsigned ii=0; ii < nodes.size(); ++ii)
      {
        double *nc = static_cast<double*>(stk::mesh::field_data(*eMesh.get_coordinates_field(), nodes[ii].entity()));
        for (unsigned jj=0; jj < 3; ++jj)
          {
            xyz[jj] += bases[ii]*nc[jj];
            for (unsigned kk=0; kk < 2; ++kk)
              {
                G[jj][kk] += dbases[ii][kk]*nc[jj];
              }
          }
      }
  }

  // interface to analytic gradient
  template<class F>
  class LinearGradient {
  public:
    std::string getName() { return "LinearGradient"; }

    void operator()(F& f, const double u[F::Size], double g[F::Size]) const
    {
      stk::mesh::BulkData& bulk = *f.m_eMesh.get_bulk_data();
      const double *uv = &u[0];
      stk::mesh::Entity face = f.m_face;
      double *xyz = &f.xyz[0];
      const double *point = f.m_point;
      double G[3][2], H[3][2][2];

      switch(bulk.bucket(face).topology().value())
        {
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::QUAD_4:

          EvaluateGregoryPatch::evaluateLinearQuad(f.m_eMesh, uv, face, xyz, G, H);
          break;

        case stk::topology::SHELL_TRI_3:
        case stk::topology::TRI_3:

          evaluateLinearTri(f.m_eMesh, uv, face, xyz, G, H);
          break;

        default:
          throw std::runtime_error("unrecognized topology");
          break;
        }

      // (x_i - y_i)(x_i - y_i)
      // gj = 2 (x_i - y_i)* dx_i_du_j
      // hjk = dgj_duk =  2( dx_i_duk * dx_i_du_j )+ 2 (x_i - y_i) * dx_i_du_jk

      for (unsigned jj = 0; jj < 2; ++jj)
        {
          double sum=0.0;
          for (unsigned ii=0; ii < 3; ++ii)
            {
              sum += 2.0*(xyz[ii] - point[ii])*G[ii][jj];
            }
          g[jj] = sum;
        }
    }
  };

  // interface to analytic Hessian
  template<class F>
  class LinearHessian {
  public:
    void operator()(F& f, const double u[F::Size], double h[F::Size][F::Size]) const
    {
      stk::mesh::BulkData& bulk = *f.m_eMesh.get_bulk_data();
      const double *uv = &u[0];
      stk::mesh::Entity face = f.m_face;
      double *xyz = &f.xyz[0];
      const double *point = f.m_point;

      double G[3][2];
      double H[3][2][2];
      switch(bulk.bucket(face).topology().value())
        {
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::QUAD_4:

          EvaluateGregoryPatch::evaluateLinearQuad(f.m_eMesh, uv, face, xyz, G, H);
          break;

        case stk::topology::SHELL_TRI_3:
        case stk::topology::TRI_3:

          evaluateLinearTri(f.m_eMesh, uv, face, xyz, G, H);
          break;

        default:
          throw std::runtime_error("unrecognized topology");
          break;
        }

      // (x_i - y_i)(x_i - y_i)
      // gj = 2 (x_i - y_i)* dx_i_du_j
      // hjk = dgj_duk =  2( dx_i_duk * dx_i_du_j )+ 2 (x_i - y_i) * dx_i_du_jk

      for (unsigned kk = 0; kk < 2; ++kk)
        {
          for (unsigned jj = 0; jj < 2; ++jj)
            {
              double sum=0.0;
              for (unsigned ii=0; ii < 3; ++ii)
                {
                  sum += 2.0*(xyz[ii] - point[ii])*H[ii][jj][kk] + 2.0*(G[ii][kk] * G[ii][jj]);
                }
              h[jj][kk] = sum;
            }
        }
    }
  };

  template<class F>
  class GPGradient {
  public:
    std::string getName() { return "GPGradient"; }

    void operator()(F& f, const double u[F::Size], double g[F::Size]) const
    {
      stk::mesh::BulkData& bulk = *f.m_eMesh.get_bulk_data();
      const double *uv = &u[0];
      stk::mesh::Entity face = f.m_face;
      double *xyz = &f.xyz[0];
      const double *point = f.m_point;
      GregoryControlPointsType *field = 0;
      if (bulk.entity_rank(face) == stk::topology::ELEMENT_RANK)
        {
          field = f.m_eMesh.m_gregory_control_points_field_shell;
        }
      else
        {
          field = f.m_eMesh.m_gregory_control_points_field;
        }
      VERIFY_OP_ON(field, !=, 0, "Null GregoryControlPointsType field");
      double *Cp = stk::mesh::field_data( *field , face);
      if (!Cp)
        {
          std::cout << "Null GregoryControlPointsType field_data for face = " << bulk.identifier(face)
                    << " parts= " << f.m_eMesh.print_entity_parts_string(face) << std::endl;
        }
      VERIFY_OP_ON(Cp, !=, 0, "Null GregoryControlPointsType field_data");

    int dim = FitGregoryPatches::MaxControlPoints();
    MDArray Cpx(Cp, 1, dim),
      Cpy(Cp+dim, 1, dim),
      Cpz(Cp+2*dim, 1, dim);

      double G[3][2];

      switch(bulk.bucket(face).topology().value())
        {
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::QUAD_4:

          xyz[0] = GregoryPatch::evalQuad(uv[0], uv[1], Cpx);
          xyz[1] = GregoryPatch::evalQuad(uv[0], uv[1], Cpy);
          xyz[2] = GregoryPatch::evalQuad(uv[0], uv[1], Cpz);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpx, G[0]);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpy, G[1]);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpz, G[2]);
          break;

        case stk::topology::SHELL_TRI_3:
        case stk::topology::TRI_3:

          xyz[0] = GregoryPatch::evalTri(uv[0], uv[1], Cpx);
          xyz[1] = GregoryPatch::evalTri(uv[0], uv[1], Cpy);
          xyz[2] = GregoryPatch::evalTri(uv[0], uv[1], Cpz);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpx, G[0]);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpy, G[1]);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpz, G[2]);
          break;

        default:
          throw std::runtime_error("unrecognized topology");
          break;
        }

      // (x_i - y_i)(x_i - y_i)
      // gj = 2 (x_i - y_i)* dx_i_du_j
      // hjk = dgj_duk =  2( dx_i_duk * dx_i_du_j )+ 2 (x_i - y_i) * dx_i_du_jk

      for (unsigned jj = 0; jj < 2; ++jj)
        {
          double sum=0.0;
          for (unsigned ii=0; ii < 3; ++ii)
            {
              sum += 2.0*(xyz[ii] - point[ii])*G[ii][jj];
            }
          g[jj] = sum;
        }
    }
  };

  // interface to analytic Hessian
  template<class F>
  class GPHessian {
  public:
    void operator()(F& f, const double u[F::Size], double h[F::Size][F::Size]) const
    {
      stk::mesh::BulkData& bulk = *f.m_eMesh.get_bulk_data();
      const double *uv = &u[0];
      stk::mesh::Entity face = f.m_face;
      double *xyz = &f.xyz[0];
      const double *point = f.m_point;
      GregoryControlPointsType *field = 0;
      if (bulk.entity_rank(face) == stk::topology::ELEMENT_RANK)
        {
          field = f.m_eMesh.m_gregory_control_points_field_shell;
        }
      else
        {
          field = f.m_eMesh.m_gregory_control_points_field;
        }
      VERIFY_OP_ON(field, !=, 0, "Null GregoryControlPointsType field");
      
      double *Cp = stk::mesh::field_data( *field , face);
      VERIFY_OP_ON(Cp, !=, 0, "Null GregoryControlPointsType field_data");

      int dim = FitGregoryPatches::MaxControlPoints();
      MDArray Cpx(Cp, 1, dim),
        Cpy(Cp+dim, 1, dim),
        Cpz(Cp+2*dim, 1, dim);

      double G[3][2];
      double H[3][2][2];
      switch(bulk.bucket(face).topology().value())
        {
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::QUAD_4:

          xyz[0] = GregoryPatch::evalQuad(uv[0], uv[1], Cpx);
          xyz[1] = GregoryPatch::evalQuad(uv[0], uv[1], Cpy);
          xyz[2] = GregoryPatch::evalQuad(uv[0], uv[1], Cpz);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpx, G[0]);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpy, G[1]);
          GregoryPatch::evalQuadGrad(uv[0], uv[1], Cpz, G[2]);
          GregoryPatch::evalQuadHessian(uv[0], uv[1], Cpx, H[0]);
          GregoryPatch::evalQuadHessian(uv[0], uv[1], Cpy, H[1]);
          GregoryPatch::evalQuadHessian(uv[0], uv[1], Cpz, H[2]);
          break;

        case stk::topology::SHELL_TRI_3:
        case stk::topology::TRI_3:

          xyz[0] = GregoryPatch::evalTri(uv[0], uv[1], Cpx);
          xyz[1] = GregoryPatch::evalTri(uv[0], uv[1], Cpy);
          xyz[2] = GregoryPatch::evalTri(uv[0], uv[1], Cpz);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpx, G[0]);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpy, G[1]);
          GregoryPatch::evalTriGrad(uv[0], uv[1], Cpz, G[2]);
          GregoryPatch::evalTriHessian(uv[0], uv[1], Cpx, H[0]);
          GregoryPatch::evalTriHessian(uv[0], uv[1], Cpy, H[1]);
          GregoryPatch::evalTriHessian(uv[0], uv[1], Cpz, H[2]);
          break;

        default:
          throw std::runtime_error("unrecognized topology");
          break;
        }

      // (x_i - y_i)(x_i - y_i)
      // gj = 2 (x_i - y_i)* dx_i_du_j
      // hjk = dgj_duk =  2( dx_i_duk * dx_i_du_j )+ 2 (x_i - y_i) * dx_i_du_jk

      for (unsigned kk = 0; kk < 2; ++kk)
        {
          for (unsigned jj = 0; jj < 2; ++jj)
            {
              double sum=0.0;
              for (unsigned ii=0; ii < 3; ++ii)
                {
                  sum += 2.0*(xyz[ii] - point[ii])*H[ii][jj][kk] + 2.0*(G[ii][kk] * G[ii][jj]);
                }
              h[jj][kk] = sum;
            }
        }
    }
  };

  // found_uv is also the initial guess
  template<class F, class Gradient, class Hessian>
  bool EvaluateGregoryPatch::
  findClosestPointInternal(const double *input_xyz, stk::mesh::Entity face, double *closest_xyz, double *found_uv)
  {
    double uv[2] = {found_uv[0], found_uv[1]};
    if (0 && !m_eMesh.in_face(face, uv, 0.0))
      {
        Math::copy_3d(closest_xyz, input_xyz);
        project_to_face_boundary(face, closest_xyz, uv);
      }
    double G[2], H[2][2];
    double tolOrig = 1.e-5;

    m_errorMsg = "";

    F func(m_eMesh, face, input_xyz);
    Gradient gradient;
    Hessian hessian;

    if (m_debug && s_enableAllDebug) std::cout << "findClosestPointInternal face= " << m_eMesh.id(face) << " topo= " << m_eMesh.topology(face) << std::endl;

    double el_ave = m_eMesh.edge_length_ave(face, m_eMesh.get_coordinates_field());
    for (unsigned iter=0; iter < 100; ++iter)
      {

        double tol = tolOrig;
        bool in_face = m_eMesh.in_face(face, uv, 1.e-4);
        if (!in_face) tol = 1.e-3;

        gradient(func, uv, G);
        double gnorm = std::sqrt(G[0]*G[0] + G[1]*G[1]);
        hessian(func, uv, H);
        double detH = (H[0][0]*H[1][1] - H[1][0]*H[0][1]);
        // if (std::fabs(detH) < 1.e-12*el_ave)
        //   {
        //     m_errorMsg = "det Hessian too small = " + toString(detH);
        //   }
        double du = -( H[1][1]*G[0] - H[0][1]*G[1])/detH;
        double dv = -(-H[1][0]*G[0] + H[0][0]*G[1])/detH;
        double dunorm = std::sqrt(du*du + dv*dv);
        if (dunorm > 10.)
          {
            m_errorMsg = " dunorm too big = " + toString(dunorm);
          }
        double maxStep = s_maxstep;
        if (dunorm > maxStep)
          {
            du = maxStep*du/dunorm;
            dv = maxStep*dv/dunorm;
          }
        dunorm = std::sqrt(du*du + dv*dv);

        if (m_debug && s_enableAllDebug) std::cout << gradient.getName() << ": iter= " << iter << " err= " << m_errorMsg << " uv= " << uv[0] << ", " << uv[1]
                               << " gnorm= " << gnorm << " dunorm= " << dunorm << " du= " << du << " dv= " << dv << " detH= " << detH
                               << " G= " << G[0] << ", " << G[1]
                               << " H= " << H[0][0] << ", " << H[0][1] << ", " << H[1][0] << ", " << H[1][1]
                                                   << " el_ave= " << el_ave << " tol*el_ave= " << tol*el_ave
                               << std::endl;
        // if (m_errorMsg.length())
        //   {
        //     if (m_debug && s_enableAllDebug) std::cout << m_errorMsg << std::endl;
        //     return true;
        //   }
        if (gnorm < tol*el_ave && dunorm < tol)
          {
            double dist = func(uv);
            for (unsigned j=0; j < 3; ++j)
              closest_xyz[j] = func.xyz[j];
            if (found_uv)
              {
                found_uv[0] = uv[0];
                found_uv[1] = uv[1];
              }
            if (m_debug && s_enableAllDebug) std::cout << "CONVERGED: iter= " << iter << " gnorm= " << gnorm << " uv= " << uv[0] << ", " << uv[1]
                                                       << " du= " << du << " dv= " << dv << " dist= " << dist << std::endl;
            return false;
          }

        uv[0] += du;
        uv[1] += dv;

      }
    m_errorMsg += ":too many iterations:";
    return true;
  }


  bool EvaluateGregoryPatch::
  findClosestPoint(const double *input_xyz, stk::mesh::Entity face, double *closest_xyz, double *found_uv, bool linearOnly,
                   stk::mesh::Entity node )
  {
    double centroid[2] = {0.5,0.5};
    stk::topology face_topo = m_eMesh.bucket(face).topology();
    VERIFY_OP_ON(percept::FitGregoryPatches::is_surface_topology(face_topo), ==, true, "bad face_topo");

    if (face_topo.value() == stk::topology::TRI_3 || face_topo.value() == stk::topology::SHELL_TRI_3)
      {
        centroid[0] = 1./3.;
        centroid[1] = 1./3.;
      }
    double tol_in_face = 1.e-4;
    std::string errFound = "";
    double uv[2] = {centroid[0], centroid[1]};
    if (m_debug && s_enableAllDebug)
      {
        std::cout << "findClosestPoint: input_xyz= " << Math::print_3d(input_xyz, 14) << std::endl;
      }

    // use linear face to get good initial guess
    bool err = true;
    bool errLin = findClosestPointInternal<FLinear, LinearGradient<FLinear>, LinearHessian<FLinear> >(input_xyz, face, closest_xyz, uv);
    if (errLin) errFound += " ||case 0: Linear not converged";
    bool in_face = m_eMesh.in_face(face, uv, tol_in_face);
    if (!in_face)
      {
        errFound += " ||case 0: Linear not in face uv= "+Math::print_2d(uv);
      }
    if (linearOnly)
      {
        if (!in_face)
          {
            double uv_loc[2]={0,0};
            Math::copy_3d(closest_xyz, input_xyz);
            project_to_face_boundary(face, closest_xyz, uv_loc);
            if (found_uv)
              {
                found_uv[0] = uv_loc[0];
                found_uv[1] = uv_loc[1];
              }
          }
        else
          {
            if (found_uv)
              {
                found_uv[0] = uv[0];
                found_uv[1] = uv[1];
              }
          }
        m_errorMsg += errFound;

        return errLin;
      }

    if (!in_face)
      {
        Math::copy_3d(closest_xyz, input_xyz);
        project_to_face_boundary(face, closest_xyz, uv);
      }

    err = findClosestPointInternal<FGregory, GPGradient<FGregory>, GPHessian<FGregory> >(input_xyz, face, closest_xyz, uv);

    if (found_uv)
      {
        found_uv[0] = uv[0];
        found_uv[1] = uv[1];
      }
    m_errorMsg += errFound;

    return err;
  }
}
