// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SmootherMetric_hpp
#define SmootherMetric_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/PerceptMesh.hpp>
#include <percept/MeshType.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/mesh/mod/smoother/sgrid_metric_helpers.hpp>

  namespace percept {

    enum CombineOp {
      COP_SUM,
      COP_MIN,
      COP_MAX
    };

    template<typename MeshType>
    class SmootherMetricBase
    {

    public:

      virtual ~SmootherMetricBase() = default;
      SmootherMetricBase(PerceptMesh *eMesh) : m_topology_data(0), m_eMesh(eMesh), m_node(), m_is_nodal(false), m_combine(COP_SUM), m_debug(false) {}

      virtual double length_scaling_power() { return 1.0; }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(typename MeshType::MTElement element, bool& valid)=0;

      virtual double grad_metric(typename MeshType::MTElement /*element*/, bool& /*valid*/, double /*grad*/[8][4]) {
        throw std::runtime_error("not impl SmootherMetric::grad_metric");
        return 0.0;
      }
      virtual bool has_gradient() { return false; }
      void set_combine_op(CombineOp combine) { m_combine= combine; }
      CombineOp get_combine_op() { return m_combine; }

      // ------ VVVVVVV --------
      virtual double grad_and_hessian(typename MeshType::MTElement /*element*/, bool& /*valid*/, double /*grad*/[8][4], double /*hess*/[8][4][8][4], const bool /*getGrad*/=false, const bool /*getHess*/=false)
      {
        throw std::runtime_error("not impl SmootherMetric::grad_and_hessian");
        return 0.0;
      }

      virtual bool has_gradient_and_hessian() { return false; }
      const typename MeshType::MTCellTopology * m_topology_data ;
      virtual void set_node(typename MeshType::MTNode node) { m_node=node; }
      typename MeshType::MTNode get_node() { return m_node; }
      bool is_nodal() { return m_is_nodal; }

    public:
      PerceptMesh *m_eMesh;
      typename MeshType::MTNode m_node; // for metrics that are node-based
      bool m_is_nodal;
      CombineOp m_combine;
      typename MeshType::MTField *m_coord_field_current;
      typename MeshType::MTField *m_coord_field_original;

      bool m_debug;
    };


    template<typename MeshType>
    class SmootherMetricImpl : public SmootherMetricBase<MeshType>
    {
    public:
      SmootherMetricImpl(PerceptMesh *eMesh);
    };

    template<>
    class SmootherMetricImpl<STKMesh> : public SmootherMetricBase<STKMesh>
    {
    public:
      SmootherMetricImpl(PerceptMesh *eMesh) : SmootherMetricBase<STKMesh>(eMesh)
      {
        m_coord_field_current   = eMesh->get_coordinates_field();
        m_coord_field_original  = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_NM1");
        VERIFY_OP_ON(m_coord_field_original, !=, 0, "m_coord_field_original");
      }
    };

    using SmootherMetric = SmootherMetricImpl<STKMesh>;

    class SmootherMetricFunction : public SmootherMetric {
    public:
      SmootherMetric *m_smoother;
      double m_nodal_edge_length_ave;

      SmootherMetricFunction(PerceptMesh* eMesh, SmootherMetric *smoother) : SmootherMetric(eMesh), m_smoother(smoother) , m_nodal_edge_length_ave(0)
      {
        m_is_nodal = true;
      }

      double nodal_edge_length_ave()
      {
        //int spatialDim = m_eMesh->get_spatial_dim();
        double nm=0.0;

        VERIFY_OP_ON(m_eMesh, !=, 0, "eMesh");
        const MyPairIterRelation node_elems(*m_eMesh, m_node, m_eMesh->element_rank() );
        double nele = 0.0;
        for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
          {
            stk::mesh::Entity element = node_elems[i_elem].entity();
            if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
              continue;
            double elem_edge_len = m_eMesh->edge_length_ave(element, m_coord_field_original);
            nm += elem_edge_len;
            nele += 1.0;
          }
        nm /= nele;
        return nm;
      }

      // size, eps, etc. (Size is 3 even for 2D)
      enum { Size = 3 };
      std::string getName() { return "SmootherMetricFunction"; }
      double eps() const
      {
        VERIFY_OP_ON(m_nodal_edge_length_ave, >, 0.0, "bad m_nodal_edge_length_ave");
        return 1.e-6*m_nodal_edge_length_ave;
      }
      int get_n() { return m_eMesh->get_spatial_dim(); }

      virtual void set_node(stk::mesh::Entity node)
      {
        m_node = node;
        VERIFY_OP_ON(m_eMesh->is_valid(m_node), ==, true, "node valid");
        m_nodal_edge_length_ave = nodal_edge_length_ave();
      }

      double operator()(const double u[Size])
      {
        VERIFY_OP_ON(Size, >=, m_eMesh->get_spatial_dim(), " bad dim");
        double usave[Size];
        double *coords = m_eMesh->field_data(m_eMesh->get_coordinates_field(), m_node);
        for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
          {
            usave[jc] = coords[jc];
            coords[jc] = u[jc];
          }

        double mm = 0.0;
        const MyPairIterRelation node_elems(*m_eMesh, m_node, m_eMesh->element_rank());
        for (unsigned ielem=0; ielem < node_elems.size(); ++ielem)
          {
            stk::mesh::Entity element = node_elems[ielem].entity();
            bool lvalid = true;
            m_smoother->m_topology_data = m_eMesh->get_cell_topology(m_eMesh->bucket(element));

            double lm = m_smoother->metric(element, lvalid);
            mm += lm;
          }
        for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
          {
            coords[jc] = usave[jc];
          }
        return mm;
      }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid)
      {
#if 0
        double *coords = m_eMesh->field_data(m_eMesh->get_coordinates_field(), m_node);
        double u[Size];
        for (int jc=0; jc < m_eMesh.get_spatial_dim(); ++jc)
          {
            u[jc] = coords[jc];
          }
        return this->operator()(u);
#endif
        return m_smoother->metric(element, valid);
      }

      //virtual double grad_metric(stk::mesh::Entity element, bool& valid, double grad[8][4]) { throw std::runtime_error("not impl"); return 0.0; }

      //virtual bool has_gradient() { return false; }


    };

    class SmootherMetricElementFunction : public SmootherMetric {
    public:

      SmootherMetric *m_smoother;
      double m_element_edge_length_ave;
      stk::mesh::Entity m_element;
      unsigned m_nnodes;
      unsigned m_ndim;
      unsigned m_ndof;
      stk::mesh::FieldBase *m_cg_lambda_field;
      stk::mesh::FieldBase *m_cg_normal_field;
      double macheps;
      double sqrt_eps;
      double cubert_eps;

      SmootherMetricElementFunction(PerceptMesh* eMesh, SmootherMetric *smoother) : SmootherMetric(eMesh), m_smoother(smoother) , m_element_edge_length_ave(0)
      {
        m_is_nodal = false;
        m_ndim = unsigned(eMesh->get_spatial_dim());
        m_ndof = m_ndim;
        m_cg_lambda_field    = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_lambda");
        m_cg_normal_field    = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");

        macheps = std::numeric_limits<double>::epsilon();
        sqrt_eps = std::sqrt(macheps);
        cubert_eps = std::pow(macheps, 1./3.);
      }

      double element_edge_length_ave()
      {
        VERIFY_OP_ON(m_eMesh, !=, 0, "eMesh");
        double emin=0.0, emax=0.0;
        double elave = m_eMesh->edge_length_ave(m_element, m_coord_field_original, &emin, &emax);
        (void)elave;
        return elave;
      }

      // size, eps, etc. (Size is 3 even for 2D)
      enum { Size = 3*8 };  // maximum size for hex elements
      std::string getName() { return "SmootherMetricElementFunction"; }
      double eps() const
      {
        VERIFY_OP_ON(m_element_edge_length_ave, >, 0.0, "bad m_element_edge_length_ave");
        //return sqrt_eps*m_element_edge_length_ave;
        return cubert_eps*m_element_edge_length_ave;
      }
      int get_n() { return m_ndof*m_nnodes; }

      virtual void set_element(stk::mesh::Entity element)
      {
        m_element = element;
        m_nnodes = m_eMesh->get_bulk_data()->num_nodes(element);
        VERIFY_OP_ON(m_eMesh->is_valid(m_element), ==, true, "element valid");
        m_element_edge_length_ave = element_edge_length_ave();
      }

      double operator()(const double *u)
      {
        int N = get_n();
        VERIFY_OP_ON(Size, >=, m_eMesh->get_spatial_dim()*m_nnodes, " bad dim");
        double *usave = new double[N];
        stk::mesh::Entity const *nodes = m_eMesh->get_bulk_data()->begin_nodes(m_element);
        double **coordsAll = new double*[m_nnodes];
        for (unsigned inode=0; inode < m_nnodes; ++inode)
          {
            stk::mesh::Entity node = nodes[inode];
            double *coords = m_eMesh->field_data(m_eMesh->get_coordinates_field(), node);
            double Lambda = (m_cg_lambda_field? *static_cast<double*>(stk::mesh::field_data( *m_cg_lambda_field, node )) : 0.0);
            //double *norm = (m_cg_normal_field? *static_cast<double*>(stk::mesh::field_data( *m_cg_normal_field, node )) : 0.0);

            coordsAll[inode] = coords;
            for (unsigned jc=0; jc < m_ndim; ++jc)
              {
                usave[inode*m_ndof + jc] = coords[jc];
                coords[jc] = u[inode*m_ndof + jc];
              }
            usave[inode*m_ndof + m_ndim] = Lambda;
          }

        bool lvalid = true;
        m_smoother->m_debug = m_debug;
        m_smoother->m_topology_data = m_eMesh->get_cell_topology(m_eMesh->bucket(m_element));

        double mm = m_smoother->metric(m_element, lvalid);

        for (unsigned inode=0; inode < m_nnodes; ++inode)
          {
            for (unsigned jc=0; jc < m_ndim; ++jc)
              {
                coordsAll[inode][jc] = usave[inode*m_ndof + jc];
              }
          }
        delete [] coordsAll;
        delete [] usave;
        return mm;
      }

      void analytic_gradient(const double *u, double *grad)
      {
        int N = get_n();
        VERIFY_OP_ON(Size, >=, m_eMesh->get_spatial_dim()*m_nnodes, " bad dim");
        double *usave = new double[N];
        stk::mesh::Entity const *nodes = m_eMesh->get_bulk_data()->begin_nodes(m_element);
        double **coordsAll = new double*[m_nnodes];
        int spatialDim = m_eMesh->get_spatial_dim();
        for (unsigned inode=0; inode < m_nnodes; ++inode)
          {
            stk::mesh::Entity node = nodes[inode];
            double *coords = m_eMesh->field_data(m_eMesh->get_coordinates_field(), node);
            double Lambda = (m_cg_lambda_field? *static_cast<double*>(stk::mesh::field_data( *m_cg_lambda_field , node )) : 0.0);
            coordsAll[inode] = coords;
            for (int jc=0; jc < spatialDim; ++jc)
              {
                usave[inode*m_ndof + jc] = coords[jc];
                coords[jc] = u[inode*m_ndof + jc];
              }
            usave[inode*m_ndof + m_ndim] = Lambda;
          }

        bool lvalid = true;
        m_smoother->m_debug = m_debug;
        m_smoother->m_topology_data = m_eMesh->get_cell_topology(m_eMesh->bucket(m_element));

        double gradAll[8][4];
        m_smoother->grad_metric(m_element, lvalid, gradAll);

        for (unsigned inode=0; inode < m_nnodes; ++inode)
          {
            for (int jc=0; jc < spatialDim; ++jc)
              {
                grad[inode*m_ndof + jc] = gradAll[inode][jc];
                coordsAll[inode][jc] = usave[inode*m_ndof + jc];
              }
            grad[inode*m_ndof + m_ndim] = gradAll[inode][m_ndim];
          }
        delete [] coordsAll;
        delete [] usave;
      }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid)
      {
        return m_smoother->metric(element, valid);
      }
    };

    struct FDTest
    {
      bool m_debug;
      double deriv[2], hess[4];
      FDTest() : m_debug(false) {}

      double eps() { return 1.e-5; }
      unsigned get_n() { return 2; }
      double operator()(double *v)
      {
        double x = v[0], y = v[1];
        deriv[0] = y + 1.0 + 2.0*x;
        deriv[1] = x + 2.0 + 2.0*y;
        hess[0*2 + 0] = 2.0;
        hess[0*2 + 1] = 1.0;
        hess[1*2 + 0] = 1.0;
        hess[1*2 + 1] = 2.0;
        return x*y + x + 2.0*y + x*x + y*y;
      }
    };

    // finite-difference gradient of given function F
    template<class F>
    class FDGradient {
    public:

      void test()
      {
        FDGradient<FDTest> fd;
        double v[2] = {0.1, 0.2};
        double g[2];
        FDTest t;
        fd(t, v, g);
        t(v);
        VERIFY_OP_ON(std::fabs(t.deriv[0] - g[0]) , < , 1.e-6, "bad g");
        VERIFY_OP_ON(std::fabs(t.deriv[1] - g[1]) , < , 1.e-6, "bad g 1");
        std::cout << "OK g= " << g[0] << " " << g[1] << std::endl;
      }

      std::string getName() { return "FDGradient"; }

      void operator()(F& f, const double *u, double *g) const
      {
        double eps = f.eps();
        unsigned n = f.get_n();
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

      void test()
      {
        FDHessian<FDTest> fd;
        double v[2] = {0.1, 0.2};
        double h[4];
        FDTest t;
        fd(t, v, h);
        t(v);
        if (std::fabs(t.hess[0] - h[0])  >  1.e-6) std::cout << "t.hess[0] = " << t.hess[0] << " h[0]= " << h[0] << std::endl;

        VERIFY_OP_ON(std::fabs(t.hess[0] - h[0]) , < , 1.e-6, "bad h");
        VERIFY_OP_ON(std::fabs(t.hess[1] - h[1]) , < , 1.e-6, "bad h 1");
        VERIFY_OP_ON(std::fabs(t.hess[2] - h[2]) , < , 1.e-6, "bad h 2");
        VERIFY_OP_ON(std::fabs(t.hess[3] - h[3]) , < , 1.e-6, "bad h 3");
        std::cout << "OK h= " << h[0] << " " << h[1] << std::endl;
        std::cout << "OK h= " << h[2] << " " << h[3] << std::endl;
      }

      void operator()(F& f, const double *u, double * h_col_major) const
      {
        double eps = f.eps();
        unsigned n = f.get_n();
        double v[n];
        double d[n]; // stepsize
        for (unsigned i = 0; i < n; ++i)
          {
            v[i] = u[i];
            d[i] = eps*(u[i] < 0.0 ? -1.0 : 1.0);
            double du = u[i] + d[i];
            d[i] = u[i] - du;
          }
        for (unsigned i = 0; i < n; ++i)
          {
            for (unsigned j = 0; j < n; ++j)
              {
                v[i] = u[i];
                v[j] = u[j];
                v[i] += d[i];
                v[j] += d[j];
                double fpp = f(v);

                v[i] = u[i];
                v[j] = u[j];
                v[i] -= d[i];
                v[j] += d[j];
                double fmp = f(v);

                v[i] = u[i];
                v[j] = u[j];
                v[i] += d[i];
                v[j] -= d[j];
                double fpm = f(v);

                v[i] = u[i];
                v[j] = u[j];
                v[i] -= d[i];
                v[j] -= d[j];
                double fmm = f(v);

                //h_col_major[i * n + j] = (fpp - fmp - fpm + fmm)/(4.0*eps*eps);
                h_col_major[i * n + j] = ((fpp+fmm) - (fmp + fpm))/(4.0*d[i]*d[j]);
                if (h_col_major[i * n + j] != h_col_major[i * n + j])
                  {
                    f.m_debug = true;
                    double fm1 = f(v);
                    std::cout << "nan: fpp= " << fpp << " fmp= " << fmp << " fpm= " << fpm << " fmm= " << fmm << " eps= " << eps << " fm1= " << fm1 << std::endl;
                  }
                v[i] = u[i];
                v[j] = u[j];

                if (0)
                  {
                    double fm2 = f(v);
                    std::cout << std::setprecision(20);
                    std::cout << "i,j= " << i << " " << j << "u= " << u[0] << " " << u[1] << " fpp= " << fpp << " fmp= " << fmp << " fpm= " << fpm << " fmm= " << fmm
                              << " fm2= " << fm2
                              << " eps= " << eps << std::endl;
                    std::cout << std::setprecision(6);
                  }

              }
          }
        // symmetrize
        if (1)
        for (unsigned i = 1; i < n; ++i)
          {
            for (unsigned j = 0; j < i; ++j)
              {
                double hh = 0.5*(h_col_major[i * n + j] + h_col_major[j * n + i] );
                h_col_major[i * n + j] = hh;
                h_col_major[j * n + i] = hh;
              }
          }
      }
    };


    template<class F>
    class FDGradWrapper
    {
    public:
      F& m_f;
      int m_index;
      FDGradWrapper(F& f) : m_f(f), m_index(0) {}

      double eps() { return m_f.eps(); }
      unsigned get_n() { return m_f.get_n(); }
      double operator()(double *v)
      {
        unsigned n = get_n();
        double grad[n];
        m_f.analytic_gradient(v, grad);
        return grad[m_index];
      }
    };

    // finite-difference Hessian of given function F when analytic gradient is known
    template<class F>
    class FDHessianFromGrad {
    public:

      F& m_f;
      typedef FDGradWrapper<F> Wrapper;
      Wrapper m_g;
      FDGradient<Wrapper> m_h;

      FDHessianFromGrad(F& f) : m_f(f) , m_g(f) {}

      void test()
      {
        FDTest t;
        FDHessianFromGrad<FDTest> fd(t);
        double v[2] = {0.1, 0.2};
        double h[4];
        fd(t, v, h);
        t(v);
        if (std::fabs(t.hess[0] - h[0])  >  1.e-6) std::cout << "t.hess[0] = " << t.hess[0] << " h[0]= " << h[0] << std::endl;

        VERIFY_OP_ON(std::fabs(t.hess[0] - h[0]) , < , 1.e-6, "bad h");
        VERIFY_OP_ON(std::fabs(t.hess[1] - h[1]) , < , 1.e-6, "bad h 1");
        VERIFY_OP_ON(std::fabs(t.hess[2] - h[2]) , < , 1.e-6, "bad h 2");
        VERIFY_OP_ON(std::fabs(t.hess[3] - h[3]) , < , 1.e-6, "bad h 3");
        std::cout << "OK h= " << h[0] << " " << h[1] << std::endl;
        std::cout << "OK h= " << h[2] << " " << h[3] << std::endl;
      }

      void operator()(F& f, const double *u, double * h_col_major)
      {
        unsigned n = m_f.get_n();
        double grad[n];
        for (unsigned i = 0; i < n; ++i)
          {
            m_g.m_index = i;
            m_h(m_g, u, grad);
            for (unsigned j = 0; j < n; ++j)
              {
                h_col_major[i * n + j] = grad[j];
              }
          }
        // symmetrize
        if (1)
          for (unsigned i = 1; i < n; ++i)
            {
              for (unsigned j = 0; j < i; ++j)
                {
                  double hh = 0.5*(h_col_major[i * n + j] + h_col_major[j * n + i] );
                  h_col_major[i * n + j] = hh;
                  h_col_major[j * n + i] = hh;
                }
            }
      }

    };

    // template<typename MeshType>
    // std::string field_name(typename MeshType::MTField *field);

    // template<>
    // std::string field_name<STKMesh>(typename STKMesh::MTField *field) {return field->name(); }


    class HexMeshSmootherMetric {
    const double m_BETA_MULT;
public:
    bool m_use_ref_mesh;
    bool m_untangling;

    HexMeshSmootherMetric(PerceptMesh * /*eMesh*/) :
            m_BETA_MULT(0.05),  m_use_ref_mesh(true), m_untangling(true) {}
                                //boolean memebrs will get altered by smoother
//private:
public:
    KOKKOS_INLINE_FUNCTION
    Double metric(double v_i_current[8][3],
    double v_i_org[8][3], bool& valid) const {

        valid = true;
        double nodal_A[8], nodal_W[8];
        Kokkos::Array<double[3][3], 8> J_A;
        Kokkos::Array<double[3][3], 8> J_W;

        sGridJacobianUtil(nodal_A, v_i_current,
                J_A);
        sGridJacobianUtil(nodal_W, v_i_org,
                J_W);

        if (m_untangling) {
            double val_untangle = 0.0;
            for (int i = 0; i < 8; i++) {

                double detAi = nodal_A[i];
                double detWi = nodal_W[i];

                if (detAi <= 0.)
                    valid = false;

                double vv = m_BETA_MULT * detWi - detAi;
                (vv < 0.0 ? vv = 0.0 : vv);
                val_untangle += vv * vv;
            }
            return val_untangle;
        }

        else{//smoothing
            double val = 0.0, val_shape = 0.0;

            double T[3][3];
            double WI[3][3];

            for (int i = 0; i < 8; i++) {
                const double detAi = nodal_A[i];
                double detWi = nodal_W[i];
                const double detWiC = detWi;

                if (detAi <= 0) {
                    valid = false;
                }

                const double * W[3]=
                        { J_W[i][0], J_W[i][1], J_W[i][2] };
                const double * A[3]=
                        { J_A[i][0], J_A[i][1], J_A[i][2] };

                double shape_metric = 0.0;
                if (detAi > 0.0) {

                    if (m_use_ref_mesh) {
                        matrix_inverse(W, detWi, WI);
                        matrix_product(A, WI, T);
                    } else {
                        T[0][0] = A[0][0];
                        T[0][1] = A[0][1];
                        T[0][2] = A[0][2];

                        T[1][0] = A[1][0];
                        T[1][1] = A[1][1];
                        T[1][2] = A[1][2];

                        T[2][0] = A[2][0];
                        T[2][1] = A[2][1];
                        T[2][2] = A[2][2];

                        //T *= 1.0/std::pow(detWi, 1./3.);
                        detWi = 1.0;
                    }
                    double d = detAi / detWi;
                    double f = matrix_sqr_Frobenius(T);

                    // || T ||^3 / (3 sqrt(3) det(T)) - 1.0
//                       if (spatialDim == 2)
//                         {
//                           // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
//                           f = f - 1.0;
//                           VERIFY_OP_ON(f, >, 0.0, "bad f");
//                           f = std::sqrt(f);
//                           double fac = 2.0;
//                           double den = fac * d;
//                           shape_metric = (f*f)/den - 1.0;
//                           shape_metric = shape_metric*detWiC;
//                           if (Base::m_debug)
//                             {
//                               std::cout << "nan: shape_metric= " << shape_metric << " f= " << f << " den= " << den << " detWi= " << detWi << " detAi= " << detAi << std::endl;
//                             }
//                         }
//                       else
                    {
                        f = std::sqrt(f);
                        const double fac = 3.0 * std::sqrt(3.0);
                        double den = fac * d;
                        shape_metric = (f * f * f) / den - 1.0;
                        shape_metric = shape_metric * detWiC;
                    }
                }
                val_shape += shape_metric;
            }
            val = val_shape;
            return val;
        }
        return 0.0;
    }//metric

    KOKKOS_INLINE_FUNCTION
    double grad_metric(  double v_i_current[8][3],
            double v_i_org[8][3], bool& valid, double grad[8][4]) const
    {
        const unsigned locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                4, 2 }, { 7, 5, 6, 3 } };

        valid = true;
        double nodal_A[8], nodal_W[8];
        Kokkos::Array<double[3][3], 8> J_A;
        Kokkos::Array<double[3][3], 8> J_W;
        double dMetric_dA[8][3][3];

        sGridJacobianUtil(nodal_A, v_i_current,
                J_A);
        sGridJacobianUtil(nodal_W, v_i_org,
                J_W);


        if(m_untangling)
        {
            double val = 0.0;

            for (int i=0; i < 8; i++)
              {
                double detAi = nodal_A[i];
                double detWi = nodal_W[i];
                if (detAi <= 0)
                  {
                    valid = false;
                  }

                double vv = m_BETA_MULT*detWi - detAi;
                double vv1 = (vv>0.0 ? vv : 0);

                if (vv > 0.0)
                  {
                    transpose_adj_dev(J_A[i],dMetric_dA[i],-2.0*vv);

//                    if (Base::m_eMesh->get_spatial_dim() == 2)
//                      {
//                        for (unsigned ii=0; ii < 3; ++ii)
//                          {
//                            wrt_A(2,ii) = 0.0;
//                            wrt_A(ii,2) = 0.0;
//                          }
//                      }
                  }
                else{
                    for(int j=0;j<3;j++)
                        for(int k=0;k<3;k++)
                        dMetric_dA[i][j][k]=0.0;
                }
                val += vv1*vv1;
              }


            double grads_fld[8][8][3];
            grad_metric_util(dMetric_dA,grads_fld ,locs_hex_dev);



            // combine into total
            for (unsigned i=0; i < 8; i++)
              for (unsigned j = 0; j < 3; j++)
                grad[i][j] = 0.0;

            for (unsigned k=0; k < 8; k++)
              for (unsigned i=0; i < 8; i++)
                for (unsigned j = 0; j < 3; j++){
                  grad[i][j] += grads_fld[k][i][j];

                }

            return val;
        }//untangling

        else{//SMOOTHING

            double WIT[3][3];
            double TAT[3][3];
            double T[3][3];
            double WI[3][3];

            double val=0.0, val_shape=0.0;

                   for (int i=0; i < 8; i++)
                     {
                       double detAi = nodal_A[i];
                       double detWi = nodal_W[i];
                       double detWiC = detWi;

                       if (detAi <= 0)
                         {
                           valid = false;
                         }
                       const double * W[3]=
                               { J_W[i][0], J_W[i][1], J_W[i][2] };

                       const double * A[3]=
                               { J_A[i][0], J_A[i][1], J_A[i][2] };

                       // frob2 = h^2 + h^2 + 1
                       // frob21 = 2 h^2
                       // f = h sqrt(2)
                       // det = h*h
                       // met = f*f / (det*2) - 1
                       // frob3 = 3 h^2
                       // f = h sqrt(3)
                       // det = h*h*h
                       // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
                       double shape_metric = 0.0;
                       if (detAi > 0) //! 1.e-15)
                         {
                           if (m_use_ref_mesh)
                             {
                               matrix_inverse(W, detWi, WI);
                               matrix_product(A, WI, T);
                             }
                           else
                             {
                               T[0][0] = A[0][0];
                               T[0][1] = A[0][1];
                               T[0][2] = A[0][2];

                               T[1][0] = A[1][0];
                               T[1][1] = A[1][1];
                               T[1][2] = A[1][0];

                               T[2][0] = A[2][0];
                               T[2][1] = A[2][1];
                               T[2][2] = A[2][2];
                               detWi = 1.0;
                               identity_dev(WI);
                             }

                           double d = detAi / detWi;
                           double f = matrix_sqr_Frobenius(T);
//                           double f0 = f;

//                           if (spatialDim == 2)
//                             {
//                               // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
//                               f = f - 1.0;
//                               f = std::sqrt(f);
//                               double fac = 2.0;
//                               double den = fac * d;
//                               shape_metric = (f*f)/den - 1.0;
//                               DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
//                               {
//                                 wrt_A = 0.5*(2.0*T  / d);
//                                 wrt_A -= 0.5*(f0 - 1.0) / (d*d)   * transpose_adj(T);
//
//                                 // now convert to wrt_A
//                                 if (m_use_ref_mesh)
//                                   wrt_A = wrt_A * transpose(WI);
//
//                                 for (unsigned ii=0; ii < 3; ++ii)
//                                   {
//                                     wrt_A(2,ii) = 0.0;
//                                     wrt_A(ii,2) = 0.0;
//                                   }
//
//                                 wrt_A = wrt_A * detWiC;
//                               }
//                             }
//                           else
                             {
                               f = std::sqrt(f);
                               const double fac = 3.0*std::sqrt(3.0);
                               double den = fac * d;
                               shape_metric = (f*f*f)/den - 1.0;
                               {
                                 double norm = f;
                                 double iden = 1.0/den;
                                 //result = norm_cube * iden - 1.0;
                                 // wrt_T...
                                 dMetric_dA[i][0][0] = (3 * norm * iden)*T[0][0];
                                 dMetric_dA[i][0][1] = (3 * norm * iden)*T[0][1];
                                 dMetric_dA[i][0][2] = (3 * norm * iden)*T[0][2];

                                 dMetric_dA[i][1][0] = (3 * norm * iden)*T[1][0];
                                 dMetric_dA[i][1][1] = (3 * norm * iden)*T[1][1];
                                 dMetric_dA[i][1][2] = (3 * norm * iden)*T[1][2];

                                 dMetric_dA[i][2][0] = (3 * norm * iden)*T[2][0];
                                 dMetric_dA[i][2][1] = (3 * norm * iden)*T[2][1];
                                 dMetric_dA[i][2][2] = (3 * norm * iden)*T[2][2];

                                 double norm_cube = norm*norm*norm;
                                 transpose_adj_dev(T, TAT);
                                 dMetric_dA[i][0][0] -= (norm_cube * iden/d) * TAT[0][0];
                                 dMetric_dA[i][0][1] -= (norm_cube * iden/d) * TAT[0][1];
                                 dMetric_dA[i][0][2] -= (norm_cube * iden/d) * TAT[0][2];

                                 dMetric_dA[i][1][0] -= (norm_cube * iden/d) * TAT[1][0];
                                 dMetric_dA[i][1][1] -= (norm_cube * iden/d) * TAT[1][1];
                                 dMetric_dA[i][1][2] -= (norm_cube * iden/d) * TAT[1][2];

                                 dMetric_dA[i][2][0] -= (norm_cube * iden/d) * TAT[2][0];
                                 dMetric_dA[i][2][1] -= (norm_cube * iden/d) * TAT[2][1];
                                 dMetric_dA[i][2][2] -= (norm_cube * iden/d) * TAT[2][2];

                                 // now convert to dMetric_dA[i]
                                 if (m_use_ref_mesh)
                                   {
                                     matrix_transpose(WI, WIT);
                                     matrix_product_mutator(dMetric_dA[i],WIT);
                                   }
                                 dMetric_dA[i][0][0] = dMetric_dA[i][0][0] * detWiC;
                                 dMetric_dA[i][0][1] = dMetric_dA[i][0][1] * detWiC;
                                 dMetric_dA[i][0][2] = dMetric_dA[i][0][2] * detWiC;

                                 dMetric_dA[i][1][0] = dMetric_dA[i][1][0] * detWiC;
                                 dMetric_dA[i][1][1] = dMetric_dA[i][1][1] * detWiC;
                                 dMetric_dA[i][1][2] = dMetric_dA[i][1][2] * detWiC;

                                 dMetric_dA[i][2][0] = dMetric_dA[i][2][0] * detWiC;
                                 dMetric_dA[i][2][1] = dMetric_dA[i][2][1] * detWiC;
                                 dMetric_dA[i][2][2] = dMetric_dA[i][2][2] * detWiC;
                               }
                             }
                         }//if (detAi > 0)
                       val_shape += shape_metric;
                     }//foreach node

                   // compute grad for all nodes
                   double grads_fld[8][8][3];
                   grad_metric_util(dMetric_dA,grads_fld ,locs_hex_dev);

                   // combine into total
                   for (int i=0; i < 8; i++)
                     for (int j = 0; j < 3; j++)
                       grad[i][j] = 0.0;

                   for (int k=0; k < 8; k++)
                     for (int i=0; i < 8; i++)
                       for (int j = 0; j < 3; j++)
                         grad[i][j] += grads_fld[k][i][j];

                   val = val_shape;
                   return val;
        }//SMOOTHING
        return 0.0;
    }
};//HexMeshSmootherMetric

    template<typename MeshType>
    class SmootherMetricUntangleImpl : public SmootherMetricImpl<MeshType>
    {
    protected:
      double m_beta_mult;

    public:

      using Base = SmootherMetricImpl<MeshType>;

      SmootherMetricUntangleImpl(PerceptMesh *eMesh);

      virtual bool has_gradient() override { return true; }

      virtual double length_scaling_power() override { return 3.0; }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(typename MeshType::MTElement element, bool& valid) override;

      /// computes metric and its gradient - see Mesquite::TShapeB1, TQualityMetric, TargetMetricUtil
      virtual double grad_metric(typename MeshType::MTElement element, bool& valid, double grad[8][4]) override
      {
        JacobianUtilImpl<MeshType> jacA, jacW;
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);
        jacW(W_, *Base::m_eMesh, element, Base::m_coord_field_original, Base::m_topology_data);
        double val=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            double vv = m_beta_mult*detWi - detAi;
            double vv1 = std::max(vv, 0.0);
            DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
            if (vv > 0.0)
              {
                wrt_A = -2.0*vv*transpose_adj(A);

                if (Base::m_eMesh->get_spatial_dim() == 2)
                  {
                    for (unsigned ii=0; ii < 3; ++ii)
                      {
                        wrt_A(2,ii) = 0.0;
                        wrt_A(ii,2) = 0.0;
                      }
                  }
              }
            else
              wrt_A.zero();
            val += vv1*vv1;
          }

        // compute grad for all nodes
        jacA.grad_metric_util( *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);

        int spatialDim = Base::m_eMesh->get_spatial_dim();

        // combine into total
        for (int i=0; i < nn; i++)
          for (int j = 0; j < spatialDim; j++)
            grad[i][j] = 0.0;

        for (int k=0; k < nn; k++)
          for (int i=0; i < nn; i++)
            for (int j = 0; j < spatialDim; j++)
              grad[i][j] += jacA.m_grad[k][i][j];

        return val;
      }

    };
    using SmootherMetricUntangle = SmootherMetricUntangleImpl<STKMesh>;

    class SmootherMetricShapeSizeOrient : public SmootherMetric
    {
    public:
      SmootherMetricShapeSizeOrient(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}
      virtual double length_scaling_power() override { return 1.0; }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA, jacW;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        DenseMatrix<3,3> Ident;
        identity(Ident);

        DenseMatrix<3,3> AI, Atmp, WAI;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            if (std::fabs(detAi) > 0) //! 1.e-10)
              {
                DenseMatrix<3,3>& W = jacW.m_J[i];
                DenseMatrix<3,3>& A = jacA.m_J[i];
                inverse(A, AI);
                product(W, AI, WAI);
                difference(WAI, Ident, Atmp);
                shape_metric = my_sqr_Frobenius(Atmp);
              }
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }
    };

    class SmootherMetricScaledJacobianNodal : public SmootherMetric
    {
    public:
      SmootherMetricScaledJacobianNodal(PerceptMesh *eMesh) : SmootherMetric(eMesh)
      { m_is_nodal=true; m_combine=COP_MAX; }
      virtual double length_scaling_power() override { return 1.0; }

//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        VERIFY_OP_ON(m_node, !=, stk::mesh::Entity(), "must set a node");
        valid = true;
        JacobianUtil jacA, jacSA, jacW;

        double SA_ = 0.0, A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacSA(SA_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        const MyPairIterRelation elem_nodes(*m_eMesh, element, m_eMesh->node_rank());
        VERIFY_OP_ON((int)elem_nodes.size(), ==, jacA.m_num_nodes, "node num mismatch");
        val_shape = 0.0;
        bool found = false;
        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            if (elem_nodes[i].entity() == m_node)
              {
                double detAi = jacA.m_detJ[i];
                double detSAi = jacSA.m_detJ[i];
                double detWi = jacW.m_detJ[i];
                if (detAi <= 0)
                  {
                    valid = false;
                  }
                double shape_metric = 0.0;
                //DenseMatrix<3,3>& A = jacA.m_J[i];
                double scale_factor = detWi;
                scale_factor = 1.0;
                double fac = 0.2;
                shape_metric = scale_factor* (detSAi > fac ? fac : detSAi);
                val_shape = shape_metric;
                //std::cout << "tmp srk i= " << i << " detAi = " << detAi << " detSAi= " << detSAi << " shape_metric= " << shape_metric << " val_shape= " << val_shape << std::endl;
                found = true;
                break;
              }
          }
        VERIFY_OP_ON(found, ==, true, "logic err");
        val = -val_shape;
        //std::cout << "tmp srk val = " << val << std::endl;
        return val;
      }
    };

    class SmootherMetricScaledJacobianElemental : public SmootherMetric
    {
    public:
      SmootherMetricScaledJacobianElemental(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}
      virtual double length_scaling_power() override { return 1.0; }

      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA, jacSA, jacW;

        double SA_ = 0.0, A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacSA(SA_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        val_shape = 0.0;
        //val_shape = std::numeric_limits<double>::max();
        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detSAi = jacSA.m_detJ[i];
            //double detWi = jacW.m_detJ[i];
            double FAi = Frobenius(jacA.m_J[i]);
            double FWi = Frobenius(jacW.m_J[i]);
            if (detAi <= 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            //DenseMatrix<3,3>& A = jacA.m_J[i];
            double scale_factor = (FAi < FWi ? FAi/FWi : (FAi > 1.e-6? FWi / FAi : 1.0));
            scale_factor = FAi/FWi;
            //scale_factor = 1.0;
            //double sign_SA = (detSAi > 0.0 ? 1.0 : -1.0);
            //double fac = 0.2;
            //shape_metric = scale_factor* (detSAi > fac ? fac : detSAi);
            shape_metric = scale_factor*detSAi;
            val_shape += shape_metric;
            //val_shape = std::min(val_shape, shape_metric);
            //std::cout << "tmp srk i= " << i << " detAi = " << detAi << " detSAi= " << detSAi << " shape_metric= " << shape_metric << " val_shape= " << val_shape << " scale_factor= " << scale_factor << " FAi= " << FAi << " FWi= " << FWi << std::endl;
          }

        val = -val_shape;
        //std::cout << "tmp srk val = " << val << std::endl;
        return val;
      }
    };

    template<typename MeshType>
    class SmootherMetricShapeB1Impl : public SmootherMetricImpl<MeshType>
    {
    protected:
      DenseMatrix<3,3> Ident;
      DenseMatrix<3,3> WI, T, WIT, TAT;
      const int spatialDim;
      bool m_use_ref_mesh;
    public:
      using Base = SmootherMetricImpl<MeshType>;
      SmootherMetricShapeB1Impl(PerceptMesh *eMesh, bool use_ref_mesh=true) : Base(eMesh), spatialDim( Base::m_eMesh->get_spatial_dim()), m_use_ref_mesh(use_ref_mesh) {
        identity(Ident);
      }

      virtual bool has_gradient() override { return true; }

      virtual double length_scaling_power() override { return 1.0; }
      virtual double metric(typename MeshType::MTElement element, bool& valid) override
      {
        valid = true;

        JacobianUtilImpl<MeshType> jacA, jacW;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);
        jacW(W_, *Base::m_eMesh, element, Base::m_coord_field_original, Base::m_topology_data);
        double val=0.0, val_shape=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            const double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            VERIFY_OP_ON(detWi, >, 0.0, "bad reference mesh");
            const double detWiC = detWi;
//            bool isNaN = false;
//            bool isInf = false;
            if (detAi <= 0)
              {
                valid = false;
              }
            // else
            //   {
            //     isNaN = (1.0/detAi != 1.0/detAi);
            //     isInf = std::isinf(1.0/detAi);
            //     if (isNaN || isInf)
            //       valid = false;
            //   }

            const DenseMatrix<3,3>& W = jacW.m_J[i];
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            // frob2 = h^2 + h^2 + 1
            // frob21 = 2 h^2
            // f = h sqrt(2)
            // det = h*h
            // met = f*f / (det*2) - 1
            // frob3 = 3 h^2
            // f = h sqrt(3)
            // det = h*h*h
            // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
            double shape_metric = 0.0;
            if (detAi > 0.0 /*&& !isNaN && !isInf*/)
              {

                if (m_use_ref_mesh)
                  {
                    inverse(W, detWi, WI);
                    product(A, WI, T);
                  }
                else
                  {
                    T = A;
                    //T *= 1.0/std::pow(detWi, 1./3.);
                    detWi = 1.0;
                  }
                double d = detAi/detWi;
                double f = my_sqr_Frobenius(T);

                // || T ||^3 / (3 sqrt(3) det(T)) - 1.0
                if (spatialDim == 2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;
                    VERIFY_OP_ON(f, >, 0.0, "bad f");
                    f = std::sqrt(f);
                    double fac = 2.0;
                    double den = fac * d;
                    shape_metric = (f*f)/den - 1.0;
                    shape_metric = shape_metric*detWiC;
                    if (Base::m_debug)
                      {
                        std::cout << "nan: shape_metric= " << shape_metric << " f= " << f << " den= " << den << " detWi= " << detWi << " detAi= " << detAi << std::endl;
                      }
                  }
                else
                  {
                    f = std::sqrt(f);
                    const double fac = 3.0*std::sqrt(3.0);
                    double den = fac * d;
                    shape_metric = (f*f*f)/den - 1.0;
                    shape_metric = shape_metric*detWiC;
                  }
              }
            val_shape += shape_metric;
            if (Base::m_debug)
              {
                std::cout << "0nan: shape_metric= " << shape_metric << " detWi= " << detWi << " detAi= " << detAi << std::endl;
              }
          }
        val = val_shape;
        return val;
      }//metric

      /// computes metric and its gradient - see Mesquite::TShapeB1, TQualityMetric, TargetMetricUtil
      virtual double grad_metric(typename MeshType::MTElement element, bool& valid, double grad[8][4]) override
      {
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        JacobianUtilImpl<MeshType> jacA, jacW;
        jacA(A_, *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);
        jacW(W_, *Base::m_eMesh, element, Base::m_coord_field_original, Base::m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            double detWiC = detWi;

            if (detAi <= 0)
              {
                valid = false;
              }
            const DenseMatrix<3,3>& W = jacW.m_J[i];
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            // frob2 = h^2 + h^2 + 1
            // frob21 = 2 h^2
            // f = h sqrt(2)
            // det = h*h
            // met = f*f / (det*2) - 1
            // frob3 = 3 h^2
            // f = h sqrt(3)
            // det = h*h*h
            // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
            double shape_metric = 0.0;
            if (detAi > 0) //! 1.e-15)
              {
                if (m_use_ref_mesh)
                  {
                    inverse(W, detWi, WI);
                    product(A, WI, T);
                  }
                else
                  {
                    T = A;
                    detWi = 1.0;
                    identity(WI);
                  }

                double d = detAi / detWi;
                double f = my_sqr_Frobenius(T);
                double f0 = f;

                if (spatialDim == 2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;
                    f = std::sqrt(f);
                    double fac = 2.0;
                    double den = fac * d;
                    shape_metric = (f*f)/den - 1.0;
                    DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
                    {
                      wrt_A = 0.5*(2.0*T  / d);
                      wrt_A -= 0.5*(f0 - 1.0) / (d*d)   * transpose_adj(T);

                      // now convert to wrt_A
                      if (m_use_ref_mesh)
                        wrt_A = wrt_A * transpose(WI);

                      for (unsigned ii=0; ii < 3; ++ii)
                        {
                          wrt_A(2,ii) = 0.0;
                          wrt_A(ii,2) = 0.0;
                        }

                      wrt_A = wrt_A * detWiC;
                    }
                  }
                else
                  {
                    f = std::sqrt(f);
                    const double fac = 3.0*std::sqrt(3.0);
                    double den = fac * d;
                    shape_metric = (f*f*f)/den - 1.0;

                    DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
                    {
                      double norm = f;
                      double iden = 1.0/den;
                      //result = norm_cube * iden - 1.0;
                      // wrt_T...
                      wrt_A = (3 * norm * iden)*T;

                      double norm_cube = norm*norm*norm;
                      transpose_adj(T, TAT);
                      wrt_A -= (norm_cube * iden/d) * TAT;

                      // now convert to wrt_A
                      if (m_use_ref_mesh)
                        {
                          transpose(WI, WIT);
                          wrt_A = wrt_A * WIT;
                        }
                      wrt_A = wrt_A * detWiC;
                    }
                  }
              }
            val_shape += shape_metric;
          }

        // compute grad for all nodes
        jacA.grad_metric_util( *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);

        // combine into total
        for (int i=0; i < nn; i++)
          for (int j = 0; j < 3; j++)
            grad[i][j] = 0.0;

        for (int k=0; k < nn; k++)
          for (int i=0; i < nn; i++)
            for (int j = 0; j < spatialDim; j++)
              grad[i][j] += jacA.m_grad[k][i][j];

        val = val_shape;
        return val;
      }

    };

    using SmootherMetricShapeB1 = SmootherMetricShapeB1Impl<STKMesh>;

    /* 2d: Frob^2(T) / 2*(|T|) - 1 */
    /* 3d: Frob^2(T) / 3*(|T|^2/3) - 1 */
    class SmootherMetricShapeMeanRatio : public SmootherMetric
    {
    protected:
      DenseMatrix<3,3> Ident;
      DenseMatrix<3,3> WI, T, WIT, TAT;
      const int spatialDim;
      bool m_use_ref_mesh;
    public:
      SmootherMetricShapeMeanRatio(PerceptMesh *eMesh, bool use_ref_mesh=true) : SmootherMetric(eMesh), spatialDim( m_eMesh->get_spatial_dim()), m_use_ref_mesh(use_ref_mesh) {
        identity(Ident);
      }

      virtual bool has_gradient() override { return true; }

      virtual double length_scaling_power() override { return 1.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;

        JacobianUtil jacA, jacW;
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            const double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            VERIFY_OP_ON(detWi, >, 0.0, "bad reference mesh");
            const double detWiC = 1.0; //detWi;
            bool isNaN = false;
            bool isInf = false;
            if (detAi <= 0)
              {
                valid = false;
              }

            const DenseMatrix<3,3>& W = jacW.m_J[i];
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            // 2d:  result = 0.5 * sqr_Frobenius(T) / d - 1;
            // 3d:  result = sqr_Frobenius(T) / (3 * det_cbrt * det_cbrt) - 1;

            double shape_metric = 0.0;
            if (detAi > 0.0 && !isNaN && !isInf)
              {
                if (m_use_ref_mesh)
                  {
                    inverse(W, detWi, WI);
                    product(A, WI, T);
                  }
                else
                  {
                    T = A;
                    detWi = 1.0;
                  }
                double d = detAi/detWi;
                double f = my_sqr_Frobenius(T);

                // || T ||^3 / (3 sqrt(3) det(T)) - 1.0
                if (spatialDim == 2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;
                    VERIFY_OP_ON(f, >, 0.0, "bad f");
                    shape_metric = 0.5*f / d - 1.0;
                    shape_metric = shape_metric*detWiC;
                    if (m_debug)
                      {
                        std::cout << "nan: shape_metric= " << shape_metric << " f= " << f << " d= " << d << " detWi= " << detWi << " detAi= " << detAi << std::endl;
                      }
                  }
                else
                  {
                    shape_metric = f / (3 * std::pow(d, 2./3.)) - 1.0;
                    shape_metric = shape_metric*detWiC;
                  }
              }
            val_shape += shape_metric;
            if (m_debug)
              {
                std::cout << "0nan: shape_metric= " << shape_metric << " detWi= " << detWi << " detAi= " << detAi << std::endl;
              }
          }
        val = val_shape;
        return val;
      }

      /// computes metric and its gradient - see Mesquite::TShapeB1, TQualityMetric, TargetMetricUtil
      virtual double grad_metric(stk::mesh::Entity element, bool& valid, double grad[8][4]) override
      {
        JacobianUtil jacA, jacW;
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            double detWiC = 1.0; //detWi;
            if (detAi <= 0)
              {
                valid = false;
              }
            const DenseMatrix<3,3>& W = jacW.m_J[i];
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            // frob2 = h^2 + h^2 + 1
            // frob21 = 2 h^2
            // f = h sqrt(2)
            // det = h*h
            // met = f*f / (det*2) - 1
            // frob3 = 3 h^2
            // f = h sqrt(3)
            // det = h*h*h
            // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
            double shape_metric = 0.0;
            if (detAi > 0) //! 1.e-15)
              {
                if (m_use_ref_mesh)
                  {
                    inverse(W, detWi, WI);
                    product(A, WI, T);
                  }
                else
                  {
                    T = A;
                    detWi = 1.0;
                    identity(WI);
                  }

                double d = detAi / detWi;
                double f = my_sqr_Frobenius(T);
                //double f0 = f;

                if (spatialDim == 2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;

                    shape_metric = 0.5*f / d - 1.0;
                    //shape_metric = (f*f)/den - 1.0;
                    DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
                    {
                      wrt_A = (T  / d) ;
                      wrt_A -=  f / (2.0*d*d) * transpose_adj(T);

                      // now convert to wrt_A
                      if (m_use_ref_mesh)
                        wrt_A = wrt_A * transpose(WI);

                      for (unsigned ii=0; ii < 3; ++ii)
                        {
                          wrt_A(2,ii) = 0.0;
                          wrt_A(ii,2) = 0.0;
                        }

                      wrt_A = wrt_A * detWiC;
                    }
                  }
                else
                  {
                    //f = std::sqrt(f);
                    //const double fac = 3.0*std::sqrt(3.0);
                    //double den = fac * d;
                    //shape_metric = (f*f*f)/den - 1.0;
                    shape_metric = f / (3 * std::pow(d, 2./3.)) - 1.0;

                    DenseMatrix<3,3>& wrt_A = jacA.m_dMetric_dA[i];
                    {
                      wrt_A = 2.0*T/(3.0*std::pow(d,2./3.));

                      transpose_adj(T, TAT);
                      wrt_A -= (2./3.*std::pow(d,-5./3.)) * TAT;

                      // now convert to wrt_A
                      if (m_use_ref_mesh)
                        {
                          transpose(WI, WIT);
                          wrt_A = wrt_A * WIT;
                        }
                      wrt_A = wrt_A * detWiC;
                    }
                  }
              }
            val_shape += shape_metric;
          }

        // compute grad for all nodes
        jacA.grad_metric_util( *m_eMesh, element, m_coord_field_current, m_topology_data);

        // combine into total
        for (int i=0; i < nn; i++)
          for (int j = 0; j < 3; j++)
            grad[i][j] = 0.0;

        for (int k=0; k < nn; k++)
          for (int i=0; i < nn; i++)
            for (int j = 0; j < spatialDim; j++)
              grad[i][j] += jacA.m_grad[k][i][j];

        val = val_shape;
        return val;
      }

    };


    /** |T-I|^2/ (2 det(T)) */

    class SmootherMetricShapeSizeOrientB1 : public SmootherMetric
    {
    public:
      SmootherMetricShapeSizeOrientB1(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}

      virtual double length_scaling_power() override { return 1.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA, jacW;

        //int spatialDim = m_eMesh->get_spatial_dim();
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        DenseMatrix<3,3> Ident;
        identity(Ident);

        DenseMatrix<3,3> WI, T;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            DenseMatrix<3,3>& W = jacW.m_J[i];
            DenseMatrix<3,3>& A = jacA.m_J[i];

            double shape_metric = 0.0;
            if (std::fabs(detAi) > 0) //!1.e-15)
              {
                inverse(W, WI);
                product(A, WI, T);

                /** |T-I|^2/ (2 det(T)) */

                double d = det(T);
                double f = my_sqr_Frobenius(T-Ident)/(2*d);
                shape_metric = f*det(W);

#if 0
                double n = Frobenius(T);
                double tau = d;
                shape_metric = n*n*n - 3*MSQ_SQRT_THREE*( std::log(tau) + 1 );

                n = Frobenius(A - W);
                shape_metric = n*n*n - det(W)*3*MSQ_SQRT_THREE*( std::log(detAi) + 1 );

                shape_metric = my_sqr_Frobenius( A - 1/detAi * transpose_adj(A) * transpose(W) * W );
#endif
              }
            val_shape += shape_metric;
            //val_shape += std::fabs(shape_metric);
            //val_shape += shape_metric*shape_metric;
          }
        val = val_shape;
        //val = val_shape*val_shape;
        return val;
      }


    };

    class SmootherMetricShapeC1 : public SmootherMetric
    {
    public:
      SmootherMetricShapeC1(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}

      virtual double length_scaling_power() override { return 1.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA, jacW;

        int spatialDim = m_eMesh->get_spatial_dim();
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        DenseMatrix<3,3> Ident;
        identity(Ident);

        DenseMatrix<3,3> WI, T;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                //valid = false;
              }
            DenseMatrix<3,3>& W = jacW.m_J[i];
            DenseMatrix<3,3>& A = jacA.m_J[i];

            // frob2 = h^2 + h^2 + 1
            // frob21 = 2 h^2
            // f = h sqrt(2)
            // det = h*h
            // met = f*f / (det*2) - 1
            // frob3 = 3 h^2
            // f = h sqrt(3)
            // det = h*h*h
            // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
            double shape_metric = 0.0;
            //            if (std::fabs(detAi) > 1.e-15)
              {
                inverse(W, WI);
                //product(A, WI, T);
                T = A;
                double d = det(T);
                double f = my_sqr_Frobenius(T);
                if (0 && spatialDim==2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;
                    f = std::sqrt(f);
                    double fac = 2.0;
                    double den = fac * d;
                    shape_metric = (f*f)/den - 1.0;
                  }
                else
                  {
                    f = std::sqrt(f);
                    double fac = 3.0*std::sqrt(3.0);
                    double den = fac * d;
                    //shape_metric = (f*f*f)/den - 1.0;
                    shape_metric = (f*f*f) - den;
                  }
                //shape_metric = std::fabs(shape_metric);
                //shape_metric = f/std::pow(den,1./3.) - 1.0;
              }
            val_shape += shape_metric;
            //val_shape += std::fabs(shape_metric);
            //val_shape += shape_metric*shape_metric;
          }
        val = val_shape;
        //val = val_shape*val_shape;
        return val;
      }


    };

    class SmootherMetricLaplace : public SmootherMetric
    {
    public:
      SmootherMetricLaplace(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}

      virtual double length_scaling_power() override { return 2.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA;
        //JacobianUtil jacW;

        double A_ = 0.0;
        //double W_ = 0.0;
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            DenseMatrix<3,3>& A = jacA.m_J[i];
            double shape_metric = 0.0;
            shape_metric = sqr_Frobenius(A);
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }

    };

    class SmootherMetricLaplaceInverseVolumeWeighted : public SmootherMetric
    {
    public:
      SmootherMetricLaplaceInverseVolumeWeighted(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}

      virtual double length_scaling_power() override { return 2.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA;
        //JacobianUtil jacW;

        double A_ = 0.0;
        //double W_ = 0.0;
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            DenseMatrix<3,3>& A = jacA.m_J[i];
            double shape_metric = 0.0;
            double f = sqr_Frobenius(A);
            shape_metric =  f*f*f / detAi;
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }

    };

    class SmootherMetricVolumetricEnergy : public SmootherMetric
    {
    public:
      SmootherMetricVolumetricEnergy(PerceptMesh *eMesh) : SmootherMetric(eMesh) {}

      virtual double length_scaling_power() override { return 6.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;
        JacobianUtil jacA;
        //JacobianUtil jacW;

        double A_ = 0.0;
        //double W_ = 0.0;
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi <= 0)
              {
                valid = false;
              }
            double shape_metric = detAi*detAi;
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }

    };

    class SmootherMetricSizeB1 : public SmootherMetric
    {
    protected:
      DenseMatrix<3,3> Ident;
      DenseMatrix<3,3> WI, T, WIT, TAT;
      const int spatialDim;
      bool m_use_ref_mesh;
    public:
      SmootherMetricSizeB1(PerceptMesh *eMesh, bool use_ref_mesh=true) : SmootherMetric(eMesh), spatialDim( m_eMesh->get_spatial_dim()), m_use_ref_mesh(use_ref_mesh) {
        identity(Ident);
      }

      virtual bool has_gradient() override { return false; }

      virtual double length_scaling_power() override { return 1.0; }
//KOKKOS_INLINE_FUNCTION
      virtual double metric(stk::mesh::Entity element, bool& valid) override
      {
        valid = true;

        JacobianUtil jacA, jacW;
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        const int nn = jacA.m_num_nodes;
        for (int i=0; i < nn; i++)
          {
            const double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            VERIFY_OP_ON(detWi, >, 0.0, "bad reference mesh");
            if (detAi <= 0)
              {
                valid = false;
              }
            const DenseMatrix<3,3>& W = jacW.m_J[i];
            const DenseMatrix<3,3>& A = jacA.m_J[i];

            // d + 1/d - 2.0
            double shape_metric = 0.0;
            if (detAi > 0.0)
              {
                if (m_use_ref_mesh)
                  {
                    inverse(W, detWi, WI);
                    product(A, WI, T);
                  }
                else
                  {
                    T = A;
                    detWi = 1.0;
                  }
                double d = detAi/detWi;

                // det(T) + 1/det(T) - 2.0
                if (1)
                  {
                    shape_metric = d + 1.0/d - 2.0;
                    shape_metric = shape_metric*detWi;
                    if (m_debug)
                      {
                        std::cout << "nan: shape_metric= " << shape_metric << " detWi= " << detWi << " detAi= " << detAi << std::endl;
                      }
                  }
              }
            val_shape += shape_metric*shape_metric;
            if (m_debug)
              {
                std::cout << "0nan: shape_metric= " << shape_metric << " detWi= " << detWi << " detAi= " << detAi << std::endl;
              }
          }
        val = val_shape;
        return val;
      }

    };

  }//percept

#endif
#endif
