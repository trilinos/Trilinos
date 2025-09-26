// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmootherConjugateGradientDef_hpp
#define ReferenceMeshSmootherConjugateGradientDef_hpp

#include <percept/PerceptUtils.hpp>

namespace percept {

template <typename MeshType>
ReferenceMeshSmootherConjugateGradientImpl<MeshType>::
ReferenceMeshSmootherConjugateGradientImpl(PerceptMesh *eMesh,
                                       STKMesh::MTSelector *stk_select,
                                       typename MeshType::MTMeshGeometry *meshGeometry,
                                       int inner_iterations,
                                       double grad_norm,
                                       int parallel_iterations)

  :  Base(eMesh, stk_select, meshGeometry, inner_iterations, grad_norm, parallel_iterations), m_max_edge_length_factor(1.0),
     m_coord_updater(this,m_eMesh,0),
     m_metric_computinator(this->m_metric,eMesh,0,0,0,0)
{}

  template<>
  Double
  ReferenceMeshSmootherConjugateGradientImpl< STKMesh >::
  nodal_edge_length_ave(stk::mesh::Entity node)
  {
    //int spatialDim = m_eMesh->get_spatial_dim();
    Double nm=0.0;

    double min=std::numeric_limits<double>::max();
    const MyPairIterRelation node_elems(*m_eMesh, node, m_eMesh->element_rank() );
    Double nele = 0.0;
    for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
      {
        stk::mesh::Entity element = node_elems[i_elem].entity();
        if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
          continue;
        double lmin=0,lmax=0;
        double elem_edge_len = m_eMesh->edge_length_ave(element, m_coord_field_original, &lmin, &lmax);
        if (lmin < min) min=lmin;
        nm += elem_edge_len;
        nele += 1.0;
      }
    nm /= nele;
    if (0)
      return min;
    return nm;
  }

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  debug_print(double /*alpha*/)
  {
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  snap_nodes
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_snap_nodes
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using This = GenericAlgorithm_snap_nodes<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;

    std::vector<typename MeshType::MTNode> nodes;
    double& dmax;

    GenericAlgorithm_snap_nodes(RefMeshSmoother *rms, PerceptMesh *eMesh, double& dm);

    void run()
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index);


  };


  template<>
  GenericAlgorithm_snap_nodes<STKMesh>::
  GenericAlgorithm_snap_nodes(RefMeshSmoother *rms, PerceptMesh *eMesh, double& dm) : m_rms(rms), m_eMesh(eMesh), dmax(dm)
  {
    coord_field = m_eMesh->get_coordinates_field();
    coord_field_current   = coord_field;

    on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = m_eMesh->get_spatial_dim();

    dmax=0.0;

    // node loop - build list...
    {
      nodes.resize(0);
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // update local and globally shared
          if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  nodes.push_back(node);
                }
            }
        }
    }
  }

  template<typename MeshType>
  void GenericAlgorithm_snap_nodes<MeshType>::
  operator()(int64_t& index)
  {
    typename MeshType::MTNode node = nodes[index];
    std::pair<bool,int> fixed = m_rms->get_fixed_flag(node);
    if (fixed.first)
      {
        return;
      }

    double coord_current[3];
    get_field<MeshType>(coord_current, spatialDim, m_eMesh, coord_field_current, node);

    double coord_project[3] = {0,0,0};
    double coord_old[3] = {0,0,0};
    for (int i=0; i < spatialDim; i++)
      {
        coord_project[i] = coord_current[i];
        coord_old[i] = coord_current[i];
      }

    if (fixed.second == MS_SURFACE)
      {
        bool keep_node_unchanged = false;
        //snap_to(node, coord_project, keep_node_unchanged);
        m_rms->snap_to(node, &coord_current[0], keep_node_unchanged);
      }

    double dm = 0.0;
    for (int i=0; i < spatialDim; i++)
      {
        dm += (coord_old[i] - coord_project[i])*(coord_old[i] - coord_project[i]);
      }
    dm = std::sqrt(dm);
    dmax = std::max(dmax, dm);
  }

  template<>
  void ReferenceMeshSmootherConjugateGradientImpl< STKMesh >::
  snap_nodes()
  {
    static int anim_step = 0;
    bool save_anim = m_eMesh->getProperty("ReferenceMeshSmootherConjugateGradientImpl.snap_nodes.save_anim") == "true";
    if (save_anim)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);
        std::string oname = "snap.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        m_eMesh->save_as(oname);
        ++anim_step;
      }

    std::string prevSetting;
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    if (m_meshGeometry && m_meshGeometry->geomKernel)
      {
        prevSetting = m_meshGeometry->geomKernel->get_property("GKGP:use_unprojected_coordinates");
        m_meshGeometry->geomKernel->set_property("GKGP:use_unprojected_coordinates", "false");
      }
#endif

    double dmax=0.0;

    GenericAlgorithm_snap_nodes<STKMesh> ga(this, m_eMesh, dmax);
    ga.run();

    if (save_anim)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);
        std::string oname = "snap.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        m_eMesh->save_as(oname);
        ++anim_step;
      }

    // FIXME - add coordinate comm
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    m_meshGeometry->geomKernel->set_property("GKGP:use_unprojected_coordinates", prevSetting);
#endif
  }

template<typename MeshType>
Double ReferenceMeshSmootherConjugateGradientImpl<MeshType>::total_metric(
        Double alpha, double /*multiplicative_edge_scaling*/, bool& valid,
        size_t* num_invalid) {
    Double mtot = Double(0.0);
    size_t n_invalid = 0;

    m_coord_updater.alpha = alpha;
    m_coord_updater.run();

    valid = true;


    if (m_eMesh->get_bulk_data()) {//if you have a stk mesh
        //for each block      //madbrew: how is the idea of an unstructured block represented on a stk mesh?
        m_metric_computinator.updateState(this->m_metric, m_eMesh, valid,
                num_invalid, mtot, n_invalid);
        m_metric_computinator.run(/*iBlock*/0);
        valid = m_metric_computinator.valid;
        mtot += m_metric_computinator.mtot;
        n_invalid +=m_metric_computinator.n_invalid;
    } //stk mesh code

    // reset coordinates
//    m_eMesh->copy_field(m_coord_updater->coord_field_current,
//            m_coord_updater->coord_field_lagged);
    // reset coordinates
    m_eMesh->copy_field("coordinates", "coordinates_lagged");
    stk::all_reduce(m_eMesh->parallel(), stk::ReduceSum<1>(&mtot));
    stk::all_reduce(m_eMesh->parallel(), stk::ReduceMin<1>(&valid));

    if (num_invalid) {
        *num_invalid = n_invalid;
        stk::all_reduce(m_eMesh->parallel(), stk::ReduceSum<1>(num_invalid));
    }

    return mtot;
}	      //total_metric

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  update_node_positions
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  update_node_positions( Double alpha)
  {
      m_coord_updater.alpha=alpha;
      m_coord_updater.run(true);

    stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & m_dmax ) );
    stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & m_dmax_relative ) );

    {
      std::vector< const typename MeshType::MTField *> fields;
      fields.push_back(m_coord_updater.coord_field);
      MTcommFields<MeshType>(fields, m_eMesh);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  line_search
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_line_search
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using This = GenericAlgorithm_line_search<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_g_field;
    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_r_field;

    bool restarted;
    double mfac_mult;
    Double alpha;
    bool extra_print;
    GenericAlgorithm_line_search(RefMeshSmoother *rms, PerceptMesh *eMesh, double mf_mult);
    void run();
  };

  template<>
  GenericAlgorithm_line_search<STKMesh>::
  GenericAlgorithm_line_search(RefMeshSmoother *rms, PerceptMesh *eMesh, double mf_mult) : m_rms(rms), m_eMesh(eMesh), mfac_mult(mf_mult)
  {
    restarted = false;
    extra_print = false;
    if (m_eMesh->getProperty("ReferenceMeshSmootherConjugateGradientImpl::line_search.extra_print") == "true")
      extra_print = true;
    cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
  }

  template<typename MeshType>
  void GenericAlgorithm_line_search<MeshType>::
  run()
  {
    const Double alpha_fac = 10.0;
    const Double tau = 0.5;
    const Double c0 = 1.e-1;
    const Double min_alpha_factor = 1.e-12;

    PerceptMesh *eMesh = m_eMesh;
    m_rms->m_alpha_0 = m_rms->get_alpha_0();

    alpha = alpha_fac*m_rms->m_alpha_0;

    bool total_valid = false;
    size_t n_invalid = 0;
    size_t* n_invalid_p = &n_invalid;
    const Double metric_0 = m_rms->total_metric( 0.0, 1.0, total_valid, n_invalid_p);
    Double metric = 0.0;

    Double sDotGrad = eMesh->nodal_field_dot("cg_s", "cg_g");
    if (sDotGrad >= 0.0)
      {
        eMesh->copy_field("cg_s", "cg_r");
        m_rms->m_alpha_0 = m_rms->get_alpha_0();
        alpha = alpha_fac*m_rms->m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot("cg_s", "cg_g");
        restarted = true;
      }
    VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad");

    const Double armijo_offset_factor = c0*sDotGrad;
    bool converged = false;
    total_valid = false;
    int liter = 0, niter = 1000;
    while (!converged && liter < niter)
      {
        metric = m_rms->total_metric(alpha, 1.0, total_valid, n_invalid_p);

        const Double mfac = alpha*armijo_offset_factor * mfac_mult;
        converged = (metric < metric_0 + mfac);
        if (m_rms->m_untangled) converged = converged && total_valid;
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_rms->m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        restarted = true;
        eMesh->copy_field("cg_s", "cg_r");
        m_rms->m_alpha_0 = m_rms->get_alpha_0();
        alpha = alpha_fac*m_rms->m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot("cg_s", "cg_g");
        RMSCG_PRINT_1("not converged, trying restart, sDotGrad new= " << sDotGrad << " m_stage= " << m_rms->m_stage << " m_iter= " << m_rms->m_iter);
        VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad 2nd time");
      }

    liter = 0;
    while (!converged && liter < niter)
      {
        metric = m_rms->total_metric(alpha, 1.0, total_valid);

        const Double mfac = alpha*armijo_offset_factor;
        converged = (metric < metric_0 + mfac);
        if (m_rms->m_untangled)
          {
            converged = converged && total_valid;
            if (total_valid && m_rms->m_dmax_relative < m_rms->gradNorm)
              converged = true;
          }

        RMSCG_PRINT_1(  "alpha 2nd time= " << alpha << " alpha_0= " << m_rms->m_alpha_0 << " sDotGrad= " << sDotGrad << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac)
                  << " m_untangled = " << m_rms->m_untangled << " m_stage= " << m_rms->m_stage
                  << " total_valid= " << total_valid );
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_rms->m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        RMSCG_PRINT_1( "WARNING: can't reduce metric 2nd time = " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor << " sDotGrad = " << sDotGrad << " alpha_0= " << m_rms->m_alpha_0 << " alpha= " << alpha);
        throw std::runtime_error("can't reduce metric");
      }
    else
      {
        Double a1 = alpha/2.;
        Double a2 = alpha;
        Double f0 = metric_0, f1 = m_rms->total_metric( a1, 1.0, total_valid), f2 = m_rms->total_metric( a2, 1.0, total_valid);
        Double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
        Double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
        if (std::fabs(den) > 1.e-10)
          {
            Double alpha_quadratic = num/den;
            if (alpha_quadratic > 1.e-10 && alpha_quadratic < 2*alpha)
              {
                Double fm=m_rms->total_metric( alpha_quadratic, 1.0, total_valid);
                if (fm < f2 && (m_rms->m_stage==0 || total_valid))
                  {
                    alpha = alpha_quadratic;
                  }
                if (fm < f2 && (m_rms->m_stage!=0 && !total_valid))
                  {
                    RMSCG_PRINT_1( "WARNING: !total_valid alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 );
                  }
              }
          }
      }

  }//run

  template<typename MeshType>
  Double
  ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  line_search(bool& restarted, double mfac_mult)
  {
    GenericAlgorithm_line_search<MeshType> gal(this, m_eMesh, mfac_mult);
    gal.run();
    restarted = gal.restarted;
    return gal.alpha;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_surface_normals
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(STK_PERCEPT_HAS_GEOMETRY)
  template<typename MeshType>
  struct GenericAlgorithm_get_surface_normals
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using This = GenericAlgorithm_get_surface_normals<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_normal_field;
    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    std::vector<double> norm;
    std::vector<typename MeshType::MTNode> nodes;

    GenericAlgorithm_get_surface_normals(RefMeshSmoother *rms, PerceptMesh *eMesh);

    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];
      double cg_normal[3];
      get_field<MeshType>(cg_normal, m_eMesh->get_spatial_dim(), m_eMesh, cg_normal_field, node);

      std::pair<bool,int> fixed = m_rms->get_fixed_flag(node);
      if (!fixed.first && (fixed.second == MS_SURFACE || fixed.second == MS_ON_BOUNDARY))
        {
          m_rms->m_meshGeometry->normal_at(m_eMesh, node, norm);
          for (int ii=0; ii < m_eMesh->get_spatial_dim(); ++ii)
            {
              cg_normal[ii] = norm[ii];
            }
          set_field<MeshType>(cg_normal, m_eMesh->get_spatial_dim(), m_eMesh, cg_normal_field, node);
        }
    }

  };

  template<>
  GenericAlgorithm_get_surface_normals<STKMesh>::
  GenericAlgorithm_get_surface_normals(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    norm.resize(3,0.0);
    cg_normal_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");
    VERIFY_OP_ON(cg_normal_field, !=, 0, "must have cg_normal_field");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    nodes.resize(0);
    const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                nodes.push_back(node);
              }
          }
      }
  }

#endif

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_surface_normals(PerceptMesh * eMesh)
  {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    GenericAlgorithm_get_surface_normals< MeshType > ga_gsn(this, eMesh);
    ga_gsn.run();
#else
    VERIFY_MSG("you must configure Percept with geometry (STK_PERCEPT_HAS_GEOMETRY=1) to use get_surface_normals");
#endif
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_edge_lengths
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<typename MeshType>
  struct GenericAlgorithm_get_edge_lengths
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;
    typename MeshType::MTField *cg_edge_length_field;
    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    std::vector<typename MeshType::MTNode> nodes;

    using This = GenericAlgorithm_get_edge_lengths<MeshType>;

    GenericAlgorithm_get_edge_lengths(RefMeshSmoother *rms, PerceptMesh *eMesh);


    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];
      double cg_edge_length[1];
      get_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);

      //if (on_locally_owned_part(node) || on_globally_shared_part(node))
      {
        Double edge_length_ave = m_rms->nodal_edge_length_ave(node);
        cg_edge_length[0] = edge_length_ave;
        set_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
      }

    }
  };

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_edge_lengths(PerceptMesh * eMesh)
  {
    if (eMesh->get_bulk_data()) {
        GenericAlgorithm_get_edge_lengths<MeshType> gae(this, eMesh);
        gae.run();
    }
  }

  template<>
  GenericAlgorithm_get_edge_lengths<STKMesh>::
  GenericAlgorithm_get_edge_lengths(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    nodes.resize(0);

    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");
    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                nodes.push_back(node);
              }
          }
      }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_alpha_0
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  template<typename MeshType>
  struct GenericAlgorithm_get_alpha_0
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_edge_length_field;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;

    Double alpha;
    bool alpha_set;

    std::vector<typename MeshType::MTNode> nodes;
    using This = GenericAlgorithm_get_alpha_0<MeshType>;

    GenericAlgorithm_get_alpha_0(RefMeshSmoother *rms, PerceptMesh *eMesh);

    void run()
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }

      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & alpha_set ) );
      if (!alpha_set)
        alpha = 1.0;

      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMin<1>( & alpha ) );
      VERIFY_OP_ON(alpha, > , 0.0, "bad alpha");

    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];

      //VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");

      double cg_edge_length[1];
      get_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
      Double edge_length_ave = cg_edge_length[0];

      bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
      VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
      bool fixed = m_rms->get_fixed_flag(node).first;
      if (fixed || isGhostNode)
        return;

      double cg_s[3];
      get_field<MeshType>(cg_s, spatialDim, m_eMesh, cg_s_field, node);
      Double sn = 0.0;
      for (int idim=0; idim < spatialDim; idim++)
        {
          sn += cg_s[idim]*cg_s[idim];
        }
      sn = std::sqrt(sn);
      if (sn > 0.0)
        {
          Double alpha_new = edge_length_ave / sn;
          if (!alpha_set)
            {
              alpha_set = true;
              alpha = alpha_new;
            }
          else if (alpha_new < alpha)
            {
              alpha = alpha_new;
            }
        }
    }
  };

  /// gets a global scale factor so that local gradient*scale is approximately the size of the local mesh edges
  /// also uses the reference mesh to compute a local scaled gradient norm for convergence checks
  template<typename MeshType>
  Double
  ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_alpha_0()
  {
      Double alpha_tot=0.0;
      if (m_eMesh->get_bulk_data()) {
        GenericAlgorithm_get_alpha_0<MeshType> ga0(this, m_eMesh);
        ga0.run();
        alpha_tot += ga0.alpha;
    }

    return alpha_tot;
  }

  template<>
  GenericAlgorithm_get_alpha_0<STKMesh>::
  GenericAlgorithm_get_alpha_0(RefMeshSmoother * rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = eMesh->get_spatial_dim();

    alpha = std::numeric_limits<double>::max();
    alpha_set = false;

    nodes.resize(0);

    {
      const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // FIXME
          if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh) && (on_locally_owned_part(**k) || on_globally_shared_part(**k)))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");
                  nodes.push_back(node);
                }
            }
        }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  check_convergence
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<typename MeshType>
  bool ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  check_convergence()
  {
    Double grad_check = gradNorm;
    bool retval=false;
    if (m_stage == 0 && (m_dnew == 0.0 || m_total_metric == 0.0))
      {
        retval = true; // for untangle
      }
    else if (m_stage == 0 && m_num_invalid == 0 && (m_dmax_relative < grad_check))
      {
        retval = true;
      }
    else if (m_num_invalid == 0 && (m_iter > 10 && (m_dmax_relative < grad_check && (m_dnew < grad_check*grad_check*m_d0 || m_grad_norm_scaled < grad_check))))
      {
        retval = true;
      }
    return retval;
  }

  template<typename MeshType>
  double ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  run_one_iteration()
  {
	stk::diag::Timer cumulative_timer(m_eMesh->getProperty("in_filename"), rootTimerStructured());
	stk::diag::Timer runOneIter("RMSCGI::run_one_iteration()",cumulative_timer);
	stk::diag::TimeBlock my_time_block(runOneIter);

    PerceptMesh *m_eMesh2 = Base::m_eMesh;

    bool total_valid=false;

    if (Base::m_iter == 0)
      {
        get_edge_lengths(m_eMesh2);
        m_eMesh2->nodal_field_set_value("cg_r", 0.0);
        m_eMesh2->nodal_field_set_value("cg_d", 0.0);
        m_eMesh2->nodal_field_set_value("cg_s", 0.0);
      }

    /// f'(x)
    get_gradient();

    /// r = -g
    m_eMesh2->nodal_field_axpby(-1.0, "cg_g", 0.0, "cg_r");
    Base::m_dold = m_eMesh2->nodal_field_dot("cg_d", "cg_d");
    Base::m_dmid = m_eMesh2->nodal_field_dot("cg_r", "cg_d");
    Base::m_dnew = m_eMesh2->nodal_field_dot("cg_r", "cg_r");
    if (Base::m_iter == 0)
      {
        Base::m_d0 = Base::m_dnew;
      }

    Double metric_check = total_metric(0.0, 1.0, total_valid);
    Base::m_total_metric = metric_check;

    if (check_convergence() || metric_check == 0.0)
      {
        return total_metric(0.0, 1.0, total_valid);
      }

    Double cg_beta = 0.0;
    if (Base::m_dold == 0.0)
      cg_beta = 0.0;
    else if (Base::m_iter > 0)
      cg_beta = (Base::m_dnew - Base::m_dmid) / Base::m_dold;

    size_t N = Base::m_num_nodes;
    if (Base::m_iter % N == 0 || cg_beta <= 0.0)
      {
        /// s = r
        m_eMesh2->copy_field("cg_s", "cg_r");
      }
    else
      {
        /// s = r + beta * s
        m_eMesh2->nodal_field_axpby(1.0, "cg_r", cg_beta, "cg_s");
      }

    m_eMesh2->copy_field("cg_d", "cg_r");
    bool restarted = false;
    Double alpha = line_search(restarted);
    Double snorm = m_eMesh2->nodal_field_dot("cg_s", "cg_s");
    Base::m_grad_norm_scaled = Base::m_alpha_0*std::sqrt(snorm)/Double(Base::m_num_nodes);

    /// x = x + alpha*s
    Base::m_alpha = alpha;
    update_node_positions(alpha);
    // check if re-snapped geometry is acceptable
    if (m_eMesh2->get_smooth_surfaces())
      {
        snap_nodes();
        if (Base::m_stage != 0)
          {
            bool total_valid_0=true;
            total_metric( 0.0, 1.0, total_valid_0);
//            std::cout<<"total metric comp"<<std::endl;
            VERIFY_OP_ON(total_valid_0, ==, true, "bad mesh after snap_node_positions...");
          }
      }

    Double tm = total_metric(0.0,1.0, total_valid);
    return tm;
  }//run_one_iteration()


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_gradient
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
  template<class MeshType>
  struct GenericAlgorithmBase_get_gradient {

    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;
    typename MeshType::MTField *cg_g_field;
    typename MeshType::MTField *cg_r_field;
    typename MeshType::MTField *cg_edge_length_field;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;

    int spatialDim;

    GenericAlgorithmBase_get_gradient(RefMeshSmoother *rms, PerceptMesh *eMesh);

  };

  template<>
  GenericAlgorithmBase_get_gradient<STKMesh>::
  GenericAlgorithmBase_get_gradient(RefMeshSmoother *rms, PerceptMesh *eMesh) :  m_rms(rms), m_eMesh(eMesh), spatialDim(eMesh->get_spatial_dim())
  {
    rms->m_scale = 1.e-10;

    coord_field = eMesh->get_coordinates_field();
    coord_field_current   = coord_field;
    cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

  }

  template<class MeshType>
  struct GenericAlgorithm_get_gradient_1 : public GenericAlgorithmBase_get_gradient<MeshType> {

    using Base = GenericAlgorithmBase_get_gradient<MeshType> ;
    using This = GenericAlgorithm_get_gradient_1<MeshType>;

    using Base::m_eMesh;
    using Base::coord_field;
    using Base::coord_field_current;
    using Base::cg_g_field;
    using Base::cg_r_field;
    using Base::cg_edge_length_field;

    using Base::on_locally_owned_part;
    using Base::on_globally_shared_part;
    using Base::spatialDim;

    using RefMeshSmoother = ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    GenericAlgorithm_get_gradient_1(RefMeshSmoother *rms, PerceptMesh *eMesh) : Base(rms, eMesh)
    {
    }

    std::vector<typename MeshType::MTElement> elements;
    std::vector<const typename MeshType::MTCellTopology *> topos;

    void MTSetTopo(int64_t index) const;

    void run() const
    {

      for (int64_t index = 0; index < (int64_t)elements.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {

      const double macheps = std::numeric_limits<double>::epsilon();
      const double sqrt_eps = std::sqrt(macheps);

      typename MeshType::MTElement element = elements[index];
      MTSetTopo(index);

      unsigned num_node = get_num_nodes<MeshType>(m_eMesh, element);
      std::vector<typename MeshType::MTNode> nodes_buffer;
      const typename MeshType::MTNode *elem_nodes = get_nodes<MeshType>(m_eMesh, element, &nodes_buffer);

      Double edge_length_ave = 0;

      const bool use_analytic_grad = true;
      double analytic_grad[8][4];
      const bool test_analytic_grad = false;
      VERIFY_OP_ON(Base::m_rms->m_metric->has_gradient(), ==, true, "bad has_gradient");
      if ((test_analytic_grad || use_analytic_grad) && Base::m_rms->m_stage >= 0 && Base::m_rms->m_metric->has_gradient())
        {
          bool gmvalid = true;
          Double gm = Base::m_rms->m_metric->grad_metric(element, gmvalid, analytic_grad);
          (void)gm;
          if ((gmvalid || Base::m_rms->m_stage == 0) && !test_analytic_grad)
            {

              for (unsigned inode=0; inode < num_node; inode++)
                {
                  typename MeshType::MTNode node = elem_nodes[ inode ];

                  bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
                  bool node_locally_owned = MTnode_locally_owned<MeshType>(m_eMesh, node);
                  bool fixed = Base::m_rms->get_fixed_flag(node).first;

                  if (fixed || isGhostNode)
                    continue;

                  VERIFY_OP_ON(Base::spatialDim, ==, spatialDim, "bad spatialDim");
                  VERIFY_OP_ON(Base::spatialDim, >=, 2, "bad spatialDim");
                  double cg_g[3];
                  get_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);

                  for (int jdim=0; jdim < spatialDim; jdim++)
                    {
                      if (node_locally_owned)
                        cg_g[jdim] += analytic_grad[inode][jdim];
                      else
                        cg_g[jdim] = 0.0;
                    }
                  set_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                }
            }
          if (!test_analytic_grad)
            return;
        }

      if (1)
        {
          // finite-different grad
          for (unsigned inode=0; inode < num_node; inode++)
            {
              typename MeshType::MTNode node = elem_nodes[ inode ];

              bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
              bool node_locally_owned = MTnode_locally_owned<MeshType>(m_eMesh, node);
              bool fixed = Base::m_rms->get_fixed_flag(node).first;
              if (fixed || isGhostNode)
                continue;

              double cg_edge_length[1];
              get_field<MeshType>(cg_edge_length, 1, Base::m_rms->m_eMesh, cg_edge_length_field, node);

              edge_length_ave = cg_edge_length[0];

              Base::m_rms->m_metric->set_node(node);
              double coord_current[3];
              get_field<MeshType>(coord_current, spatialDim, Base::m_rms->m_eMesh, coord_field_current, node);
              double cg_g[3];
              get_field<MeshType>(cg_g, spatialDim, Base::m_rms->m_eMesh, cg_g_field, node);

              Double eps1 = sqrt_eps*edge_length_ave;

              double gsav[3]={0,0,0};
              for (int idim=0; idim < spatialDim; idim++)
                {
                  Double cc = coord_current[idim];
                  coord_current[idim] += eps1;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
                  bool pvalid=false, mvalid=false;
                  Double mp = Base::m_rms->m_metric->metric(element, pvalid);
                  const bool second_order = true;
                  if (second_order)
                    coord_current[idim] -= 2.0*eps1;
                  else
                    coord_current[idim] = cc;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
                  Double mm = Base::m_rms->m_metric->metric(element, mvalid);
                  coord_current[idim] = cc;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
                  Double dd = 0.0;
                  if (second_order)
                    {
                      dd = (mp - mm)/(2*eps1);
                    }
                  else
                    {
                      dd = (mp - mm)/(eps1);
                    }
                  gsav[idim] = dd;

                  if (node_locally_owned)
                    {
                      cg_g[idim] += dd;
                    }
                  else
                    {
                      cg_g[idim] = 0.0;
                    }
                  set_field<MeshType>(cg_g, spatialDim, idim, Base::m_rms->m_eMesh, cg_g_field, node);
                }
              if (test_analytic_grad)
                {
                  double fd = std::max(std::fabs( gsav[0]),std::fabs( gsav[1]));
                  if (spatialDim==3) fd = std::max(fd, std::fabs( gsav[2]));
                  double ag = std::max(std::fabs( analytic_grad[inode][0]),std::fabs( analytic_grad[inode][1] ));
                  if (spatialDim==3) ag = std::max(ag, std::fabs( analytic_grad[inode][2] ));
                  double diff = std::fabs(ag-fd);
                  double comp = (ag+fd)/2.0;
                  if (comp > 1.e-6 && diff > 1.e-8 && diff > 1.e-3*comp)
                    {
                      std::cout << "analytic_grad= " << analytic_grad[inode][0] << " " << analytic_grad[inode][1] << " "
                                << " fd_grad= " << gsav[0] << " " << gsav[1]
                                << " diff = " << analytic_grad[inode][0]-gsav[0] << " " << analytic_grad[inode][1] -gsav[1]
                                << " ag= " << ag << " fd= " << fd
                                << std::endl;
                    }

                }
            }
        }
    }
  };


  template<>
  void GenericAlgorithm_get_gradient_1<STKMesh>::
  MTSetTopo(int64_t index) const
  {
    Base::m_rms->m_metric->m_topology_data = topos[index];
  }

  template<>
  GenericAlgorithm_get_gradient_1<STKMesh>::
  GenericAlgorithm_get_gradient_1(ReferenceMeshSmootherConjugateGradientImpl<STKMesh> *rms, PerceptMesh *eMesh)
    : GenericAlgorithmBase_get_gradient<STKMesh>(rms, eMesh)
  {
    eMesh->nodal_field_set_value("cg_g", 0.0);

    // get elements
    if (1)
      {
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                Base::m_rms->m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];

                    if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                      continue;
                    elements.push_back(element);
                    topos.push_back(m_eMesh->get_cell_topology(bucket));
                  }
              }
          }
      }
  }

  template<typename MeshType>
  struct GenericAlgorithm_get_gradient_2 : public GenericAlgorithmBase_get_gradient<MeshType>
  {
    std::vector<typename MeshType::MTNode> nodes;

    using Base = GenericAlgorithmBase_get_gradient<MeshType>;
    using This = GenericAlgorithm_get_gradient_2<MeshType>;

    using Base::m_eMesh;
    using Base::spatialDim;
    using Base::cg_g_field;
    using Base::cg_r_field;

    GenericAlgorithm_get_gradient_2(ReferenceMeshSmootherConjugateGradientImpl<MeshType> *rms, PerceptMesh *eMesh);

    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      // project deltas to surface
      if (m_eMesh->get_smooth_surfaces())
        {
          typename MeshType::MTNode node = nodes[index];

          std::pair<bool,int> fixed = Base::m_rms->get_fixed_flag(node);
          if (!fixed.first)
            {
              if (fixed.second == MS_SURFACE)
                {
                  double cg_g[3];
                  get_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                  Base::m_rms->project_delta_to_tangent_plane(node, cg_g);
                  set_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                }
            }
        }
    }
  };

  template<>
  GenericAlgorithm_get_gradient_2<STKMesh>::
  GenericAlgorithm_get_gradient_2(ReferenceMeshSmootherConjugateGradientImpl<STKMesh> *rms, PerceptMesh *eMesh) :  GenericAlgorithmBase_get_gradient<STKMesh>(rms, eMesh)
  {
    nodes.clear();

    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                nodes.push_back(bucket[i_node]);
              }
          }
      }
  }

//  /// fills cg_g_field with f'(x)
  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_gradient()
  {
      if (m_eMesh->get_bulk_data()) {//if you have a stk mesh
          GenericAlgorithm_get_gradient_1<MeshType> ga_1(this, m_eMesh);
          ga_1.run();

          //Double gnorm = m_eMesh->nodal_field_dot(ga_1.cg_g_field, ga_1.cg_g_field);
          //std::cout << "gnorm= " << gnorm << std::endl;


          std::vector<const typename MeshType::MTField *> fields_0(1, ga_1.cg_g_field);
          MTsum_fields<MeshType>(fields_0, m_eMesh);//madbrew: for now doesn't do anything for sgrids

          m_eMesh->copy_field("cg_r", "cg_g");

          GenericAlgorithm_get_gradient_2<MeshType> ga_2(this, m_eMesh);
          ga_2.run(); //madbrew: for now this one doesn't do anything for sgrids because  Base::m_rms->project_delta_to_tangent_plane(node, cg_g)    has no sgrdi impl

          {
            std::vector<  const  typename MeshType::MTField *> fields;
            fields.push_back(ga_1.cg_g_field);
            fields.push_back(ga_1.cg_r_field);
            MTcommFields<MeshType>(fields, m_eMesh); //madbrew: for now doesn't do anything for sgrids
          }
        }//stk mesh code
      } //get_gradient()
}//percept

#endif
