#ifndef stk_adapt_TestLocalRefinerTri_N_3_EdgeBasedAnisotropic_hpp
#define stk_adapt_TestLocalRefinerTri_N_3_EdgeBasedAnisotropic_hpp

#include <stk_adapt/IEdgeAdapter.hpp>

#include <stk_mesh/base/GetBuckets.hpp>

namespace stk {
  namespace adapt {

    void exact_hessian(const int id, const double *xyz, double *hess, const int spatial_dim)
    {
      switch (id) 
	{
	case 0: 
	  {
	    // simple poly with max Hess at x=1
	    const double x = xyz[0];
	    hess[0] = 200*x*x;
	    hess[1] = 0;
	    hess[2] = 0;
	    hess[3] = 0;
	  }
	  break;
	case 1: 
	  {
	    // simple poly with max Hess at y=1
	    const double y = xyz[1];
	    hess[0] = 0;
	    hess[1] = 0;
	    hess[2] = 0;
	    hess[3] = 200*y*y;
	  }
	  break;
	case 2: 
	  {
	    // fully 2d case
	    const double x = xyz[0];
	    const double y = xyz[1];
	    hess[0] = 200*x*x;
	    hess[1] = 0;
	    hess[2] = 0;
	    hess[3] = 200*y*y;
	  }
	  break;
	case 3: 
	  {
	    // full tensor 2d case
	    // f = 100*(y-x+1)*(y-x-1) = 100*y*y+100*y*y+(-200)*x*y + linear terms
	    hess[0] = 200;
	    hess[1] = -200;
	    hess[2] = -200;
	    hess[3] = 200;
	  }
	  break;
	default: 
	  {
	    std::ostringstream msg;
	    msg << " exact_hessian unknown id = " << id << "\n";
	    throw new std::runtime_error(msg.str());
	  }
	}
    }

    void interp_nodal_hessian(
      const int hess_id,
      PerceptMesh& eMesh, 
      VectorFieldType* nodal_hessian_field)
    {
      const int spatial_dim = eMesh.get_spatial_dim();

      stk::mesh::Selector selector = 
	eMesh.get_fem_meta_data()->locally_owned_part() | 
	eMesh.get_fem_meta_data()->globally_shared_part();

      std::vector<stk::mesh::Bucket*> buckets;
      stk::mesh::get_buckets( selector, eMesh.get_bulk_data()->buckets( eMesh.node_rank() ), buckets );
      
      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {

	stk::mesh::Bucket & bucket = **k ;
	
	const unsigned num_nodes_in_bucket = bucket.size();	

	for (unsigned i = 0; i < num_nodes_in_bucket; i++) {

	  stk::mesh::Entity& node = bucket[i];

	  const double *coords = stk::mesh::field_data( *eMesh.get_coordinates_field() , node);
	  double *hess = stk::mesh::field_data( *nodal_hessian_field , node);
	    
	  exact_hessian(hess_id, coords, hess, spatial_dim);
	}
      }
    }
    
    /**
     * A test implementation as a use case for EdgeBasedAnisotropic
     */
    class TestLocalRefinerTri_N_3_EdgeBasedAnisotropic : public IEdgeAdapter
    {
    public:
      TestLocalRefinerTri_N_3_EdgeBasedAnisotropic(
        percept::PerceptMesh& eMesh, 
	UniformRefinerPatternBase & bp, 
	VectorFieldType * nodal_hessian_field,
	stk::mesh::FieldBase *proc_rank_field=0);

    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int mark(
        const stk::mesh::Entity& element, 
	unsigned which_edge, 
	stk::mesh::Entity & node0, 
	stk::mesh::Entity & node1,
	double *coord0, 
	double *coord1, 
	std::vector<int>* existing_edge_marks) ;

      VectorFieldType * m_nodal_hessian_field;
      
      const double m_upper_bound, m_lower_bound;
      const double length_min, length_max;
    };

    TestLocalRefinerTri_N_3_EdgeBasedAnisotropic::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic(
      percept::PerceptMesh& eMesh, 
      UniformRefinerPatternBase &  bp, 
      VectorFieldType * nodal_hessian_field,
      stk::mesh::FieldBase *proc_rank_field) 
      : 
      IEdgeAdapter(eMesh, bp, proc_rank_field),
      m_nodal_hessian_field(nodal_hessian_field),
      m_upper_bound(sqrt(2.0)),
      m_lower_bound(1./sqrt(2.0)),
      length_min(0.01),
      // for coarsen only => use length_max longer than any edge length
      length_max(2.0*sqrt(2.0)) // diagonal of 2x2 square
    {}

    int 
    TestLocalRefinerTri_N_3_EdgeBasedAnisotropic::mark(
      const stk::mesh::Entity& element, 
      unsigned which_edge, 
      stk::mesh::Entity & node0, 
      stk::mesh::Entity & node1,
      double *coord0, double *coord1, 
      std::vector<int>* existing_edge_marks)
    {
      int mark=0;

      const int spatial_dim = m_eMesh.get_spatial_dim();

      std::vector<double> hessian(spatial_dim*spatial_dim, 0.0);

      // Hessian from midpoint interp of nodal Hessian field
      const double *hess0 = stk::mesh::field_data( *m_nodal_hessian_field , node0);
      const double *hess1 = stk::mesh::field_data( *m_nodal_hessian_field , node1);
      for (int d=0; d<spatial_dim*spatial_dim; d++) {
 	hessian[d] += 0.5*(hess0[d]+hess1[d]);
      }

      std::vector<double> edge_vec(spatial_dim, 0.0);
      double length = 0;
      for (int d=0; d<spatial_dim; d++) {
	edge_vec[d] = coord0[d] - coord1[d];
	length += edge_vec[d]*edge_vec[d];
      }
      length = sqrt(length);

      // Test 1: modified eigenvalues of (diag) Hessian
      /*
      const double local_error_tol = 0.001;

      // HACK assuming diag hessian for now
      std::vector<double> lambda(spatial_dim, 0.0);
      for (int d=0; d<spatial_dim; d++) {
	lambda[d] = std::min(std::max(avg_hessian[d]/local_error_tol, 
				      1./(length_max*length_max)), 
			     1./(length_min*length_min));
      }

      // calc metric, length at midpoint
      double local_metric2 = 0;
      // HACK assuming diag hessian
      for (int d=0; d<spatial_dim; d++) {
	local_metric2 += lambda[d] * edge_vec[d] * edge_vec[d];
      }

      const double metric = sqrt(local_metric2) * length_max;
      */

      // Test 2: scaled, modified metric
      /*
      double local_metric = 0;
      double local_length = 0;
      for (int d1=0; d1<spatial_dim; d1++) {
	for (int d2=0; d2<spatial_dim; d2++) {
	  local_metric += avg_hessian[d1+d2*spatial_dim]*edge_vec[d1]*edge_vec[d2];
	}
	local_length += edge_vec[d1]*edge_vec[d1];
      }
      local_metric = sqrt(local_metric);
      local_length = sqrt(local_length);
      */
      // Test 2.1: metric is ratio of current to optimal mesh size
      /*
      double inv_h_opt = local_metric / local_length;

      // rescale to maintain min/max mesh size
      // length_min < h_opt < length_max
      //  => 1/length_max < inv_h_opt < 1/length_min
      // TEST for refinement only replace length_max with local_length
      //inv_h_opt = std::min(std::max(inv_h_opt, 1./length_max), 1./length_min);
      inv_h_opt = std::min(std::max(inv_h_opt, 1./local_length), 1./length_min);
      
      const double metric = 1. / (local_length * inv_h_opt);
      */

      // Test 2.2: mark if local_metric > 0
      /*
      const double metric = (local_metric > 0) ? 2*m_Rup : 0.5*m_Rup;
      */
      
      double seminorm = 0;
      for (int d1=0; d1<spatial_dim; d1++) {
	for (int d2=0; d2<spatial_dim; d2++) {
	  seminorm += hessian[d1+d2*spatial_dim]*edge_vec[d1]*edge_vec[d2];
	}
      }

      if ( seminorm < 0 ) {
	throw std::logic_error("hessian matrix is non SPD!!");
      }
      seminorm = sqrt(seminorm);
      
      //   ideally the seminorm should = (h / h_opt)
      //   to enforce bounds on h_opt, we bound the seminorm
      //   (h / length_max) <= seminorm <= (h / length_min)
      double metric = std::min(std::max(seminorm, 
					length / length_max), 
			       length / length_min);
      
      // TEST simplest case: metric > Rup
      if (metric > m_upper_bound) {
	mark |= DO_REFINE;
      }
//       else if (metric < m_lower_bound) {
// 	mark |= DO_UNREFINE;
//       }
//       else {
// 	mark |= DO_NOTHING;
//       }

      return mark;
    }
  }
}

#endif
