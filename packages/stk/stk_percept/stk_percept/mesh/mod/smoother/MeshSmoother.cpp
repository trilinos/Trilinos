#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__)

#include <stk_percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother.hpp>
#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>
//#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother2.hpp>
//#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother3.hpp>
#include <stk_percept/mesh/mod/smoother/SmootherMetric.hpp>

#include "mpi.h"

#define DEBUG_PRINT 0

namespace stk {
  namespace percept {



    /// preferred for parallel
    int MeshSmoother::parallel_count_invalid_elements(PerceptMesh *eMesh)
    {
      SmootherMetricUntangle utm(eMesh);
      stk::mesh::FieldBase *coord_field_current   = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_original  = eMesh->get_field("coordinates_NM1");
      JacobianUtil jacA, jacW;

      double detA_min = std::numeric_limits<double>::max();
      double detW_min = std::numeric_limits<double>::max();
      double shapeA_max = 0.0;
      double shapeW_max = 0.0;
      const bool get_mesh_diagnostics = false;

      int num_invalid=0;
      // element loop
      {
        stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmoother::select_bucket(**k, eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                const CellTopologyData* topology_data = eMesh->get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    bool valid=true;
                    if (get_mesh_diagnostics)
                      {
                        double A_ = 0.0, W_ = 0.0; // current and reference detJ
                        jacA(A_, *eMesh, element, coord_field_current, topology_data);
                        jacW(W_, *eMesh, element, coord_field_original, topology_data);

                        for (int i=0; i < jacA.m_num_nodes; i++)
                          {
                            double detAi = jacA.m_detJ[i];
                            double detWi = jacW.m_detJ[i];
                            DenseMatrix<3,3>& W = jacW.m_J[i];
                            DenseMatrix<3,3>& A = jacA.m_J[i];
                            if (detAi <= 0.)
                              {
                                valid = false;
                              }
                            detA_min = std::min(detA_min, detAi);
                            detW_min = std::min(detW_min, detWi);
                            double frobAi = std::sqrt(my_sqr_Frobenius(A));
                            double frobWi = std::sqrt(my_sqr_Frobenius(W));
                            double shapeAi = std::abs(frobAi*frobAi*frobAi/(3*std::sqrt(3.)*detAi) - 1.0);
                            double shapeWi = std::abs(frobWi*frobWi*frobWi/(3*std::sqrt(3.)*detWi) - 1.0);
                            shapeA_max = std::max(shapeA_max, shapeAi);
                            shapeW_max = std::max(shapeW_max, shapeWi);
                          }
                      }
                    else
                      {
                        utm.metric(element, valid);
                      }
                    if (!valid)
                      ++num_invalid;
                  }
              }
          }
      }
      stk::all_reduce( MPI_COMM_WORLD, stk::ReduceSum<1>( &num_invalid ) );
      if (get_mesh_diagnostics)
        {
          stk::all_reduce( MPI_COMM_WORLD, stk::ReduceMin<1>( &detA_min ) );
          stk::all_reduce( MPI_COMM_WORLD, stk::ReduceMin<1>( &detW_min ) );
          stk::all_reduce( MPI_COMM_WORLD, stk::ReduceMax<1>( &shapeA_max ) );
          stk::all_reduce( MPI_COMM_WORLD, stk::ReduceMax<1>( &shapeW_max ) );
          if (eMesh->get_rank() == 0)
            {
              std::cout << "P[0] detA_min= " << detA_min << " detW_min= " << detW_min
                        << " shapeA_max= " << shapeA_max << " shapeW_max= " << shapeW_max << std::endl;
            }
        }
      return num_invalid;
    }


    int MeshSmoother::count_invalid_elements()
    {

      int num_invalid = 0;
      throw std::runtime_error("not implemented");
      return num_invalid;
    }

    int MeshSmoother::
    classify_node(stk::mesh::Entity node, size_t& curveOrSurfaceEvaluator) const
    {
      int dof =0;
      if (m_meshGeometry)
        dof = m_meshGeometry->classify_node(node, curveOrSurfaceEvaluator);
      else
        dof = m_eMesh->get_spatial_dim();
      return dof;
    }

    bool MeshSmoother:: select_bucket(stk::mesh::Bucket& bucket, PerceptMesh *eMesh)
    {
      const CellTopologyData * cell_topo_data = eMesh->get_cell_topology(bucket);
      shards::CellTopology cell_topo(cell_topo_data);

      if (cell_topo.getKey() == shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()->key
          || cell_topo.getKey() == shards::getCellTopologyData<shards::ShellTriangle<3> >()->key
          || cell_topo.getKey() == shards::getCellTopologyData<shards::ShellLine<2> >()->key
          || cell_topo.getKey() == shards::getCellTopologyData<shards::Beam<2> >()->key)
        {
          return false;
        }
      else
        {
          return true;
        }
    }

    std::pair<bool,int> MeshSmoother::get_fixed_flag(stk::mesh::Entity node_ptr)
    {
      int dof = -1;
      std::pair<bool,int> ret(true,MS_VERTEX);
      //if the owner is something other than the top-level owner, the node
      // is on the boundary; otherwise, it isn't.
      bool& fixed = ret.first;
      int& type = ret.second;
      if (m_boundarySelector)
        {
          if ((*m_boundarySelector)(node_ptr))
            {
              fixed=true;
              type=MS_ON_BOUNDARY;
            }
          else
            {
              fixed=false;
              type = MS_NOT_ON_BOUNDARY;
            }
        }
      else
        {
          if (m_meshGeometry)
            {
              size_t curveOrSurfaceEvaluator;
              dof = m_meshGeometry->classify_node(node_ptr, curveOrSurfaceEvaluator);
              //std::cout << "tmp srk classify node= " << node_ptr->identifier() << " dof= " << dof << std::endl;
              // vertex
              if (dof == 0)
                {
                  fixed=true;
                  type=MS_VERTEX;
                }
              // curve (for now we hold these fixed)
              else if (dof == 1)
                {
                  fixed=true;
                  type=MS_CURVE;
                  //fixed=false;   // FIXME
                }
              // surface - also fixed
              else if (dof == 2)
                {
                  //fixed=false;
                  fixed=true;
                  if (m_eMesh->get_smooth_surfaces())
                    {
                      fixed = false;
                      //std::cout << "tmp srk found surface node unfixed= " << node_ptr.identifier() << std::endl;
                    }
                  type=MS_SURFACE;
                  if (DEBUG_PRINT) std::cout << "tmp srk found surface node unfixed= " << node_ptr.identifier() << std::endl;
                }
              // interior/volume - free to move
              else
                {
                  fixed=false;
                  type=MS_VOLUME;
                }
            }
          else
            {
              fixed=false;
              type=MS_VOLUME;
            }
        }
      if (DEBUG_PRINT) std::cout << "tmp srk classify node= " << node_ptr.identifier() << " dof= " << dof << " fixed= " << fixed << " type= " << type << std::endl;

      return ret;
    }

    void MeshSmoother::run( bool always_smooth, int debug)
    {
#ifdef USE_CALLGRIND
      CALLGRIND_START_INSTRUMENTATION
        CALLGRIND_TOGGLE_COLLECT
#endif
      PerceptMesh *eMesh = m_eMesh;

      int num_invalid = parallel_count_invalid_elements(eMesh);
      if (!m_eMesh->get_rank())
        std::cout << "\ntmp srk MeshSmoother num_invalid before= " << num_invalid
                      << (num_invalid ? " WARNING: invalid elements exist before  smoothing" :
                          (!always_smooth ? "WARNING: no smoothing requested since always_smooth=false" : " "))
                      << std::endl;
      //if (num_invalid) throw std::runtime_error("MeshSmoother can't start from invalid mesh...");

      if (always_smooth)
        {
          //int  msq_debug             = debug; // 1,2,3 for more debug info

          //bool do_untangle_only = false;
          std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk innerIter= " << innerIter << " parallelIterations= " << parallelIterations << std::endl;
          this->run_algorithm();

          //if (!m_eMesh->get_rank())

          num_invalid = parallel_count_invalid_elements(eMesh);
          //if (!m_eMesh->get_rank())
          std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk MeshSmoother num_invalid after= " << num_invalid << " "
                    << (num_invalid ? " ERROR still have invalid elements after smoothing" :
                        " SUCCESS: smoothed and removed invalid elements ")
                    << std::endl;
          MPI_Barrier( MPI_COMM_WORLD );
          std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk after barrier" << std::endl;
        }


#ifdef USE_CALLGRIND
      CALLGRIND_TOGGLE_COLLECT
        CALLGRIND_STOP_INSTRUMENTATION
#endif
    }


    void MeshSmoother::
    project_delta_to_tangent_plane(stk::mesh::Entity node, double *delta)
    {
      if (!m_meshGeometry)
        {
          return;
        }

      std::vector<double> normal(3,0.0);
      m_meshGeometry->normal_at(m_eMesh, node, normal);
      double dot=0.0;
      for (int i = 0; i < m_eMesh->get_spatial_dim(); i++)
        {
          dot += delta[i]*normal[i];
        }
      for (int i = 0; i < m_eMesh->get_spatial_dim(); i++)
        {
          delta[i] -= dot*normal[i];
        }
    }

    //! Modifies "coordinate" so that it lies on the
    //! domain to which "node" is constrained.
    //! The node determines the domain.  The coordinate
    //! is the proposed new position on that domain.
    void MeshSmoother::
    snap_to(stk::mesh::Entity node_ptr,
            double *coordinate, bool reset) const 
    {

      if (!m_meshGeometry) return;
      stk::mesh::FieldBase* field = m_eMesh->get_coordinates_field();
      double *f_data = PerceptMesh::field_data(field, node_ptr);
      double f_data_save[3] = {f_data[0], f_data[1], 0};
      if (m_eMesh->get_spatial_dim() > 2) f_data_save[2] = f_data[2];

      f_data[0] = coordinate[0];
      f_data[1] = coordinate[1];
      if (m_eMesh->get_spatial_dim() > 2) 
        f_data[2] = coordinate[2];

      static std::vector<stk::mesh::Entity > nodes(1);
      nodes[0] = node_ptr;
      m_meshGeometry->snap_points_to_geometry(m_eMesh, nodes);
      coordinate[0] = f_data[0];
      coordinate[1] = f_data[1];
      if (m_eMesh->get_spatial_dim() > 2)
        coordinate[2] = f_data[2];

#if 0
      //if (node_ptr->identifier() == 584)
      if (0)
        {
          std::cout << "tmp snap_to: node= " << node_ptr->identifier() << " orig= " 
                    << f_data_save[0] << " "
                    << f_data_save[1] << " "
                    << f_data_save[2] << " "
                    << " new= " 
                    << f_data[0] << " "
                    << f_data[1] << " "
                    << f_data[2] << " diff1= " 
                    << (std::fabs(f_data[0] - f_data_save[0])+
                        std::fabs(f_data[1] - f_data_save[1])+
                        std::fabs(f_data[2] - f_data_save[2]))
                    << std::endl;
        }
#endif

      // reset the node
      if (reset)
        {
          f_data[0] = f_data_save[0];
          f_data[1] = f_data_save[1];
          if (m_eMesh->get_spatial_dim() > 2)
            f_data[2] = f_data_save[2];
        }

    }


  }
}


#endif
