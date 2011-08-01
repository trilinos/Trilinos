/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <boost/lexical_cast.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/Refiner.hpp>

#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/BeamFixture.hpp>
#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>
#include <stk_percept/fixtures/WedgeFixture.hpp>

#include <stk_adapt/EdgeBasedMarkerPredicate.hpp>
#include <stk_adapt/FieldBasedMarkerPredicate.hpp>
#include <stk_adapt/PredicateBasedMarker.hpp>
#include <stk_adapt/PredicateBasedEdgeMarker.hpp>

#include <stk_percept/function/ElementOp.hpp>

namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {
      const std::string input_files_loc="./input_files/";
      const std::string output_files_loc="./output_files/";

      static std::string post_fix[4] = {"np0", "np1", "np2", "np3"};

#define EXTRA_PRINT 0

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void normalize(double input_normal[3], double normal[3])
      {
        double sum = std::sqrt(input_normal[0]*input_normal[0]+
                               input_normal[1]*input_normal[1]+
                               input_normal[2]*input_normal[2]);
        normal[0] = input_normal[0] / sum;
        normal[1] = input_normal[1] / sum;
        normal[2] = input_normal[2] / sum;
      }

      static void normalize(double input_output_normal[3])
      {
        normalize(input_output_normal, input_output_normal);
      }

      static double distance(double c0[3], double c1[3])
      {
        return std::sqrt((c0[0]-c1[0])*(c0[0]-c1[0]) + (c0[1]-c1[1])*(c0[1]-c1[1]) + (c0[2]-c1[2])*(c0[2]-c1[2]) );
      }

      static void difference(double v01[3], double c0[3], double c1[3])
      {
        v01[0] = c0[0] - c1[0];
        v01[1] = c0[1] - c1[1];
        v01[2] = c0[2] - c1[2];
      }
      static double dot(double c0[3], double c1[3])
      {
        return c0[0]*c1[0] + c0[1]*c1[1] + c0[2]*c1[2];
      }

      static double plane_dot_product(double plane_point[3], double plane_normal[3], double point[3])
      {
        double normal[3]={0,0,0};
        normalize(plane_normal, normal);
        double dot = 0.0;
        for (int i = 0; i < 3; i++)
          {
            dot += (point[i] - plane_point[i])*normal[i];
          }
        return dot;
      }


#if 0
      static void project_point_to_plane(double plane_point[3], double plane_normal[3], double point[3], double projected_point[3])
      {
        double normal[3]={0,0,0};
        normalize(plane_normal, normal);
        double dot = plane_dot_product(plane_point, plane_normal, point);
        for (int i = 0; i < 3; i++)
          {
            projected_point[i] = plane_point[i] + (point[i] - dot*normal[i]);
          }
      }
#endif

      class SetRefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetRefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {
        }

        virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          double plane_point[3] = {2,0,0};
          double plane_normal[3] = {1, .5, -.5};

          const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          double *f_data = PerceptMesh::field_data_entity(field, element);
          VectorFieldType* coordField = m_eMesh.getCoordinatesField();
                
          bool found = false;
          for (unsigned inode=0; inode < num_node-1; inode++)
            {
              mesh::Entity & node_i = * elem_nodes[ inode ].entity();
              double *coord_data_i = PerceptMesh::field_data(coordField, node_i);

              for (unsigned jnode=inode+1; jnode < num_node; jnode++)
                {
                  mesh::Entity & node_j = * elem_nodes[ jnode ].entity();
                  double *coord_data_j = PerceptMesh::field_data(coordField, node_j);

                  double dot_0 = plane_dot_product(plane_point, plane_normal, coord_data_i);
                  double dot_1 = plane_dot_product(plane_point, plane_normal, coord_data_j);

                  // if edge crosses the plane...
                  if (dot_0*dot_1 < 0)
                    {
                      found=true;
                      break;
                    }
                }
            }
          if (found)
            f_data[0] = 1.0;
          else
            f_data[0] = 0.0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      class SetUnrefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetUnrefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {}
        virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)
        {
          const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          double *f_data = PerceptMesh::field_data_entity(field, element);
          VectorFieldType* coordField = m_eMesh.getCoordinatesField();
                
          bool found = true;
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();
              double *coord_data = PerceptMesh::field_data(coordField, node);

              if (coord_data[0] < 0.0 || coord_data[1] < 0.0) // || coord_data[2] > 1.1)
                {
                  found=false;
                  break;
                }
            }
          if (found)
            f_data[0] = -1.0;
          else
            f_data[0] = 0.0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_FieldBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.addField("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* unrefine_field = eMesh.addField("unrefine_field", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            SetRefineField set_ref_field(eMesh);
            eMesh.elementOpLoop(set_ref_field, refine_field);

            SetUnrefineField set_unref_field(eMesh);
            eMesh.elementOpLoop(set_unref_field, unrefine_field);
            
            eMesh.saveAs( output_files_loc+"local_tet_N_5_FieldBased_0_"+post_fix[p_size]+".e");

            stk::mesh::Selector univ_selector(eMesh.getFEM_meta_data()->universal_part());

            PredicateBasedMarker<ElementFieldBasedRefinePredicate, ElementFieldBasedUnrefinePredicate>
              breaker(ElementFieldBasedRefinePredicate(univ_selector, refine_field, 0.0),
                      ElementFieldBasedUnrefinePredicate(univ_selector, unrefine_field, 0.0),
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 3; ipass++)
              {
                eMesh.elementOpLoop(set_ref_field, refine_field);
                eMesh.elementOpLoop(set_unref_field, unrefine_field);

                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.getRank() << "] done... ipass= " << ipass << std::endl;
                eMesh.saveAs(output_files_loc+"local_tet_N_5_FieldBased_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.saveAs(output_files_loc+"local_tet_N_5_FieldBased_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tet_N_5_FieldBased_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.saveAs( output_files_loc+"local_tet_N_5_FieldBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      struct MyEdgeBasedRefinePredicate : public EdgeBasedMarkerPredicate {
      
        MyEdgeBasedRefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
          EdgeBasedMarkerPredicate(selector, field, tolerance) {}

        /// Return true for refine, false for ignore
        bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                        double *coord0, double *coord1, std::vector<int>& existing_edge_marks)
        {
          if (m_selector(element))
            {
              double plane_point[3] = {2,0,0};
              //double plane_normal[3] = {1, .5, -.5};
              double plane_normal[3] = {1, 0, 0};

              double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
              double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

              // if edge crosses the plane...
              if (dot_0*dot_1 < 0)
                {
                  return true;
                }
            }
          return false;
        }

        //double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        //return m_selector(entity) && fdata[0] > 0;
      };

      struct MyEdgeBasedUnrefinePredicate : public EdgeBasedMarkerPredicate {
      
        MyEdgeBasedUnrefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
          EdgeBasedMarkerPredicate(selector, field, tolerance) {}

        /// Return true for unrefine, false for ignore
        bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                        double *coord0, double *coord1)
        {
          if (m_selector(element))
            {
              // do something
            }
          return false;
        }

      };

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.addField("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* unrefine_field = eMesh.addField("unrefine_field", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);

                SetUnrefineField set_unref_field(eMesh);
                eMesh.elementOpLoop(set_unref_field, unrefine_field);
              }

            eMesh.saveAs( output_files_loc+"local_tet_N_5_EdgeBased_0_"+post_fix[p_size]+".e");

            stk::mesh::Selector univ_selector(eMesh.getFEM_meta_data()->universal_part());

            PredicateBasedEdgeMarker<MyEdgeBasedRefinePredicate, MyEdgeBasedUnrefinePredicate>
              breaker(MyEdgeBasedRefinePredicate(univ_selector, refine_field, 0.0),
                      MyEdgeBasedUnrefinePredicate(univ_selector, unrefine_field, 0.0),
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 3; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);
                //eMesh.elementOpLoop(set_unref_field, unrefine_field);

                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.getRank() << "] done... ipass= " << ipass << std::endl;
                eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.saveAs( output_files_loc+"local_tet_N_5_EdgeBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      static double shock_function(double x)
      {
        // normalize by width
        double width = 1./15.0;
        x /= width;
        return std::tanh(x);
      }

      struct ShockBasedRefinePredicate : public EdgeBasedMarkerPredicate {
      
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase * m_nodal_refine_field;

        ShockBasedRefinePredicate(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh, stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
          EdgeBasedMarkerPredicate(selector, field, tolerance), m_eMesh(eMesh),m_nodal_refine_field(nodal_refine_field) {}

        /// Return true for refine, false for ignore
        bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                        double *coord0, double *coord1, std::vector<int>& existing_edge_marks)
        {
          if (m_selector(element))
            {
              double plane_point[3] = {2,0,0};
              //double plane_normal[3] = {1, .5, -.5};
              double plane_normal[3] = {1, 0, 0};

#if 0
              double proj_pt_0[3]={0,0,0};
              double proj_pt_1[3]={0,0,0};
              project_point_to_plane(plane_point, plane_normal, coord0, proj_pt_0);
              project_point_to_plane(plane_point, plane_normal, coord1, proj_pt_1);
#endif

              double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
              double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

              double v01[3] = {0,0,0};
              difference(v01, coord1, coord0);
              normalize(v01);
              normalize(plane_normal);
              double v01dotn = std::abs(dot(v01, plane_normal));

              double d01p = std::abs(dot_0)+std::abs(dot_1);
              dot_0 = shock_function(dot_0);
              dot_1 = shock_function(dot_1);

              double *fd0 = m_eMesh.field_data(m_nodal_refine_field, node0);
              double *fd1 = m_eMesh.field_data(m_nodal_refine_field, node1);
              fd0[0] = dot_0;
              fd1[0] = dot_1;
              
              double d01 = distance(coord0, coord1);

              if ( (1 + 0*v01dotn + 0*(d01p/(d01+1.e-10)))*std::abs(dot_0 - dot_1) > 0.1)  // 0.2
                {
                  return true;
                }
            }
          return false;
        }

        //double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        //return m_selector(entity) && fdata[0] > 0;
      };

      struct ShockBasedUnrefinePredicate : public EdgeBasedMarkerPredicate {
      
        ShockBasedUnrefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
          EdgeBasedMarkerPredicate(selector, field, tolerance) {}

        /// Return true for unrefine, false for ignore
        bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                        double *coord0, double *coord1)
        {
          if (m_selector(element))
            {
              // do something
            }
          return false;
        }

      };

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_shock)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.addField("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.addField("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* unrefine_field = eMesh.addField("unrefine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.addField("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);

                SetUnrefineField set_unref_field(eMesh);
                eMesh.elementOpLoop(set_unref_field, unrefine_field);
              }

            eMesh.saveAs( output_files_loc+"local_tet_N_5_EdgeBased_shock_0_"+post_fix[p_size]+".e");

            stk::mesh::Selector univ_selector(eMesh.getFEM_meta_data()->universal_part());

            PredicateBasedEdgeMarker<ShockBasedRefinePredicate, ShockBasedUnrefinePredicate>
              breaker(ShockBasedRefinePredicate(nodal_refine_field, eMesh, univ_selector, refine_field, 0.0),
                      ShockBasedUnrefinePredicate(univ_selector, unrefine_field, 0.0),
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 3; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);
                //eMesh.elementOpLoop(set_unref_field, unrefine_field);

                std::cout << "P[" << eMesh.getRank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.getRank() << "] done... ipass= " << ipass << std::endl;
                eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_ipass"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.getRank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.saveAs(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_ipass"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.saveAs( output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }



    }
  }
}
