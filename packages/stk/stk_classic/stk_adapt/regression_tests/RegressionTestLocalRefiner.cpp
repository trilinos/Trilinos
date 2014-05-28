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
#include <stk_util/diag/PrintTable.hpp>

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

#include <stk_adapt/IEdgeBasedAdapterPredicate.hpp>
#include <stk_adapt/IElementBasedAdapterPredicate.hpp>
#include <stk_adapt/PredicateBasedElementAdapter.hpp>
#include <stk_adapt/PredicateBasedEdgeAdapter.hpp>

#include <stk_percept/function/ElementOp.hpp>

namespace stk
{
  namespace adapt
  {
    namespace regression_tests
    {

#if !defined(__IBMCPP__)
      bool DO_TESTS=true;
#else
      // maybe it takes too long on dawn so we turn off these long-running tests...
      bool DO_TESTS=false;
#endif

      bool LARGE_TEST_ONLY=false;

#include "RegressionTestFileLoc.hpp"

//       const std::string input_files_loc="./input_files/";
//       const std::string output_files_loc="./output_files/";

      static std::string post_fix[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};

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
          VectorFieldType* coordField = m_eMesh.get_coordinates_field();
                
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
          VectorFieldType* coordField = m_eMesh.get_coordinates_field();
                
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

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);

            eMesh.commit();

            SetRefineField set_ref_field(eMesh);
            eMesh.elementOpLoop(set_ref_field, refine_field);

            SetUnrefineField set_unref_field(eMesh);
            //eMesh.elementOpLoop(set_ref_field, refine_field);

            eMesh.save_as( output_files_loc+"local_tet_N_5_ElementBased_0_"+post_fix[p_size]+".e");

            ElementRefinePredicate erp(0, refine_field, 0.0);

            PredicateBasedElementAdapter<ElementRefinePredicate>
              breaker(erp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 3; ipass++)
              {
                eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                eMesh.elementOpLoop(set_unref_field, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_ElementBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      struct MyEdgeBasedRefinePredicate : public IEdgeBasedAdapterPredicate {
      
        MyEdgeBasedRefinePredicate(stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
          IEdgeBasedAdapterPredicate(selector, field, tolerance) {}

        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                       double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark = 0;
          if (0 == m_selector || (*m_selector)(element))
            {
              double plane_point[3] = {2,0,0};
              double plane_normal[3] = {1, 0, 0};

              double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
              double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

              // if edge crosses the plane...
              if (dot_0 * dot_1 < 0)
                {
                  mark |= DO_REFINE;
                }

              if (coord0[1] < 0 && coord1[1] < 0)
                {
                  mark |= DO_UNREFINE;
                }
            }

          return mark;
        }
      };

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_0_"+post_fix[p_size]+".e");

            MyEdgeBasedRefinePredicate mrp(0, refine_field, 0.0);

            PredicateBasedEdgeAdapter<MyEdgeBasedRefinePredicate>
              breaker(mrp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 5; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            //breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_"+post_fix[p_size]+".e");

#if 1
            for (int iunref_pass=0; iunref_pass < 5; iunref_pass++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      struct PlaneShock 
      {
        double plane_point_init[3]; // = {2 + shock_displacement,0,0};
        double plane_normal[3]; // = {1, 0, 0};
        double plane_point[3]; // = {2 + shock_displacement,0,0};

        PlaneShock()
        {
          plane_point_init[0]=0;
          plane_point_init[1]=0;
          plane_point_init[2]=0;
          plane_point[0]=0;
          plane_point[1]=0;
          plane_point[2]=0;
          plane_normal[0]=1;
          plane_normal[1]=0;
          plane_normal[2]=0;
        }

        void setCurrentPlanePoint(double shock_displacement)
        {
          normalize(plane_normal);
          plane_point[0] = plane_point_init[0] + shock_displacement*plane_normal[0];
          plane_point[1] = plane_point_init[1] + shock_displacement*plane_normal[1];
          plane_point[2] = plane_point_init[2] + shock_displacement*plane_normal[2];
        }
        
      };

      static double shock_width = 1./5.0;

      static double shock_function(double x)
      {
        // normalize by width
        double width = shock_width; // 1./5.0;
        x /= width;
        return std::tanh(x);
      }

      
      static double shock_diff(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh,
                               stk::mesh::Entity & node0, stk::mesh::Entity & node1, double *coord0, double *coord1, PlaneShock& shock, double shock_displacement)
      {
        shock.setCurrentPlanePoint(shock_displacement);
        double *plane_point = shock.plane_point;
        double *plane_normal = shock.plane_normal;

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

        if (nodal_refine_field)
          {
            double *fd0 = eMesh.field_data(nodal_refine_field, node0);
            double *fd1 = eMesh.field_data(nodal_refine_field, node1);
            fd0[0] = dot_0;
            fd1[0] = dot_1;
          }
              
        double d01 = distance(coord0, coord1);

        return (1 + 0*d01 + 0*d01p + 0*v01dotn)*std::abs(dot_0 - dot_1);
      }

      struct ShockBasedRefinePredicate : public IEdgeBasedAdapterPredicate {
      
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase * m_nodal_refine_field;
        PlaneShock m_shock;
        double m_shock_displacement;
        double m_shock_diff_criterion;

        ShockBasedRefinePredicate(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh, stk::mesh::Selector* selector, stk::mesh::FieldBase *field, double tolerance,
                                  PlaneShock shock, double shock_displacement=0, double shock_diff_criterion=0.4) :
          IEdgeBasedAdapterPredicate(selector, field, tolerance), m_eMesh(eMesh),m_nodal_refine_field(nodal_refine_field), m_shock(shock), m_shock_displacement(shock_displacement),
          m_shock_diff_criterion(shock_diff_criterion) {}


        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                        double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark=0;
          if (0 == m_selector || (*m_selector)(element))
            {

              // refine check
              double d01 = shock_diff(m_nodal_refine_field, m_eMesh, node0, node1, coord0, coord1, m_shock, m_shock_displacement);
              if ( d01 > m_shock_diff_criterion)  // 0.05, 0.2
                {
                  mark |= DO_REFINE;
                }

              // unrefine check
              d01 = shock_diff(0, m_eMesh, node0, node1, coord0, coord1, m_shock, m_shock_displacement);
              if ( d01 <= m_shock_diff_criterion/2.0 )  // 0.05, 0.2
                {
                  mark |= DO_UNREFINE;
                }

            }
          return mark;
        }

        //double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        //return m_selector(entity) && fdata[0] > 0;
      };


      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_shock)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            if (0)
              {
                SetRefineField set_ref_field(eMesh);
                eMesh.elementOpLoop(set_ref_field, refine_field);
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_shock_0_"+post_fix[p_size]+".e");

            PlaneShock shock;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            for (int ipass=0; ipass < 3; ipass++)
              {
                //eMesh.elementOpLoop(set_ref_field, refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_ipass_"+toString(ipass)+"_"+post_fix[p_size]+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_"+post_fix[p_size]+".e");

#if 0
            for (int iunref_pass=0; iunref_pass < 4; iunref_pass++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
                ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                breaker.unrefineTheseElements(elements_to_unref);
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_ipass_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
              }

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_shock_1_unref_"+post_fix[p_size]+".e");
#endif
            // end_demo
          }

      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./5.0;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            // start_demo_local_refiner_break_tet_to_tet_2

            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 2.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, 0.4);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            double delta_shock_displacement = 0.2;
            double shock_displacement = -2.0;  
            int num_ref_passes = 2;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_ref_istep_ipass_"+toString(istep)+"_"+toString(ipass)+"_"+post_fix[p_size]+".e");
                  }

                //breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << " unref list size= " << elements_to_unref.size() << std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_unref_istep_ipass_"+toString(istep)+"_"+toString(iunref_pass)+"_"+post_fix[p_size]+".e");
                  }
                
                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"final_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                breaker.unrefineAll();
              }
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }

      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock)
      {
        //if (1) return;
        const bool do_full_demo = false;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 3;  // 10 for stress testing
            bool save_intermediate=true;
            do_moving_shock_test(num_time_steps, save_intermediate);
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test_cyl_sidesets(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        shock_width = 1./25.0;
        double shock_diff_criterion = 0.04;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_tet4_0.e");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2;
            double shock_displacement = 0.0;  
            int num_ref_passes = 3;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                //breaker.deleteParentElements();
                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_"+post_fix[p_size]+".e");
                //exit(123);

                //breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    breaker.unrefineTheseElements(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }
                
//                 breaker.deleteParentElements();
                 eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_"+post_fix[p_size]+".e");
                 //exit(123);

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"cyl_sidesets_moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                breaker.unrefineAll();
                eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix[p_size]+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
              }

            if (0)
              breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"cyl_sidesets_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_cyl_sidesets)
      {
        const bool do_full_demo = false;
        if (LARGE_TEST_ONLY || !DO_TESTS) return;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_cyl_sidesets(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 1;  // 10 for stress testing
            do_moving_shock_test_cyl_sidesets(num_time_steps);
          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test_large_test(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //         shock_width = 1./50.; //1./25.0;
        //         shock_diff_criterion = 0.1;
        shock_width = 1./500.; //1./25.0;
        double shock_diff_criterion = 0.1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 2 || p_size == 8)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"large_testh1tet.g");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* refine_field       = eMesh.add_field("refine_field", eMesh.element_rank(), scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            //eMesh.delete_side_sets();

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"large_scale_moving_shock_"+post_fix[p_size]+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = -0.07;
            shock.plane_point_init[2] = -0.35;
            shock.plane_normal[0] = 0;
            shock.plane_normal[1] = .25/std::sqrt(.25*.25+1.);
            shock.plane_normal[2] = 1./std::sqrt(.25*.25+1.);

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, &univ_selector, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp, 
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            //breaker.setIgnoreSideSets(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2*0.08*(10./((double)num_time_steps));
            double shock_displacement = 0.0;  
            int num_ref_passes = 2;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                eMesh.save_as(output_files_loc+"large_scale_moving_shock_"+post_fix[p_size]+"_ref.e");
                //exit(123);

                breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }
                
                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"large_scale_moving_shock_"+post_fix[p_size]+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"large_scale_final_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                breaker.unrefineAll();
              }
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"large_scale_final_unrefed_moving_shock_"+post_fix[p_size]+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      STKUNIT_UNIT_TEST(regr_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_large_test)
      {
        const bool do_full_demo = false;
        if (!LARGE_TEST_ONLY || !DO_TESTS) return;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_large_test(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 30;  // 10 for stress testing
            do_moving_shock_test_large_test(num_time_steps, true);
            //int num_time_steps = 3;  // 10 for stress testing
            //do_moving_shock_test_large_test(num_time_steps);
          }
      }

    }
  }
}
