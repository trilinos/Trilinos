/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
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
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>

namespace stk
{
  namespace percept
  {
    namespace regression_tests
    {

#define EXTRA_PRINT 0

      static double pressure_value = 123.4;

#if 1
      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
      STKUNIT_UNIT_TEST(perceptMesh, open_new_close_PerceptMesh)
      {
        EXCEPTWATCH;

        // start_demo_open_new_close_PerceptMesh
        PerceptMesh eMesh(3u);
        eMesh.newMesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

        mesh::FieldBase* pressure_field = eMesh.addField("pressure", stk::mesh::fem::FEMMetaData::NODE_RANK, scalarDimension);
        eMesh.addField("velocity", stk::mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
        eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();

        // create a field function from the new pressure field
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        //if (eMesh.getRank()== 0) eMesh.printFields("Pressure");
        //exit(1);

        // here we could evaluate this field function
        double x=0.123, y=0.234, z=0.345, time=0.0;
        std::cout << "P[" << eMesh.getRank() << "] "
                  << "before write ff_pressure = " << eval(x,y,z,time, ff_pressure) << std::endl;

        //evalPrint(x, y, z, time, ff_pressure);

        double pval = eval(x, y, z, time, ff_pressure);
        EXPECT_DOUBLE_EQ(pval, pressure_value);

        eMesh.saveAs("./output_files/cube_with_pressure.e");
        eMesh.close();

        // end_demo

        // start_demo_open_new_close_PerceptMesh_1
        // open the file we previously saved with the new fields
        eMesh.openReadOnly("./input_files/cube_with_pressure.e");

        // get the pressure field
        pressure_field = eMesh.getField("pressure");

        // FIXME
        std::vector< const mesh::FieldBase * > sync_fields( 1 , pressure_field );
        mesh::communicate_field_data( eMesh.getBulkData()->shared_aura() , sync_fields );
        // FIXME

        //if (1 || eMesh.getRank()== 0) eMesh.printFields("Pressure");
     
        FieldFunction ff_pressure_1("ff_pressure", pressure_field, eMesh, 3, 1);
        ff_pressure_1.addAlias("P");
        StringFunction sf_pressure("P");
        std::cout << "P[" << eMesh.getRank() << "] "
                  << "after read ff_pressure = " << eval(x,y,z,time, ff_pressure_1) << std::endl;

        // a point-source at the origin
#define EXACT_SOL log(sqrt( x*x + y*y + z*z) + 1.e-10)        

        StringFunction sf_exact_solution(EXPAND_AND_QUOTE(EXACT_SOL), Name("sf_exact_solution"), 3, 1);
        StringFunction sf_error = sf_exact_solution - sf_pressure;
        
        std::cout << "P[" << eMesh.getRank() << "] "
                  << "sf_pressure = " << eval(x,y,z,time, sf_pressure) << std::endl;
        //!evalPrint(x,y,z,time, sf_error);
        std::cout << "P[" << eMesh.getRank() << "] "
                  << "sf_error = " << eval(x,y,z,time, sf_error) << std::endl;
        double val_cpp = EXACT_SOL - pressure_value;
        double val_sf  = eval(x,y,z,time, sf_error);
        EXPECT_DOUBLE_EQ(val_sf, val_cpp);
        // end_demo

      }

#endif

      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
      STKUNIT_UNIT_TEST(perceptMesh, open_new_close_PerceptMesh_2)
      {
        EXCEPTWATCH;

        double x=0.123, y=0.234, z=0.345, time=0.0;

        // start_demo_open_new_close_PerceptMesh_2
        PerceptMesh eMesh(3u);
        // open the file we previously saved with the new fields
        eMesh.openReadOnly("./input_files/cube_with_pressure.e");

        eMesh.printInfo("Info after reading mesh");

        //eMesh.printFields();
        mesh::FieldBase *f_coords = eMesh.getField("coordinates");

        // create a field function from the existing coordinates field
        FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

        // here we could evaluate this field function
        evalVec3Print(x, y, z, time, ff_coords);

        // get the pressure field
        mesh::FieldBase* pressure_field = eMesh.getField("pressure");

        // FIXME
        std::vector< const mesh::FieldBase * > sync_fields( 1 , pressure_field );
        mesh::communicate_field_data( eMesh.getBulkData()->shared_aura() , sync_fields );
        // FIXME

        //double * pdata = eMesh.node_field_data(pressure_field, 1);

        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);
        ff_pressure.addAlias("P");
        StringFunction sf_pressure("P");

        // a point-source at the origin
#define EXACT_SOL log(sqrt( x*x + y*y + z*z) + 1.e-10)        

        StringFunction sf_exact_solution(EXPAND_AND_QUOTE(EXACT_SOL), Name("sf_exact_solution"), 3, 1);
        StringFunction sf_error = sf_exact_solution - sf_pressure;
        
        evalPrint(x,y,z,time, sf_error);
        double val_cpp = EXACT_SOL - pressure_value;
        double val_sf  = eval(x,y,z,time, sf_error);
        EXPECT_DOUBLE_EQ(val_sf, val_cpp);
        // end_demo

      }

#define EXPECT_CATCH(expression, name)                                  \
      {                                                                 \
        bool didCatch_ ## name = false;                                 \
        try {                                                           \
          expression ;                                                  \
        }                                                               \
        catch ( const std::exception & X ) {                            \
          std::cout << # name << "  expected to catch this exception: " << X.what() << std::endl; \
          didCatch_ ## name = true;                                     \
        }                                                               \
        catch( ... ) {                                                  \
          std::cout << " Caught unknown exception"                      \
                    << std::endl ;                                      \
          std::cout.flush();                                            \
          didCatch_ ## name = false;                                    \
        }                                                               \
        EXPECT_TRUE(didCatch_ ## name);                                 \
      }

      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
      STKUNIT_UNIT_TEST(perceptMesh, open_new_close_PerceptMesh_3)
      {
        EXCEPTWATCH;

        // start_demo_open_new_close_PerceptMesh_3
        PerceptMesh eMesh(3u);
        eMesh.newMesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

        mesh::FieldBase* pressure_field = eMesh.addField("pressure", stk::mesh::fem::FEMMetaData::NODE_RANK, scalarDimension);
        eMesh.addField("velocity", stk::mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
        eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();

        EXPECT_CATCH( eMesh.commit() , commit_again);

        // create a field function from the new pressure field
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        // save
        eMesh.saveAs("./output_files/cube_with_pressure_3.e");
        eMesh.close();

        EXPECT_CATCH( eMesh.printInfo("bad", 1) , mesh_closed_try_print);

        // end_demo

      }


      //======================================================================================================================      
      //======================================================================================================================      
      //======================================================================================================================      
      STKUNIT_UNIT_TEST(perceptMesh, open_new_reopen_PerceptMesh)
      {
        EXCEPTWATCH;

        // start_demo_open_new_reopen_PerceptMesh
        PerceptMesh eMesh(3u);
        eMesh.newMesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

        eMesh.addField("pressure", stk::mesh::fem::FEMMetaData::NODE_RANK, scalarDimension);
        eMesh.addField("velocity", stk::mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
        eMesh.addField("element_volume", eMesh.element_rank(), scalarDimension);

        eMesh.commit();

        /// reopen the mesh to allow for more fields to be added - note that this involves a db write/read operation
        eMesh.reopen("./output_files/optional_temp_filename.e");
        mesh::FieldBase* momentum_field = eMesh.addField("momentum", stk::mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
        eMesh.commit();

        // create a field function from the new pressure field
        mesh::FieldBase *pressure_field = eMesh.getField("pressure");
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        // set the momentum field
        std::vector<double> momentum_values(3);
        momentum_values[0] = 2034.5;
        momentum_values[1] = 2134.5;
        momentum_values[2] = 2234.5;
        ConstantFunctionVec initMomentumValues(momentum_values, "initMomVal");
        
        // create a field function from the new momentum field
        FieldFunction ff_momentum("ff_momentum", momentum_field, eMesh, 3, 3);
        ff_momentum.interpolateFrom(initMomentumValues);

        // save
        eMesh.saveAs("./output_files/cube_with_pressure_and_momentum.e");
        eMesh.close();


        // end_demo

      }


    }//    namespace unit_tests
  }//  namespace percept
}// namespace stk

