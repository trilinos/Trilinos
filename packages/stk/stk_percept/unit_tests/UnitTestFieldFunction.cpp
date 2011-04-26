/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <Parallel.hpp>
#include <math.h>

namespace stk {
namespace percept {
namespace unit_tests {

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, fieldFunction_demo_1_0_0)
{
  EXCEPTWATCH;

  // start_demo_fieldFunction_1
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
  eMesh.commit();
  eMesh.printInfo("fieldFunction_demo_1_0_0", 2);

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

  // here we could evaluate this field function
  double x=0.123, y=0.234, z=0.345, time=0.0;
  evalVec3Print(x, y, z, time, ff_coords);
  // end_demo

}

STKUNIT_UNIT_TEST(function, fieldFunction_read_print)
{
  EXCEPTWATCH;
  // just reads a mesh file and prints some info about the meta data

  bool printInfo = false;

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.newMeshReadOnly(PerceptMesh::GMeshSpec(config_mesh));
	
  mesh::fem::FEMMetaData& metaData = *eMesh.getFEM_meta_data();

  const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

  unsigned nparts = parts.size();
  if (printInfo) std::cout << "Number of parts = " << nparts << std::endl;

  // here's where we can add parts
  // ...
  // ... then we would have to commit the metaData

  const stk::mesh::FieldVector & fields =  metaData.get_fields();
  unsigned nfields = fields.size();
  if (printInfo)
  {
    std::cout << "Number of fields = " << fields.size() << std::endl;
    for (unsigned ifld = 0; ifld < nfields; ifld++)
    {
      stk::mesh::FieldBase *field = fields[ifld];
      if (printInfo) std::cout << "Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
      if (printInfo) std::cout << *field << std::endl;
      unsigned nfr = field->restrictions().size();
      if (printInfo) std::cout << " number of field restrictions= " << nfr << std::endl;
      for (unsigned ifr = 0; ifr < nfr; ifr++)
      {
        const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
        mesh::Part& frpart = metaData.get_part(fr.ordinal());
        if (printInfo) std::cout << " field restriction " << ifr << " stride[0] = " << fr.stride[0] <<
                         " type= " << fr.type() << " ord= " << fr.ordinal() <<
                         " which corresponds to Part= " << frpart.name() << std::endl;
      }
    }
  }
}


//=============================================================================
//=============================================================================
//=============================================================================

#define EXPR_COORD_MAG (sqrt(x*x + y*y + z*z))

class CheckCoordMag : public GenericFunction
{
public:
  bool m_error;
  std::string m_name;
  CheckCoordMag(std::string name="") : m_error(false), m_name(name) {}
  virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
  {
    double x = domain(0);
    double y = domain(1);
    double z = domain(2);
    double v = EXPR_COORD_MAG;
    double cmag_field_node = codomain(0);
    //STKUNIT_EXPECT_DOUBLE_EQ(v, cmag_field_node);
    if (fabs(v-cmag_field_node) > 1.e-6)
    {
      std::cout << "CheckCoordMag:: " << m_name <<
        " v= " << v << " x= " << x << " y= " << y << " z= "<< z << " cmag_field_node= " << cmag_field_node << std::endl;
      Util::pause(true, "cmag_field_node");
      STKUNIT_ASSERT_NEAR(v, cmag_field_node, 1.e-9);
      m_error = true;
    }
  }

};

//=============================================================================
//=============================================================================
//=============================================================================
STKUNIT_UNIT_TEST(function, fieldFunction_demo_1)
{
  EXCEPTWATCH;

  
  {
    stk::io::util::Gmesh_STKmesh_Fixture gms(MPI_COMM_WORLD, "3x3x3|bbox:0,0,0,1,1,1");
    std::cout << "gms= " << &gms << std::endl;
  }

  std::cout << "gms= end"  << std::endl;

  // start_demo_fieldFunction_1
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
  eMesh.commit();

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

  // here we could evaluate this field function
  double x=0.123, y=0.234, z=0.345, time=0.0;
  evalVec3Print(x, y, z, time, ff_coords);
  // end_demo

}

STKUNIT_UNIT_TEST(function, fieldFunction_demo_2)
{
  EXCEPTWATCH;

  // start_demo_fieldFunction_2
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec("3x3x3|bbox:0,0,0,1,1,1")); // create a 3x3x3 hex mesh in the unit cube

  // add a new field
  // NOTE: we have to create the fields here before committing the mesh
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.addField("coords_mag_field", mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
  eMesh.commit();

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  // get the new field created by PerceptMesh
  stk::mesh::FieldBase* coords_mag_field = eMesh.getField("coords_mag_field");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh,  3, 3);
  evalVec3Print(0.1,0.1,0.1,0.0, ff_coords);

  // create a StringFunction to define the magnitude of the coordinates
  StringFunction coords_mag_sf( "sqrt(x*x + y*y + z*z)" , Name("coords_mag_sf"), 3, 1);
  double x=0.123, y=0.234, z=0.345;
  double vv = std::sqrt(x*x + y*y + z*z);            // evaluate the expression in C++
  double v1 = eval(x, y, z, 0, coords_mag_sf);       // evaluate the analytic expression
  STKUNIT_ASSERT_NEAR(vv, v1, 1.e-9);                          // the two results should be the same

  // Interpolate the analytic field defined by "coords_mag_sf" to the field we created to hold the coordinate magnitude field
  // 1. create a field function to represent the new coordinate magnitude field, and interpolate the string function to its nodes
  FieldFunction coords_mag_field_function("coords_mag_field_function", coords_mag_field, eMesh, 3, 1);

  coords_mag_field_function.interpolateFrom(coords_mag_sf);

  // We can now write the model with the new coordinates magnitude field to an Exodus file
  eMesh.saveAs("./cube_hex8_withCoordMag_out.e");
  // end_demo

  // start_demo_fieldFunction_3

  // tell Percept that we want to refer to the ff_coords FieldFunction by a simple alias "mc"
  ff_coords.addAlias("mc");

  // define a new StringFunction that does the same thing as coords_mag_sf, evaluates the coordinate magnitudes
  StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), 3, 1);
  // end_demo

}

STKUNIT_UNIT_TEST(function, fieldFunction_readMesh_createField_interpolateFrom)
{
  EXCEPTWATCH;
  // more i/o; adding a field; writing the resulting mesh; creating a FieldFunction, invoking interpolateFrom

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec(config_mesh));
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.addField("coords_mag_field", mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
  eMesh.commit();

  unsigned p_rank = eMesh.getBulkData()->parallel_rank();
  //unsigned p_size = eMesh.getBulkData()->parallel_size();
  Util::setRank(p_rank);

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  /// get the new field created by readModelCreateOptionalFields()
  stk::mesh::FieldBase* coords_mag_field = eMesh.getField("coords_mag_field");
  VERIFY_OP_ON(coords_mag_field, !=, 0, "TEST::function::fieldFunction_readMesh_createField_interpolateFrom: null coords_mag_field");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3, FieldFunction::SIMPLE_SEARCH );

  /// here we could evaluate this field function
  if (0)
  {
    if (p_rank == 0)
    {
      std::cout << "TEST::function::fieldFunction_readMesh_createField_interpolateFrom eval ff_coords=" << std::endl;
      evalVec3Print(0.1, 0.2, 0.3, 0.0, ff_coords);
    }
  }

  StringFunction coords_mag_sf( EXPAND_AND_QUOTE(EXPR_COORD_MAG) , Name("coords_mag_sf"), 3, 1);

  /// create a field function to represent the new coordinate magnitude field, and interpolate the string function to its nodes
  FieldFunction coords_mag_field_function("coords_mag_field_function", coords_mag_field, eMesh, 3, 3, FieldFunction::SIMPLE_SEARCH );
  coords_mag_field_function.interpolateFrom(coords_mag_sf);

  /// check that the coordinates mag field is set correctly
  {
    EXCEPTWATCH;
    CheckCoordMag checkCoordMag;
    //if (!p_rank) std::cout << "checkCoordMag..." << std::endl;
    eMesh.nodalOpLoop(checkCoordMag, coords_mag_field);
    //if (!p_rank) std::cout << "checkCoordMag...done" << std::endl;
    STKUNIT_EXPECT_FALSE(checkCoordMag.m_error);
  }

  try {
    ff_coords.addAlias("mc");
    StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), Dimensions(3), Dimensions(1));

    double tol1 = 1.e-12;

    {
      MDArray vv = evalVec3(0.1, 0.2, 0.3, 0.0, ff_coords);

      STKUNIT_ASSERT_NEAR(vv(0), 0.1, tol1);
      STKUNIT_ASSERT_NEAR(vv(1), 0.2, tol1);
      STKUNIT_ASSERT_NEAR(vv(2), 0.3, tol1);
    }

    {
      double vv = eval(0.1, 0.2, 0.3, 0.0, sfcm);
      double v_expect = std::sqrt(0.1*0.1 + 0.2*0.2 + 0.3*0.3);
      STKUNIT_ASSERT_NEAR(vv, v_expect, tol1);
    }

    coords_mag_field_function.interpolateFrom(sfcm);
    CheckCoordMag checkCoordMag1(std::string(EXPAND_AND_QUOTE(__FILE__))+": "+toString(__LINE__));
    if (!p_rank) std::cout << "checkCoordMag1..." << std::endl;
    eMesh.nodalOpLoop(checkCoordMag1, coords_mag_field);
    if (!p_rank) std::cout << "checkCoordMag1...done" << std::endl;
    STKUNIT_EXPECT_FALSE(checkCoordMag1.m_error);
  }
  catch ( const std::exception * X ) {
    std::cout << "  unexpected exception: " << X->what() << std::endl;
    exit(123);
  }
  catch ( const std::exception & X ) {
    std::cout << "  unexpected exception: " << X.what() << std::endl;
    exit(124);
  }
  catch( ... ) {
    std::cout << "  ... exception" << std::endl;
    exit(125);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

enum {NPTS = 4};
static double testpoints[NPTS][4] = {
  {0.1234,     0.5678,    0.9,    0.812   },
  {0.1234e-3,  0.5678e-5, 0.97,   0.01    },
  {0.101,      0.02,      0.1020, 0.0122  },
  {0.003,      0.89,      0.01,   0.5     }
};

STKUNIT_UNIT_TEST(function, fieldFunction_multiplePoints)
{
  EXCEPTWATCH;
  std::cout << "TEST::function::fieldFunction_multiplePoints" <<  std::endl;
  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec(config_mesh));
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.addField("coords_mag_field", mesh::fem::FEMMetaData::NODE_RANK, vectorDimension);
  eMesh.commit();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );
  MDArray val1 = evalVec3(0.2, 0.3, 0.4, 0.0, ff_coords);
  std::cout << "val1= \n" << val1 << std::endl;

  MDArray points(NPTS, 3);
  MDArray output(NPTS, 3);
  MDArray output_expect(NPTS, 3);

  //StringFunction sf1("x+y*z");
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];
    points(ipts, 0) = x;
    points(ipts, 1) = y;
    points(ipts, 2) = z;
    //points(ipts, 3) = t;

    //std::cout << "field_op: ipts= " << ipts << std::endl;

    MDArray vec = evalVec3(x, y, z, t, ff_coords);
    STKUNIT_EXPECT_NEAR(vec(0), x, fabs(1.e-5*x));
    STKUNIT_EXPECT_NEAR(vec(1), y, fabs(1.e-5*y));
    STKUNIT_EXPECT_NEAR(vec(2), z, fabs(1.e-5*z));
    output_expect(ipts, 0) = x;
    output_expect(ipts, 1) = y;
    output_expect(ipts, 2) = z;
  }
  std::cout << "field_op: NPTS= " << NPTS << std::endl;
  //         ff_coords.setDomainDimensions(Dimensions(NPTS,3));
  //         ff_coords.setCodomainDimensions(Dimensions(NPTS,3));
  ff_coords.setDomainDimensions(Dimensions(3));
  ff_coords.setCodomainDimensions(Dimensions(3));
  ff_coords(points, output, 0.0);
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    STKUNIT_EXPECT_NEAR(output(ipts, 0), output_expect(ipts, 0), 1.e-5*(fabs(output_expect(ipts,0))) );
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, fieldFunction_point_eval_verify)
{
  EXCEPTWATCH;
  /// test evaluation of field function at a point

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec(config_mesh));

  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// here we evaluate this field function
  MDArray val1 = evalVec3(0.2, 0.3, 0.4, 0.0, ff_coords);
  //std::cout << "eval = \n" << val1 << std::endl;

  stk::mesh::BulkData& bulkData = *eMesh.getBulkData();


  bool didCatch = false;
  try {
    // evaluate a point that is known to be outside the domain
    MDArray val10 = evalVec3(1.2, 1.3, 1.4, 0.0, ff_coords);
  }
  catch ( const std::exception & X ) {
    std::cout << "  expected to catch this exception: " << X.what() << std::endl;
    didCatch = true;
  }
  catch( ... ) {
    std::cout << "  P:" << bulkData.parallel_rank()
              << " Caught unknown exception"
              << std::endl ;
    std::cout.flush();
    didCatch = false;
  }
  STKUNIT_EXPECT_TRUE(didCatch);

  //double value = eval(1.2, 2.3, 3.4, 0.0, ff_coords);
  MDArray pts(3);
  MDArray output_pts(3);
  pts(0) = 0.2; pts(1) = 0.3; pts(2) = 0.4; //pts(3) = 0.0;
  ff_coords(pts, output_pts);
  STKUNIT_ASSERT_NEAR(pts(0), output_pts(0), 1.e-9);
  STKUNIT_ASSERT_NEAR(pts(1), output_pts(1), 1.e-9);
  STKUNIT_ASSERT_NEAR(pts(2), output_pts(2), 1.e-9);
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, fieldFunction_point_eval_timing)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  /// test evaluation of field function at a point

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.newMesh(PerceptMesh::GMeshSpec(config_mesh));

  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  //unsigned p_rank = eMesh.getBulkData()->parallel_rank();
  unsigned p_size = eMesh.getBulkData()->parallel_size();
  // FIXME
  if (p_size > 1) return;
  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.getField("coordinates");

  for (unsigned iSearchType = 0; iSearchType < 2; iSearchType++)
  {
    /// create a field function from the existing coordinates field
    FieldFunction::SearchType search_type = (iSearchType == 0 ? FieldFunction::SIMPLE_SEARCH : FieldFunction::STK_SEARCH);
    //std::cout <<  "P[" << Util::getRank() <<  "] search_type = " << search_type << " = " << typeid(search_type).name() << std::endl;
    FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                            Dimensions(3), Dimensions(3), search_type
                            );

    // The first point that is evaluated fires the setup of the stk::search data structure (oct-tree, bih-tree)
    double t1st =  stk::wall_time();
    MDArray val11 = evalVec3(0.2, 0.3, 0.4, 0.0, ff_coords);
    val11 = evalVec3(0.2, 0.3, 0.4, 0.0, ff_coords);
    t1st = stk::wall_dtime(t1st);

    // timings
    unsigned numIter = 10000;
    //unsigned numIter = 1000;
    MDArray pts(3);
    MDArray output_pts(3);

    // ensure the same set of random data each run
    Teuchos::ScalarTraits<double>::seedrandom(12345);

    double tstart =  stk::wall_time();
    for (unsigned iter = 0; iter < numIter; iter++)
    {
      double rnd = Teuchos::ScalarTraits<double>::random();
      pts(0) = (rnd+1.0)/2.0;
      rnd = Teuchos::ScalarTraits<double>::random();
      pts(1) = (rnd+1.0)/2.0;
      rnd = Teuchos::ScalarTraits<double>::random();
      pts(2) = (rnd+1.0)/2.0;

      // FIXME
      //!! pts(0) = 0.2; pts(1) = 0.3; pts(2)= 0.4;
      // FIXME
      ff_coords(pts, output_pts, 0.0);
#if 0
      STKUNIT_EXPECT_DOUBLE_EQ(pts(0), output_pts(0));
      STKUNIT_EXPECT_DOUBLE_EQ(pts(1), output_pts(1));
      STKUNIT_EXPECT_DOUBLE_EQ(pts(2), output_pts(2));
#endif
    }

    double total_time = stk::wall_dtime(tstart);
    if (1 || EXTRA_PRINT) std::cout
                            << "TEST::function::fieldFunction_point_eval_timing: "
                            << " for search_type= " << (iSearchType==0?"SIMPLE_SEARCH":"STK_SEARCH")<< "\n"
                            << "    time for 1st eval=  " << t1st << "\n"
                            << "    for " << numIter << " iterations, evaluating field(x,y,z) time = " << total_time  << "\n"
                            << "    average per point lookup and eval time = " << (total_time/((double)numIter)) << std::endl;
  }
  //std::cout << "P[" << Util::getRank() <<  "] TEST::function::fieldFunction_point_eval_timing done " << std::endl;
}

#if 0
int main()
{
  StringFunction sf1("x + y + z + t");
  sf1(xyz, out);

  StringFunction sf2("x - y");
  sfx("x");
  sfy("y");
  sfxy("x-y");
  sfxy1== sfx-sfy;

  StringFunction sf21dif = sf2-sf1;

  mesh::fem::FEMMetaData m;
  //...  setup field, etc.
  FieldFunction ff1("ff1", part1, field_1);  // can be nodal or elemental

  ff1.interpolateFrom(sf1);

  FieldFunction ff2(ff1);  // copy
  FieldFunction ff3 = ff1;  // copy

  Function ff21diff = ff2-ff1;  // should be zero
  ZeroFunction zero_func;

  //assert( ff21diff == zero_func );
  random_probe_assert(ff21diff, zero_func, 100);

  DifferenceFunction dfsf21(sf1,sf2);
  //assert( sf21dif == dfsf21 );
  random_probe_assert( sf21dif, dfsf21, 100);

  // check copy
  ff2.interpolateFrom(sf2);
  ff2.assertMostlyEqual(ff1);
}
#endif

} // namespace unit_tests
} // namespace percept
} // namespace stk

