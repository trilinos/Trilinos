/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_percept/Percept.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/norm/Norm.hpp>
#include <stk_percept/math/Math.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/fixtures/Fixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace stk {
namespace percept {
namespace unit_tests {


#if 0
static stk::diag::Writer &
dw()
{
  //static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);
  int dw_enabled = 1;
  static stk::diag::Writer s_diagWriter(std::cout.rdbuf(), dw_enabled);

  s_diagWriter.setPrintMask(LOG_NORM+LOG_ALWAYS);
  return s_diagWriter;
}
#endif

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================
struct LocalFixture
{
  PerceptMesh eMesh;
  int bogus_init;
  mesh::fem::FEMMetaData& metaData;
  mesh::BulkData& bulkData;
  mesh::FieldBase *coords_field;
  StringFunction sfx;
  ConstantFunction sfx_res;

  LocalFixture(size_t num_xyz = 4, size_t num_y=0, size_t num_z=0, bool sidesets=false) : eMesh(3u), bogus_init(init(num_xyz, num_y, num_z, sidesets)),
                                                                     metaData(*eMesh.get_fem_meta_data()), bulkData(*eMesh.get_bulk_data()),
                                                                     coords_field( metaData.get_field<mesh::FieldBase>("coordinates") ),
                                                                     sfx("x", Name("sfx"), Dimensions(3), Dimensions(1) ),
                                                                     sfx_res (0.0, "sfx_res")

  {

  }

  int init(size_t num_xyz, size_t num_y_arg, size_t num_z_arg, bool sidesets=false)
  {
    // Need a symmetric mesh around the origin for some of the tests below to work correctly (i.e. have analytic solutions)
    const size_t num_x = num_xyz;
    const size_t num_y = num_y_arg? num_y_arg : num_xyz;
    const size_t num_z = num_z_arg? num_z_arg : num_xyz;
    std::string config_mesh =
      Ioss::Utils::to_string(num_x) + "x" +
      Ioss::Utils::to_string(num_y) + "x" +
      Ioss::Utils::to_string(num_z) + "|bbox:-0.5,-0.5,-0.5,0.5,0.5,0.5";
    if (sidesets) config_mesh += "|sideset:xXyYzZ";
	
    eMesh.new_mesh(GMeshSpec(config_mesh));
    eMesh.commit();
    return 1;
  }

};

//=============================================================================
//=============================================================================
//=============================================================================

//stk::diag::WriterThrowSafe _write_throw_safe(dw());
//dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));
//dw().setPrintMask(LOG_NORM+LOG_ALWAYS);

STKUNIT_UNIT_TEST(norm, volume)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  dw().m(LOG_NORM) << "TEST::norm::volume " << stk::diag::dendl;

  LocalFixture fix(3,3,12);
  mesh::fem::FEMMetaData& metaData = fix.metaData;
  mesh::BulkData& bulkData = fix.bulkData;
  PerceptMesh& eMesh = fix.eMesh;

  mesh::FieldBase *coords_field = fix.coords_field;

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// the function to be integrated - here it is just the identity, and when integrated should produce the volume
  ConstantFunction identity(1.0, "identity");

  /// A place to hold the result.
  /// This is a "writable" function (we may want to make this explicit - StringFunctions are not writable; FieldFunctions are
  /// since we interpolate values to them from other functions).
  ConstantFunction sqrt_volume(0.0, "sqrt_volume");

  /// Create the operator that will do the work
  Norm<2> l2Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  /// get the l2 norm of identity
  l2Norm(identity, sqrt_volume);

  Teuchos::ScalarTraits<double>::seedrandom(12345);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, sqrt_volume.getValue());

  int niter=1;
  for (int iter=0; iter < niter; iter++)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, eval(Math::random01(), Math::random01(), Math::random01(), 0.0, sqrt_volume));
  }

  //// rotate the mesh

#define DO_IO_TESTING 1

#if DO_IO_TESTING
  if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_original_out.e ..." << std::endl;
  eMesh.save_as("./gmesh_hex8_original_out.e");
  if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_original_out.e done" << std::endl;
#endif

  if (1)
  {
    Math::Matrix rmx = Math::rotationMatrix(0, 30);
    Math::Matrix rmy = Math::rotationMatrix(1, -45);
    Math::Matrix rmz = Math::rotationMatrix(2, 30);
    Math::Matrix rm;
    rm =  rmy * rmz;
    rm =  rmx * rm;
    MeshTransformer meshRotate(rm);
    eMesh.nodalOpLoop(meshRotate, coords_field);

    // for testing
#if DO_IO_TESTING
    if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_rotated_out.e ..." << std::endl;
    eMesh.save_as("./gmesh_hex8_rotated_out.e");
    if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_rotated_out.e done" << std::endl;
#endif

    l2Norm(identity, sqrt_volume);

    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, sqrt_volume.getValue());
  }

  //// scale the mesh
  if (1)
  {

    double scx = M_PI;
    double scy = M_E;
    double scz = std::sqrt(3.0);
    double sc=scx*scy*scz;
    Math::Matrix smx = Math::scalingMatrix(0, scx);
    Math::Matrix smy = Math::scalingMatrix(1, scy);
    Math::Matrix smz = Math::scalingMatrix(2, scz);
    Math::Matrix sm;
    sm =  smy * smz;
    sm =  smx * sm;
    //std::cout << "sm= " << sm << std::endl;
    MeshTransformer meshScale(sm);
    eMesh.nodalOpLoop(meshScale, coords_field);

    // for testing
#if DO_IO_TESTING
    if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_scaled_out.e ..." << std::endl;
    eMesh.save_as("./gmesh_hex8_scaled_out.e");
    if (1 || EXTRA_PRINT) std::cout << "TEST.norm.volume: writing gmesh_hex8_scaled_out.e done" << std::endl;
#endif

    l2Norm(identity, sqrt_volume);

    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(std::sqrt(sc), sqrt_volume.getValue());
  }

  //Function coords_l2_norm = ff_coords.norm_l2();

#if 0
  /// OPTION 1 - separate classes for operations on fields
  Norm_L2<VectorField> coords_l2_norm("coords_l2_norm");
  FieldFunction result = coords_l2_norm(ff_coords);


  /// OPTION 1a - just showing generic base class
  /// A FunctionOperator takes a FieldFunction and returns a FieldFunction
  FunctionOperator& coords_l2_norm_op = Norm_L2<VectorField>("coords_l2_norm");
  FieldFunction coords_l2_norm = coords_l2_norm_op(ff_coords);

  /// OPTION 2 - simple set of methods on FieldFunction or Function (depending on genericity)
  FieldFunction coords_l2_norm                           = ff_coords.norm_l2();
  FieldFunction coords_gradient_tensor_field_nodewise    = ff_coords.nodewise_gradient();
  FieldFunction coords_mag_nodewise_field                = ff_coords.nodewise_magnitude();
  FieldFunction coords_gradient_tensor_field_elementwise = ff_coords.elementwise_gradient();
#endif
  /// Note: need to create new fields each time, which requires a change to the meta data

}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(norm, surface_area)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  bool sidesets=true;
  LocalFixture fix(3,3,12,sidesets);
  //mesh::fem::FEMMetaData& metaData = fix.metaData;
  mesh::BulkData& bulkData = fix.bulkData;
  PerceptMesh& eMesh = fix.eMesh;
  //eMesh.save_as("junk.123.e");

  mesh::FieldBase *coords_field = fix.coords_field;

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// the function to be integrated - here it is just the identity, and when integrated should produce the area of faces
  ConstantFunction identity(1.0, "identity");

  /// A place to hold the result.
  /// This is a "writable" function (we may want to make this explicit - StringFunctions are not writable; FieldFunctions are
  /// since we interpolate values to them from other functions).
  ConstantFunction sqrt_area(0.0, "sqrt_area");

  /// Create the operator that will do the work
  stk::mesh::Part *surface_part = eMesh.get_part("surface_1");
  stk::mesh::Selector selector(*surface_part);
  bool is_surface_norm = true;
  Norm<2> l2Norm(bulkData, &selector, TURBO_NONE, is_surface_norm);
  /// get the l2 norm of identity
  l2Norm(identity, sqrt_area);

  Teuchos::ScalarTraits<double>::seedrandom(12345);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, sqrt_area.getValue());

  int niter=1;
  for (int iter=0; iter < niter; iter++)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, eval(Math::random01(), Math::random01(), Math::random01(), 0.0, sqrt_area));
  }

  //// rotate the mesh

#if DO_IO_TESTING
  eMesh.save_as("./gmesh_hex8_area_original_out.e");
#endif

  if (1)
  {
    Math::Matrix rmx = Math::rotationMatrix(0, 30);
    Math::Matrix rmy = Math::rotationMatrix(1, -45);
    Math::Matrix rmz = Math::rotationMatrix(2, 30);
    Math::Matrix rm;
    rm =  rmy * rmz;
    rm =  rmx * rm;
    MeshTransformer meshRotate(rm);
    eMesh.nodalOpLoop(meshRotate, coords_field);

    // for testing
#if DO_IO_TESTING
    eMesh.save_as("./gmesh_hex8_area_rotated_out.e");
#endif

    l2Norm(identity, sqrt_area);

    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(1.0, sqrt_area.getValue());
  }

  //// scale the mesh
  if (1)
  {

    double scx = M_PI;
    double scy = M_PI;
    double scz = M_PI;
    double sc=scx*scy;
    Math::Matrix smx = Math::scalingMatrix(0, scx);
    Math::Matrix smy = Math::scalingMatrix(1, scy);
    Math::Matrix smz = Math::scalingMatrix(2, scz);
    Math::Matrix sm;
    sm =  smy * smz;
    sm =  smx * sm;
    //std::cout << "sm= " << sm << std::endl;
    MeshTransformer meshScale(sm);
    eMesh.nodalOpLoop(meshScale, coords_field);

    // for testing
#if DO_IO_TESTING
    eMesh.save_as("./gmesh_hex8_area_scaled_out.e");
#endif

    l2Norm(identity, sqrt_area);
    
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(std::sqrt(sc), sqrt_area.getValue());
  }

  //Function coords_l2_norm = ff_coords.norm_l2();

}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(norm, string_function)
{
  EXCEPTWATCH;
  //stk::diag::WriterThrowSafe _write_throw_safe(dw());
  //dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));
  //dw().setPrintMask(LOG_NORM+LOG_ALWAYS);

  dw().m(LOG_NORM) << "TEST.norm.string_function " << stk::diag::dendl;

  LocalFixture     fix(4);
  mesh::fem::FEMMetaData&        metaData     = fix.metaData;
  mesh::BulkData&        bulkData     = fix.bulkData;
  PerceptMesh&        eMesh     = fix.eMesh;
  mesh::FieldBase*       coords_field = fix.coords_field;
  StringFunction   sfx          = fix.sfx;
  ConstantFunction sfx_res      = fix.sfx_res;

  /// Create the operator that will do the work
  /// get the l2 norm
  Norm<2> l2Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  l2Norm(sfx, sfx_res);

  double sfx_expect = std::sqrt(0.25/3.);
  Teuchos::ScalarTraits<double>::seedrandom(12345);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res.getValue());

  int niter=1;
  for (int iter=0; iter < niter; iter++)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, eval(Math::random01(), Math::random01(), Math::random01(), 0.0, sfx_res));
  }

  /// the function to be integrated:  (Integral[ abs(x), dxdydz]) =?= (2 * |x|^2/2 @ [0, 0.5]) ==> .25)
  Norm<1> l1Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  l1Norm(sfx, sfx_res);

  sfx_expect = 0.25;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());

  /// the function to be integrated:  (Max[ x^2+y^3+z^4, dxdydz]) =?= (@ [-0.5, 0.5]^3 ) ==> .5^2+.5^3+.5^4)
  StringFunction sfmax("x^2 + y^3 + z^4", Name("sfmax"), Dimensions(3), Dimensions(1) );
  Norm<-1> lInfNorm(bulkData, &metaData.universal_part(), TURBO_NONE);
  lInfNorm.setCubDegree(10);
  lInfNorm(sfmax, sfx_res);
  double sf1=eval(.5,.5,.5,0.0, sfmax);
  sfx_expect = 0.5*0.5 + 0.5*0.5*0.5 + 0.5*0.5*0.5*0.5;
  std::cout << "sfmax= " << sf1 << " sfx_expect= " << sfx_expect << " sfx_res= " << sfx_res.getValue() << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());

  /// the function to be integrated:  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest1.py)
  StringFunction sfxyz("x*y*z", Name("sfxyz"), Dimensions(3), Dimensions(1) );
  l2Norm(sfxyz, sfx_res);
  sfx_expect = 0.0240562612162344;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());

  /// the function to be integrated (but over a rotated domain):  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py)
  /// now rotate the mesh
  Math::Matrix rmz = Math::rotationMatrix(2, 30);
  Math::Matrix rm = rmz;
  MeshTransformer meshRotate(rm);
  eMesh.nodalOpLoop(meshRotate, coords_field);

  l2Norm(sfxyz, sfx_res);
  sfx_expect = 0.0178406008037016;
  // NOTE: we need extra quadrature accuracy to reproduce this result (cubDegree==4 in IntegratedOp almost gets it right)
  //   for now, we are satisfied with 2-3 digits
  //STKUNIT_EXPECT_DOUBLE_EQ(sfx_res.getValue(), sfx_expect);
  if (std::fabs(sfx_res.getValue()-sfx_expect) > 0.01*sfx_expect)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    STKUNIT_EXPECT_TRUE(false);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

/// This test uses a back door to the function that passes in the element to avoid the lookup of the element when the
///  StringFunction contains references to FieldFunctions
void TEST_norm_string_function_turbo_verify_correctness(TurboOption turboOpt)
{
  EXCEPTWATCH;
  //stk::diag::WriterThrowSafe _write_throw_safe(dw());
  //dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));
  //dw().setPrintMask(LOG_NORM+LOG_ALWAYS);

  // not ready yet, just return
  //if(1) return;

  dw().m(LOG_NORM) << "TEST.norm.string_function " << stk::diag::dendl;

  LocalFixture     fix(4);
  mesh::fem::FEMMetaData&        metaData     = fix.metaData;
  mesh::BulkData&        bulkData     = fix.bulkData;
  PerceptMesh&        eMesh     = fix.eMesh;
  mesh::FieldBase*       coords_field = fix.coords_field;
  StringFunction   sfx          = fix.sfx;
  ConstantFunction sfx_res      = fix.sfx_res;

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// the function to be integrated:  sqrt(Integral[x^2, dxdydz]) =?= sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)
  //StringFunction sfx("x", Name("sfx"), Dimensions(3), Dimensions(1) );

  ff_coords.add_alias("mc");
  //StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), Dimensions(3), Dimensions(1));
  StringFunction sfx_mc("mc[0]", Name("sfx_mc"), Dimensions(3), Dimensions(1) );
  StringFunction sfx_mc1("mc[0]", Name("sfx_mc1"), Dimensions(3), Dimensions(1) );

  if (1)
  {


    double x = -0.49+0.98*Math::random01();
    double y = -0.49+0.98*Math::random01();
    double z = -0.49+0.98*Math::random01();

    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(eval(x,y,z,0.0, sfx), eval(x,y,z,0.0, sfx_mc));
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX(eval(.034,0,0,0.0, sfx), eval(.034,0,0,0.0, sfx_mc));
  }

  /// A place to hold the result.
  /// This is a "writable" function (we may want to make this explicit - StringFunctions are not writable; FieldFunctions are
  /// since we interpolate values to them from other functions).
  ConstantFunction sfx_res_turbo(0.0, "sfx_res_turbo");
  ConstantFunction sfx_res_slow(0.0, "sfx_res_slow");
  ConstantFunction sfx_res_fast(0.0, "sfx_res_fast");
  ConstantFunction sfx_res_bucket(0.0, "sfx_res_bucket");

  // STATE m_element....

  /// Create the operator that will do the work
  /// get the l2 norm
  Norm<2> l2Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  l2Norm(sfx, sfx_res);

  Norm<2> l2Norm_turbo(bulkData, &metaData.universal_part(), turboOpt);
  Norm<2> l2Norm_turbo_bucket(bulkData, &metaData.universal_part(), TURBO_BUCKET);
  Norm<2> l2Norm_turbo1(bulkData, &metaData.universal_part(), turboOpt);
  l2Norm_turbo(sfx, sfx_res_turbo);

  l2Norm_turbo1(sfx_mc1, sfx_res_fast);
  l2Norm_turbo_bucket(sfx_mc1, sfx_res_bucket);

  l2Norm(sfx_mc, sfx_res_slow);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(eval(.023,0,0,0.0, sfx), eval(.023,0,0,0.0, sfx_mc));

  double sfx_expect = std::sqrt(0.25/3.);
  std::cout << "sfx_expect= " << sfx_expect << std::endl;

  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_bucket.getValue());

  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res.getValue());

  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_turbo.getValue());
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_slow.getValue());
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_fast.getValue()); //!
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res.getValue());

  //Util::pause(true, "13a");
  //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_res_turbo.getValue(), sfx_expect);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(sfx_expect, sfx_res_turbo.getValue(),  1.e-8);
  //Util::pause(true, "13");

  /// the function to be integrated:  (Integral[ abs(x), dxdydz]) =?= (2 * |x|^2/2 @ [0, 0.5]) ==> .25)
  Norm<1> l1Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  l1Norm(sfx, sfx_res);

  Norm<1> l1Norm_turbo(bulkData, &metaData.universal_part(), TURBO_NONE);
  l1Norm_turbo(sfx, sfx_res_turbo);
  l1Norm_turbo(sfx_mc, sfx_res_fast);
  l1Norm(sfx_mc, sfx_res_slow);

  sfx_expect = 0.25;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());

  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_slow.getValue());
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_fast.getValue());

  //// ----- here
  /// the function to be integrated:  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest1.py)
  StringFunction sfxyz("x*y*z", Name("sfxyz"), Dimensions(3), Dimensions(1) );
  l2Norm(sfxyz, sfx_res);
  l2Norm_turbo(sfxyz, sfx_res_turbo);
  sfx_expect = 0.0240562612162344;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());


  /// the function to be integrated (but over a rotated domain):  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py)
  /// now rotate the mesh
  Math::Matrix rmz = Math::rotationMatrix(2, 30);
  Math::Matrix rm = rmz;
  MeshTransformer meshRotate(rm);
  eMesh.nodalOpLoop(meshRotate, coords_field);

  l2Norm(sfxyz, sfx_res);
  l2Norm_turbo(sfxyz, sfx_res_turbo);
  sfx_expect = 0.0178406008037016;
  // NOTE: we need extra quadrature accuracy to reproduce this result (cubDegree==4 in IntegratedOp almost gets it right)
  //   for now, we are satisfied with 3 digits
  //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_res.getValue(), sfx_expect);
  if (std::fabs(sfx_res.getValue()-sfx_expect) > 0.01*sfx_expect)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    STKUNIT_EXPECT_TRUE(false);
  }
  if (std::fabs(sfx_res_turbo.getValue()-sfx_expect) > 0.01*sfx_expect)
  {
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());
    STKUNIT_EXPECT_TRUE(false);
  }
}

STKUNIT_UNIT_TEST(norm, string_function_turbo_verify_correctness_element)
{
  EXCEPTWATCH;

  TEST_norm_string_function_turbo_verify_correctness(TURBO_ELEMENT);
}

STKUNIT_UNIT_TEST(norm, string_function_turbo_verify_correctness_bucket)
{
  EXCEPTWATCH;

  TEST_norm_string_function_turbo_verify_correctness(TURBO_BUCKET);
}

//=============================================================================
//=============================================================================
//=============================================================================

/// This test uses a back door to the function that passes in the element to avoid the lookup of the element when the
///  StringFunction contains references to FieldFunctions
void TEST_norm_string_function_turbo_timings(TurboOption turboOpt)
{
  EXCEPTWATCH;
  //stk::diag::WriterThrowSafe _write_throw_safe(dw());
  //dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));
  //dw().setPrintMask(LOG_NORM+LOG_ALWAYS);

  // not ready yet, just return
  //if(1) return;

  dw().m(LOG_NORM) << "TEST.norm.string_function " << stk::diag::dendl;

  /// create a meta data/bulk data empty pair
  PerceptMesh eMesh(3u);

  if (1)
  {
    // Need a symmetric mesh around the origin for some of the tests below to work correctly (i.e. have analytic solutions)
    const size_t nxyz = 4;
    const size_t num_x = nxyz;
    const size_t num_y = nxyz;
    const size_t num_z = nxyz;
    std::string config_mesh =
      Ioss::Utils::to_string(num_x) + "x" +
      Ioss::Utils::to_string(num_y) + "x" +
      Ioss::Utils::to_string(num_z) + "|bbox:-0.5,-0.5,-0.5,0.5,0.5,0.5";
	
    eMesh.new_mesh(GMeshSpec(config_mesh));

    eMesh.commit();
  }

  mesh::fem::FEMMetaData& metaData = *eMesh.get_fem_meta_data();
  mesh::BulkData& bulkData = *eMesh.get_bulk_data();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  mesh::FieldBase *coords_field = metaData.get_field<mesh::FieldBase>("coordinates");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// the function to be integrated:  sqrt(Integral[x^2, dxdydz]) =?= sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)
  StringFunction sfx("x", Name("sfx"), Dimensions(3), Dimensions(1) );

  ff_coords.add_alias("mc");
  //StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), Dimensions(3), Dimensions(1));
  StringFunction sfx_mc("mc[0]", Name("sfx_mc"), Dimensions(3), Dimensions(1) );

  /// the function to be integrated:  sqrt(Integral[x^2, dxdydz]) =?= sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)

  /// A place to hold the result.
  /// This is a "writable" function (we may want to make this explicit - StringFunctions are not writable; FieldFunctions are
  /// since we interpolate values to them from other functions).
  ConstantFunction sfx_res(0.0, "sfx_res");
  ConstantFunction sfx_res_turbo(0.0, "sfx_res_turbo");
  ConstantFunction sfx_res_slow(0.0, "sfx_res_slow");
  ConstantFunction sfx_res_fast(0.0, "sfx_res_fast");

#define COL_SEP "|"
#define EXPR_CELL_WIDTH (80)

#define TIME_IT2(expr_none,expr_turbo,msg,topt)                         \
  {                                                                     \
    double TURBO_NONE_time    = 0;                                      \
    double TURBO_ON_time = 0;                                           \
    TIME_IT(expr_none,TURBO_NONE_time);                                 \
    TIME_IT(expr_turbo,TURBO_ON_time);                                  \
    if (1) std::cout << msg << #topt << "for expression= " << QUOTE(expr_none) << " timings= " << std::endl; \
    if (1) std::cout << "TURBO_NONE_time= " << TURBO_NONE_time << " "   \
                     << ( turboOpt==TURBO_ELEMENT?"TURBO_ELEMENT_time":"TURBO_BUCKET_time") <<"= " << TURBO_ON_time \
                     << " ratio= " << TURBO_NONE_time/TURBO_ON_time << std::endl; \
  }

  int numIter = 1;
  for (int iter = 0; iter < numIter; iter++)
  {
    /// Create the operator that will do the work
    /// get the l2 norm
    Norm<2> l2Norm      (bulkData, &metaData.universal_part(), TURBO_NONE);
    Norm<2> l2Norm_turbo(bulkData, &metaData.universal_part(), turboOpt);

    //double TURBO_ELEMENT_time=0;
    TIME_IT2(l2Norm(sfx, sfx_res); , l2Norm_turbo(sfx, sfx_res_turbo);, "Should be the same turboOpt= ", turboOpt );
    TIME_IT2(l2Norm(sfx_mc, sfx_res_slow); , l2Norm_turbo(sfx, sfx_res_fast); , "StringFunction with ref to FF, slow vs. fast" ,turboOpt );

    double sfx_expect = std::sqrt(0.25/3.);
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL( sfx_expect, sfx_res_turbo.getValue(), 1.e-8);
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_slow.getValue());
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_fast.getValue());

    /// the function to be integrated:  (Integral[ abs(x), dxdydz]) =?= (2 * |x|^2/2 @ [0, 0.5]) ==> .25)
    Norm<1> l1Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
    l1Norm(sfx, sfx_res);

    Norm<1> l1Norm_turbo(bulkData, &metaData.universal_part(), TURBO_NONE);
    l1Norm_turbo(sfx, sfx_res_turbo);

    sfx_expect = 0.25;
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());

    /// the function to be integrated:  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest1.py)
    StringFunction sfxyz("x*y*z", Name("sfxyz"), Dimensions(3), Dimensions(1) );

    TIME_IT2(l2Norm(sfxyz, sfx_res); ,  l2Norm_turbo(sfxyz, sfx_res_turbo); ,  "should be the same", turboOpt );

    sfx_expect = 0.0240562612162344;
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());

    /// the function to be integrated (but over a rotated domain):  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py)
    /// now rotate the mesh
    Math::Matrix rmz = Math::rotationMatrix(2, 30);
    Math::Matrix rm = rmz;
    MeshTransformer meshRotate(rm);
    eMesh.nodalOpLoop(meshRotate, coords_field);

    TIME_IT2( l2Norm(sfxyz, sfx_res); , l2Norm_turbo(sfxyz, sfx_res_turbo); , "should be the same", turboOpt );

    sfx_expect = 0.0178406008037016;
    // NOTE: we need extra quadrature accuracy to reproduce this result (cubDegree==4 in IntegratedOp almost gets it right)
    //   for now, we are satisfied with 3 digits
    //STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_res.getValue(), sfx_expect);
    if (std::fabs(sfx_res.getValue()-sfx_expect) > 0.01*sfx_expect)
    {
      STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
      STKUNIT_EXPECT_TRUE(false);
    }
    if (std::fabs(sfx_res_turbo.getValue()-sfx_expect) > 0.01*sfx_expect)
    {
      STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());
      STKUNIT_EXPECT_TRUE(false);
    }
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(norm, string_function_turbo_timings)
{
  EXCEPTWATCH;
  TEST_norm_string_function_turbo_timings(TURBO_ELEMENT);
}
STKUNIT_UNIT_TEST(norm, string_function_turbo_timings_bucket)
{
  EXCEPTWATCH;
  TEST_norm_string_function_turbo_timings(TURBO_BUCKET);
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(norm, field_function)
{
  EXCEPTWATCH;
  //stk::diag::WriterThrowSafe _write_throw_safe(dw());
  //dw().setPrintMask(dw_option_mask.parse(vm["dw"].as<std::string>().c_str()));
  //dw().setPrintMask(LOG_NORM+LOG_ALWAYS);

  dw().m(LOG_NORM) << "TEST.norm.field_function " << stk::diag::dendl;

  LocalFixture     fix(4);
  mesh::fem::FEMMetaData&        metaData     = fix.metaData;
  mesh::BulkData&        bulkData     = fix.bulkData;
  //PerceptMesh&        eMesh     = fix.eMesh;
  mesh::FieldBase*       coords_field = fix.coords_field;

  /// Create the operator that will do the work
  /// get the l2 norm
  Norm<2> l2Norm(bulkData, &metaData.universal_part(), TURBO_NONE);

  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  std::vector<double> vals(3,0.0);
  ConstantFunctionVec sfx_res_vec(vals, "sfx_res_vec");

  Util::setFlag(0, false);
  l2Norm(ff_coords, sfx_res_vec);

  double sfx_expect = std::sqrt(0.25/3.);
  //std::cout << "tmp vec" << sfx_res_vec.getValue() << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX(sfx_expect, sfx_res_vec.getValue()[0]);

  /// the function to be integrated:  (Integral[ abs(x), dxdydz]) =?= (2 * |x|^2/2 @ [0, 0.5]) ==> .25)
  Norm<1> l1Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
  l1Norm(ff_coords, sfx_res_vec);

  sfx_expect = 0.25;
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_vec.getValue()[0]);

  Util::setFlag(0, false);
}

}
}
}
