#include "stk_util/diag/Timer.hpp"
#include "gtest/gtest.h"
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_MeshSpecs.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Akri_Refinement.hpp>
#include <Akri_AuxMetaData.hpp>
#include "Akri_TransitionElementEdgeMarker.hpp"
#include <Akri_Unit_LogRedirecter.hpp>
#include <Akri_ContourElement.hpp>
#include <Akri_OutputUtils.hpp>

namespace krino
{

class ConstrainedRedistance : public StkMeshTriFixture
{
public:
  ConstrainedRedistance()
      : logfile(),
        my_ls(LevelSet::build(mMesh.mesh_meta_data(), "LS", sierra::Diag::sierraTimer())),
        myRefinement(mMesh.mesh_meta_data(), &this->get_aux_meta().active_part(), sierra::Diag::sierraTimer())
  {
    my_ls.setup(); // registers field and sets field refs on the object

    //for potential UMR refinement in test
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    stk::mesh::FieldBase & field =
        meta.declare_field<int>(stk::topology::ELEMENT_RANK, myElementMarkerFieldName, 1);
    myElementMarkerField = FieldRef(field);
    int init_val = 1;
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), 1, 1, &init_val);
  }

  void mark_nonparent_elements()
  {
    for (auto && bucket :
        mMesh.get_buckets(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part()))
    {
      const auto markerValue =
          myRefinement.is_parent(*bucket) ? Refinement::RefinementMarker::NOTHING : Refinement::RefinementMarker::REFINE;
      auto * elemMarker = field_data<int>(myElementMarkerField, *bucket);
      for (size_t i = 0; i < bucket->size(); ++i)
        elemMarker[i] = static_cast<int>(markerValue);
    }
  }

  void do_refinement(int num_refine)
  {
    const TransitionElementEdgeMarker edgeMarker(mMesh, myRefinement, myElementMarkerField.name());

    for (int i = 0; i < num_refine; i++)
    {
      mark_nonparent_elements();
      myRefinement.do_refinement(edgeMarker);
    }
  }

  void build_single_regular_tri()
  {
    RegularTri meshSpec;
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }

  void build_quad_split_4tri()
  {
    QuadSplit4Tri meshSpec;
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1});
  }

  void create_scaled_linear_x_dist_field(const double & mult)
  {
    const FieldRef dField = my_ls.get_distance_field();
    const FieldRef xField = my_ls.get_coordinates_field();

    const stk::mesh::Selector selector = stk::mesh::selectField(dField);
    stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

    for (auto && bucket : buckets)
    {
      for (auto && node : *bucket)
      {
        double * dist = field_data<double>(dField, node);
        double * x = field_data<double>(xField, node);
        stk::math::Vector3d coords(x, 2);

        *dist = mult * coords[0]; // function of x only
      }
    }
  }

  void create_scaled_circle_dist_field(const double & circle_radius, const double & mult)
  {
    const FieldRef dField = my_ls.get_distance_field();
    const FieldRef xField = my_ls.get_coordinates_field();

    const stk::mesh::Selector selector =
        stk::mesh::selectField(dField); //& this->get_aux_meta().active_locally_owned_selector();
    stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

    double x_c = 0;// -0.01;
    double y_c = 0;//-0.01;
    for (auto && bucket : buckets)
    {
      for (auto && node : *bucket)
      {
        double * dist = field_data<double>(dField, node);
        double * x = field_data<double>(xField, node);
        stk::math::Vector3d coords(x, 2);
        const double r = std::sqrt((coords[0]-x_c) * (coords[0]-x_c) + (coords[1]-y_c) * (coords[1]-y_c));

        *dist = mult*(r - circle_radius);
      }
    }
  }

    void create_scaled_sinusoidal_circle_dist_field(const double & circle_radius, const double amp=0.15, const double period = 8)
  {
    const FieldRef dField = my_ls.get_distance_field();
    const FieldRef xField = my_ls.get_coordinates_field();

    const stk::mesh::Selector selector =
        stk::mesh::selectField(dField); //& this->get_aux_meta().active_locally_owned_selector();
    stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

    double x_c = 0;// -0.01;
    double y_c = 0;//-0.01;
    for (auto && bucket : buckets)
    {
      for (auto && node : *bucket)
      {
        double * dist = field_data<double>(dField, node);
        double * x = field_data<double>(xField, node);
        stk::math::Vector3d coords(x, 2);
        const double r = std::sqrt((coords[0]-x_c) * (coords[0]-x_c) + (coords[1]-y_c) * (coords[1]-y_c));

        const auto theta = std::atan2(coords[1],coords[0]);
        const double pr = amp*circle_radius*sin(period*theta);
        *dist = r - (circle_radius + pr);
      }
    }
  }

  std::vector<double> calc_vol_per_elem()
  {
    const auto isoField = my_ls.get_distance_field();
    const auto coordsField = my_ls.get_coordinates_field();

    stk::mesh::Selector active_field_selector =
        stk::mesh::selectField(isoField) & this->get_aux_meta().active_locally_owned_selector();
    std::vector<stk::mesh::Entity> elements;
    stk::mesh::get_selected_entities(active_field_selector, mMesh.buckets(stk::topology::ELEMENT_RANK), elements);

    std::vector<double> vols;
    for (auto && elem : elements)
    {
      ContourElement contourElement(mMesh, elem, coordsField, isoField, 0.);
      vols.push_back(contourElement.volume());
    }
    return vols;
  }

  std::vector<double> calc_neg_vol_per_elem()
  {
    const auto isoField = my_ls.get_distance_field();
    const auto coordsField = my_ls.get_coordinates_field();

    stk::mesh::Selector active_field_selector =
        stk::mesh::selectField(isoField) & this->get_aux_meta().active_locally_owned_selector();
    std::vector<stk::mesh::Entity> elements;
    stk::mesh::get_selected_entities(active_field_selector, mMesh.buckets(stk::topology::ELEMENT_RANK), elements);

    const double avgElemSize = my_ls.compute_average_edge_length();

    std::vector<double> elem_neg_vols;
    for (auto && elem : elements)
    {
      ContourElement contourElement(mMesh, elem, coordsField, isoField, 0.);
      contourElement.compute_subelement_decomposition(avgElemSize);
      const double negVol = contourElement.compute_signed_volume(-1);
      elem_neg_vols.push_back(negVol);
    }
    return elem_neg_vols;
  }

  double calc_negative_volume()
  {
    const auto isoField = my_ls.get_distance_field();
    const auto coordsField = my_ls.get_coordinates_field();

    stk::mesh::Selector active_field_selector =
        stk::mesh::selectField(isoField) & this->get_aux_meta().active_locally_owned_selector();
    std::vector<stk::mesh::Entity> elements;
    stk::mesh::get_selected_entities(active_field_selector, mMesh.buckets(stk::topology::ELEMENT_RANK), elements);

    const double avgElemSize = my_ls.compute_average_edge_length();

    double negVol = 0.;
    for (auto && elem : elements)
    {
      ContourElement contourElement(mMesh, elem, coordsField, isoField, 0.);
      contourElement.compute_subelement_decomposition(avgElemSize);
      negVol += contourElement.compute_signed_volume(-1);
    }

    const double localNegVol = negVol;
    stk::all_reduce_sum(mMesh.parallel(), &localNegVol, &negVol, 1);

    return negVol;
  }

  void test_local_volume_conservation(const std::vector<double> & beforeVol, const std::vector<double> & afterVol, const std::vector<double> & totalVol, const double goldTol)
  {
    ASSERT_EQ(beforeVol.size(), afterVol.size());
    ASSERT_EQ(beforeVol.size(), totalVol.size());

    for (int unsigned i = 0; i < beforeVol.size(); i++)
    {
      if(beforeVol[i] == 0 && afterVol[i] == 0) continue;

      const double relError = abs(afterVol[i] - beforeVol[i])/totalVol[i];

      EXPECT_LE(relError, goldTol) << "VOLUME CHANGED FOR i = " << i << "  before_vol = " << beforeVol[i]
          << "  after_vol = " << afterVol[i] << " relative error = " <<  relError << std::endl;
    }
  }

protected:
  LogRedirecter logfile;
  LevelSet & my_ls;
  Refinement myRefinement;
  std::string myElementMarkerFieldName{"ELEMENT_MARKER"};
  FieldRef myElementMarkerField;

};

namespace
{

} // namespace

TEST_F(ConstrainedRedistance, CreateDistanceField)
{
  build_single_regular_tri();
  create_scaled_linear_x_dist_field(1.0);
}

TEST_F(ConstrainedRedistance, DoRedistanceDistancedField)
{
  build_single_regular_tri();

  create_scaled_linear_x_dist_field(1.0);

  my_ls.constrained_redistance(true);

  const FieldRef dField = my_ls.get_distance_field();
  const FieldRef xField = my_ls.get_coordinates_field();

  const stk::mesh::Selector selector = stk::mesh::selectField(dField);
  stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

  for (auto && bucket : buckets)
  {
    for (auto && node : *bucket)
    {
      double * dist = field_data<double>(dField, node);
      double * x = field_data<double>(xField, node);
      EXPECT_NEAR(*dist, *x, 1e-8);
    }
  }
}

TEST_F(ConstrainedRedistance, DoRedistanceDistancedFieldRefine)
{
  build_single_regular_tri();
  do_refinement(3);

  create_scaled_linear_x_dist_field(1.0);

  my_ls.constrained_redistance(true);

  const FieldRef dField = my_ls.get_distance_field();
  const FieldRef xField = my_ls.get_coordinates_field();

  const stk::mesh::Selector selector = stk::mesh::selectField(dField);
  stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

  for (auto && bucket : buckets)
  {
    for (auto && node : *bucket)
    {
      double * dist = field_data<double>(dField, node);
      double * x = field_data<double>(xField, node);
      EXPECT_NEAR(*dist, *x, 1e-8);
    }
  }
}


TEST_F(ConstrainedRedistance, DoRedistanceNonDistanceField)
{
  build_single_regular_tri();
  create_scaled_linear_x_dist_field(2.0);

  my_ls.constrained_redistance(true);

  const FieldRef dField = my_ls.get_distance_field();
  const FieldRef xField = my_ls.get_coordinates_field();

  const stk::mesh::Selector selector = stk::mesh::selectField(dField);
  stk::mesh::BucketVector const & buckets = mMesh.get_buckets(stk::topology::NODE_RANK, selector);

  for (auto && bucket : buckets)
  {
    for (auto && node : *bucket)
    {
      double * dist = field_data<double>(dField, node);
      double * x = field_data<double>(xField, node);
      EXPECT_NEAR(*dist, *x, 1e-8);
    }
  }
}

TEST_F(ConstrainedRedistance, DoRedistanceCircleDistField)
{
  const double radius = 0.45;
  build_quad_split_4tri();

  do_refinement(3);

  create_scaled_circle_dist_field(radius, 5.);

  const double negVolInitial = calc_negative_volume();
  
  EXPECT_NEAR(negVolInitial, M_PI*radius*radius,1e-2);

  my_ls.redistance();

  const double negVolFinal = calc_negative_volume();

  EXPECT_NEAR(negVolInitial, 0.630540, 1e-6);
  EXPECT_NEAR(negVolFinal, 0.628500, 1e-6);
}


TEST_F(ConstrainedRedistance, DoConstrainedRedistanceCircleDistField)
{
  const double radius = 0.45;
  build_quad_split_4tri();

  do_refinement(3);

  create_scaled_circle_dist_field(radius, 5.);
  
  const double negVolInitial = calc_negative_volume();

  EXPECT_NEAR(negVolInitial, M_PI*radius*radius, 1e-2);

  my_ls.constrained_redistance(true, 1e-3);

  const double negVolFinal = calc_negative_volume();

  EXPECT_NEAR(negVolInitial, negVolFinal, 1e-5);
}

TEST_F(ConstrainedRedistance, DoLocalConstrainedRedistanceCircleDistField)
{
  const double radius = 0.45;
  build_quad_split_4tri();

  do_refinement(3);

  create_scaled_circle_dist_field(radius, 5.);

  const double negVolInitial = calc_negative_volume();
  const auto beforeVol = calc_neg_vol_per_elem();
  const auto totalVol = calc_vol_per_elem();

  my_ls.locally_conserved_redistance();

  const double negVolFinal = calc_negative_volume();
  const auto afterVol = calc_neg_vol_per_elem();

  EXPECT_NEAR(negVolInitial, negVolFinal, 1e-5) << "initial vol = " << negVolInitial << " final_vol = " << negVolFinal << std::endl;
  test_local_volume_conservation(beforeVol, afterVol, totalVol, 1.e-3);
}

TEST_F(ConstrainedRedistance, DoRedistanceSinusoidalCircleDistField)
{
  const bool doWriteOutput = false;

  const double radius = 0.25;
  build_quad_split_4tri();

  do_refinement(4);

  create_scaled_sinusoidal_circle_dist_field(radius);

  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output_dist.e", 1, 0.0);
  
  const double negVolInitial = calc_negative_volume();

  EXPECT_NEAR(negVolInitial, M_PI*radius*radius,1e-2);

  my_ls.redistance();

  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output_dist.e", 2, 1.0, stk::io::APPEND_RESULTS);

  const double negVolFinal = calc_negative_volume();

  EXPECT_NEAR(negVolInitial, 0.197918, 1e-6);
  EXPECT_NEAR(negVolFinal, 0.196678, 1e-6);
}

TEST_F(ConstrainedRedistance, DoGlobalConstrainedRedistanceSinusoidalCircleDistField)
{
  const bool doWriteOutput = false;

  const double radius = 0.25;
  build_quad_split_4tri();

  do_refinement(4);

  create_scaled_sinusoidal_circle_dist_field(radius);

  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output_global.e", 1, 0.0);
  
  const double negVolInitial = calc_negative_volume();

  my_ls.constrained_redistance();

  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output_global.e", 2, 1.0, stk::io::APPEND_RESULTS);

  const double negVolFinal = calc_negative_volume();

  const double goldGlobalTol = 1.e-3; // Surprisingly poor tolerance because sign change constraint limits correction.  This may be worth exploring further.
  EXPECT_NEAR(negVolInitial, negVolFinal, goldGlobalTol) << "initial vol = " << negVolInitial << " final_vol = " << negVolFinal << std::endl;
  std::cout << logfile.get_log() << std::endl;
}

TEST_F(ConstrainedRedistance, DoLocalConstrainedRedistanceSinusoidalCircleDistField)
{
  const bool doWriteOutput = false;

  const double radius = 0.25;
  build_quad_split_4tri();

  do_refinement(4);

  create_scaled_sinusoidal_circle_dist_field(radius);

  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output_original.e", 1, 0.0);
  if (doWriteOutput)
    output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output.e", 1, 0.0);
  
  const auto beforeVol = calc_neg_vol_per_elem();
  const auto totalVol = calc_vol_per_elem();
  const double negVolInitial = calc_negative_volume();

  const int maxRedistanceSteps = 1;
  for (int iter=0; iter<maxRedistanceSteps; ++iter)
  {
    my_ls.locally_conserved_redistance();

    if (doWriteOutput)
      output_composed_mesh_with_fields(mMesh, get_aux_meta().active_part(), "output.e", 2+iter, 1.0*(iter+1), stk::io::APPEND_RESULTS);
  }

  const auto afterVol = calc_neg_vol_per_elem();
  const double negVolFinal = calc_negative_volume();

  const double goldGlobalTol = 1.e-5;
  EXPECT_NEAR(negVolInitial, negVolFinal, goldGlobalTol) << "initial vol = " << negVolInitial << " final_vol = " << negVolFinal << std::endl;

  const double goldLocalTol = 2.5e-3; // Is this all we can expect?
  test_local_volume_conservation(beforeVol, afterVol, totalVol, goldLocalTol);

}

} // namespace krino
