// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LevelSet.hpp>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_ClosestPointRedistance.hpp>
#include <Akri_Compute_Surface_Distance.hpp>
#include <Akri_ContourElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_TypeDefs.hpp>
#include <Akri_Facet.hpp>
#include <Akri_FacetsFromSides.hpp>
#include <Akri_IC_Alg.hpp>
#include <Akri_Fast_Marching.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <Akri_Phase_Support.hpp>

#include <math.h>
#include <map>
#include <iomanip>
#include <cstdio>

#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_FastIterativeMethod.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_OrientedSideNodes.hpp>
#include <Akri_Surface_Manager.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_PatchInterpolator.hpp>
#include <Akri_Sign.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

bool all_nodes_have_field_data(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity entity, const stk::mesh::FieldBase& field)
{
  const unsigned nnodes = stk_bulk.num_nodes(entity);
  const stk::mesh::Entity* nodes = stk_bulk.begin_nodes(entity);
  for (unsigned i = 0; i < nnodes; ++i )
  {
    if (field_bytes_per_entity(field, nodes[i]) == 0)
    {
      return false;
    }
  }
  return true;
}

stk::mesh::BulkData & LevelSet::mesh()
{
  return my_meta.mesh_bulk_data();
}
const stk::mesh::BulkData & LevelSet::mesh() const
{
  return my_meta.mesh_bulk_data();
}

stk::mesh::MetaData & LevelSet::meta()
{
  return my_meta;
}
const stk::mesh::MetaData & LevelSet::meta() const
{
  return my_meta;
}

AuxMetaData & LevelSet::aux_meta()
{
  return my_aux_meta;
}
const AuxMetaData & LevelSet::aux_meta() const
{
  return my_aux_meta;
}

bool LevelSet::has_IC_surfaces()
{
  IC_Alg& ic_alg = get_IC_alg();
  return ic_alg.numberSurfaces() > 0;
}

BoundingBox LevelSet::get_IC_surface_bounding_box()
{
  IC_Alg& ic_alg = get_IC_alg();
  return ic_alg.get_surface_bounding_box();
}

IC_Alg& LevelSet::get_IC_alg()
{
  if (!my_IC_alg)
  {
    my_IC_alg = std::make_unique<IC_Alg>(*this);
  }
  return *my_IC_alg;
}

void LevelSet::setup(stk::mesh::MetaData & meta)
{
  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
  {
    ls->setup();
  }
}

void LevelSet::post_commit_setup(stk::mesh::MetaData & meta)
{
  const double max_elem_size = compute_maximum_element_size(meta.mesh_bulk_data(), meta.universal_part());

  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
  {
    ls->setup();
    if (ls->my_narrow_band_multiplier > 0.)
    {
      ls->narrow_band_size(ls->my_narrow_band_multiplier * max_elem_size);
    }
    else if (ls->narrow_band_size() > 0.)
    {
      const double narrow_band = ls->narrow_band_size();
      STK_ThrowErrorMsgIf(!(narrow_band > max_elem_size),
          "Currently, narrow_band_size must be greater than the maximum element size of " << max_elem_size << std::endl
          << "in order to avoid unintentional accuracy degradation. If this feature is needed, please contact krino developers.");
    }
  }
}

void LevelSet::set_flag_is_transient(stk::mesh::MetaData & meta)
{
  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
    ls->set_flag_is_transient();
}

void LevelSet::setup()
{
  register_fields();

  facets = FacetedSurfaceBase::build(my_meta.spatial_dimension());

  // initializes CDFEM_Support
  if (!meta().is_commit())
  {
    CDFEM_Support::get(my_meta);
  }
}

void LevelSet::register_nodal_levelset_field_and_set_ref()
{
  const unsigned numStates = (myIsTransient || (my_redistance_method == FAST_ITERATIVE)) ? 2 : 1;

  const FieldType & type_double  = FieldType::REAL;

  if(krinolog.shouldPrint(LOG_DEBUG))
    krinolog << "KRINO: Registering levelset field with name '" << get_isovar_field_name() << "'." << stk::diag::dendl;

  FieldRef isovarRef = aux_meta().declare_field(get_isovar_field_name(), type_double, stk::topology::NODE_RANK, numStates);

  const bool cdfem_is_active = krino::CDFEM_Support::is_active(meta());
  if (cdfem_is_active)
  {
    Phase_Support & phase_support = Phase_Support::get(meta());
    for (auto * partPtr : meta().get_mesh_parts())
    {
      if (partPtr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
          phase_support.level_set_is_used_by_nonconformal_part(get_identifier(), &phase_support.find_nonconformal_part(*partPtr)))
      {
        stk::mesh::put_field_on_mesh(isovarRef.field(), *partPtr, 1, nullptr);
      }
    }
  }
  else
  {
    stk::mesh::put_field_on_mesh(isovarRef.field(), meta().universal_part(), 1, nullptr);
  }

  set_isovar_field(isovarRef);
}

//--------------------------------------------------------------------------------
void LevelSet::register_fields()
{
  if (aux_meta().has_field(stk::topology::NODE_RANK, get_isovar_field_name()))
  {
    set_isovar_field( aux_meta().get_field(stk::topology::NODE_RANK, get_isovar_field_name()) );
  }
  else if (aux_meta().has_field(stk::topology::ELEMENT_RANK, get_isovar_field_name()))
  {
    STK_ThrowRequireMsg(trackIsoSurface, "Krino currently only support element isovar fields for level set-based death.");
    set_isovar_field( aux_meta().get_field(stk::topology::ELEMENT_RANK, get_isovar_field_name()) );
  }
  else
  {
    if (aux_meta().using_fmwk())
    {
      // fmwk region -> assume field is owned by other region and should already have been registered
      if (trackIsoSurface)
      {
        ThrowRuntimeError("Level set-based death variable '" << get_isovar_field_name() << "' should already be registered.");
      }
      else
      {
        ThrowRuntimeError("ERROR: level set field " << get_isovar_field_name() << " is not registered for level set " << name() << ".  "
          << "This can be caused by an incorrect level set field.  "
          << "Or, in aria, there could be a conflict between the specified subindex and the named species.  "
          << "If so, try using a subindex greater than the number of species for your level set.");
      }
    }
    else
    {
      // Non fmwk region (ie krino) -> assume that we will own the level set field so register it now
      register_nodal_levelset_field_and_set_ref();
    }
  }

  if (!myTimeOfArrivalBlockSpeedsByName.empty())
  {
    myTimeOfArrivalBlockSpeeds.resize(my_meta.get_parts().size());
    for (auto entry : myTimeOfArrivalBlockSpeedsByName)
    {
      if (my_aux_meta.has_part(entry.first))
      {
        myTimeOfArrivalBlockSpeeds[my_aux_meta.get_part(entry.first).mesh_meta_data_ordinal()] = entry.second;
      }
      else
      {
        STK_ThrowErrorMsgIf(true, "Could not find block " << entry.first << " when setting speed for computing time-of-arrival.");
      }
    }
  }

  if (!my_time_of_arrival_element_speed_field_name.empty())
  {
    const bool hasSpeedField = aux_meta().has_field(stk::topology::ELEMENT_RANK, my_time_of_arrival_element_speed_field_name);
    STK_ThrowErrorMsgIf(!hasSpeedField, "Could not find element speed field " << my_time_of_arrival_element_speed_field_name << " for computing time-of-arrival.");
    myTimeOfArrivalElementSpeedField = aux_meta().get_field(stk::topology::ELEMENT_RANK, my_time_of_arrival_element_speed_field_name);
    STK_ThrowRequireMsg(myTimeOfArrivalBlockSpeeds.empty(), "Speed for time-of-arrival calculation should be specified via element speed or block speed (not both).");
  }

  const bool usingLocallyConservedRedistancing = true; // where should this be?
  if (usingLocallyConservedRedistancing)
  {
    const FieldType & type_double  = FieldType::REAL;
    myDistanceCorrectionNumerator = aux_meta().register_field( "DistanceCorrectionNumerator", type_double, stk::topology::NODE_RANK, 1, 1, meta().universal_part() );
    myDistanceCorrectionDenominator = aux_meta().register_field( "DistanceCorrectionDenominator", type_double, stk::topology::NODE_RANK, 1, 1, meta().universal_part() );
  }
}

//-----------------------------------------------------------------------------------

void
LevelSet::set_time_of_arrival_block_speed(const std::string & blockName, const double blockSpeed)
{
  std::string lowerBlockName = blockName;
  std::transform(lowerBlockName.begin(), lowerBlockName.end(), lowerBlockName.begin(), ::tolower);
  auto entry = myTimeOfArrivalBlockSpeedsByName.find(lowerBlockName);
  STK_ThrowRequireMsg(entry == myTimeOfArrivalBlockSpeedsByName.end(), "Speed for block " << blockName << " specified more than once.");
  myTimeOfArrivalBlockSpeedsByName[lowerBlockName] = blockSpeed;
}

//-----------------------------------------------------------------------------------
void
LevelSet::write_facets(void)
{
  /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::facets_exoii(void)"); /* %TRACE% */
  if (stk::is_true_on_all_procs(mesh().parallel(), (!facets || facets->size() == 0)))
    return;
  const std::string fileBaseName = "facets_" + name();
  if (2 == mesh().mesh_meta_data().spatial_dimension())
    krino::write_facets(facets->get_facets_2d(), fileBaseName, my_facetFileIndex++, mesh().parallel());
  else
    krino::write_facets(facets->get_facets_3d(), fileBaseName, my_facetFileIndex++, mesh().parallel());
}

//-----------------------------------------------------------------------------------
void
LevelSet::compute_surface_distance(const double narrowBandSize, const double farFieldValue)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::compute_surface_distance(void)"); /* %TRACE% */

  stk::mesh::Selector surface_selector = selectUnion(my_compute_surface_distance_parts);

  if (narrowBandSize > 0.0)
  {
    my_narrow_band_size = narrowBandSize;
  }

  const stk::mesh::Field<double>& coords = reinterpret_cast<const stk::mesh::Field<double>&>(get_coordinates_field().field());
  const stk::mesh::Field<double>& dist = reinterpret_cast<const stk::mesh::Field<double>&>(get_isovar_field().field());

  Compute_Surface_Distance::calculate(
     mesh(),
     get_timer(),
     coords,
     dist,
     surface_selector,
     my_narrow_band_size,
     farFieldValue);

  // output for time 0 is from old
  if (get_isovar_field().number_of_states() > 1)
  {
    stk::mesh::field_copy(get_isovar_field(), get_isovar_field().field_state(stk::mesh::StateOld));
  }

}

//-----------------------------------------------------------------------------------

void
LevelSet::advance_semilagrangian_using_velocity_string_expression(const double timeN, const double timeNp1)
{
  STK_ThrowRequireMsg(myInterfaceVelocity.size() == spatial_dimension, "Did not find interface velocity expression.  Was it provided?");
  const auto velAtClosestPt = build_extension_velocity_at_closest_point_using_string_expressions(myInterfaceVelocity);

  const stk::mesh::Selector activeFieldSelector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(get_isovar_field());
  const BoundingBox paddedNodeBBox = compute_padded_node_bounding_box_for_semilagrangian_using_string_velocity_expressions(mesh(), activeFieldSelector, timeN, timeNp1, get_coordinates_field(), myInterfaceVelocity, *facets);

  advance_semilagrangian(get_coordinates_field(), timeN, timeNp1, paddedNodeBBox, velAtClosestPt);
}

void
LevelSet::advance_semilagrangian(const FieldRef coordsField, const double timeN, const double timeNp1, const BoundingBox & paddedNodeBBox, const ExtensionVelocityFunction & velAtClosestPt)
{
  krinolog << "Evolving level set for " << name() << " using semilagrangian..." << stk::diag::dendl;
  stk::mesh::field_copy(get_isovar_field().field_state(stk::mesh::StateOld), get_isovar_field()); // 0th order predictor needed for preserving sign with narrow_band

  std::unique_ptr<FacetedSurfaceBase> oldFacets = FacetedSurfaceBase::build(my_meta.spatial_dimension());
  oldFacets->swap(*facets);

  auto oldExtV = build_extension_velocity_using_velocity_at_closest_point(*oldFacets, velAtClosestPt);

  const stk::mesh::Selector activeFieldSelector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(get_isovar_field());
  const double avgEdgeLength = compute_average_edge_length();

  oldFacets->prepare_to_compute(paddedNodeBBox, my_narrow_band_size);

  if (mySemiLagrangianAlg == NON_ADAPTIVE_SINGLE_STEP)
  {
    calc_single_step_nonadaptive_semilagrangian_nodal_distance_and_build_facets(mesh(), activeFieldSelector, timeN, timeNp1, coordsField, get_isovar_field(), oldExtV, my_narrow_band_size, avgEdgeLength, *oldFacets, *facets);
  }
  else if (mySemiLagrangianAlg == ADAPTIVE_SINGLE_STEP)
  {
    calc_single_step_semilagrangian_nodal_distance_and_build_facets(mesh(), activeFieldSelector, timeN, timeNp1, coordsField, get_isovar_field(), oldExtV, my_narrow_band_size, avgEdgeLength, *oldFacets, *facets);
  }
  else if (mySemiLagrangianAlg == ADAPTIVE_PREDICTOR_CORRECTOR)
  {
    std::unique_ptr<FacetedSurfaceBase> facetsPred = FacetedSurfaceBase::build(my_meta.spatial_dimension());
    predict_semilagrangian_nodal_distance_and_build_facets(mesh(), activeFieldSelector, timeN, timeNp1, coordsField, get_isovar_field(), oldExtV, my_narrow_band_size, avgEdgeLength, *oldFacets, *facetsPred);

    facetsPred->prepare_to_compute(paddedNodeBBox, my_narrow_band_size);
    auto predictExtV = build_extension_velocity_using_velocity_at_closest_point(*facetsPred, velAtClosestPt);
    correct_semilagrangian_nodal_distance_and_build_facets(mesh(), activeFieldSelector, timeN, timeNp1, coordsField, get_isovar_field(), oldExtV, predictExtV, my_narrow_band_size, avgEdgeLength, *oldFacets, *facets);
  }
  else
  {
    STK_ThrowRequireMsg(false, "Unrecognized Semi-Lagrangian algorithm " << mySemiLagrangianAlg);
  }

  if (krinolog.shouldPrint(LOG_FACETS))
  {
    write_facets();
  }
}

//-----------------------------------------------------------------------------------
void
LevelSet::initialize(stk::mesh::MetaData & meta)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::initialize(void)"); /* %TRACE% */

  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
    ls->initialize(0.);
}

//-----------------------------------------------------------------------------------
void
LevelSet::initialize(const double time)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::initialize(void)"); /* %TRACE% */

  if (trackIsoSurface)
  {
    return;
  }

  sync_all_fields_to_host();

  if (!my_compute_surface_distance_parts.empty())
  {
    compute_surface_distance();
    return;
  }

  krinolog << "Initializing levelset " << name() << "..." << stk::diag::dendl;

  /* process analytic surfaces */
  if (my_IC_alg)
  {
    my_IC_alg->execute(time);
  }

  if (compute_time_of_arrival())
  {
    fast_methods_redistance(my_meta.universal_part(), true);
  }
  else if (my_perform_initial_redistance)
  {
    constrained_redistance();
  }

  // Offset initialized LS if requested
  if (my_ic_offset != 0.0) increment_distance(my_ic_offset, false);

  // Scale initialized LS if requested
  if (my_ic_scale != 1.0) scale_distance(my_ic_scale);

  if (get_isovar_field().number_of_states() > 1)
    stk::mesh::field_copy(get_isovar_field(), get_isovar_field().field_state(stk::mesh::StateOld));
}

void
LevelSet::build_initial_facets(stk::mesh::MetaData & meta)
{
  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
    ls->build_initial_facets(0.);

  if (krinolog.shouldPrint(LOG_FACETS))
  {
    for (auto&& ls : surfaceManager.get_levelsets())
      ls->write_facets();
  }
}

bool LevelSet::can_create_adaptive_initial_facets_from_initial_surfaces_because_initial_distance_is_solely_from_initial_surfaces() const
{
  if (!my_compute_surface_distance_parts.empty() ||
      compute_time_of_arrival() ||
      my_perform_initial_redistance ||
      my_ic_offset != 0.0 ||
      my_ic_scale != 1.0 ||
      my_IC_alg->numberCalculators() > 0)
    return false;
  return true;
}

void LevelSet::build_initial_facets(const double time)
{
  if (get_isovar_field().entity_rank() != stk::topology::NODE_RANK)
  {
    krinolog << "Cannot build initial facets for invalid or non-nodal level set field " << get_isovar_field().name() << stk::diag::dendl;
    return;
  }

  const bool buildAdaptiveFacets = (mySemiLagrangianAlg == ADAPTIVE_SINGLE_STEP || mySemiLagrangianAlg == ADAPTIVE_PREDICTOR_CORRECTOR) &&
    can_create_adaptive_initial_facets_from_initial_surfaces_because_initial_distance_is_solely_from_initial_surfaces();

  Composite_Surface & initSurfaces = my_IC_alg->get_surfaces(); // Note that it is assumed that this fn is called right after initialize and therefore initSurfaces.prepare_to_compute has already been called.

  const stk::mesh::Selector activeFieldSelector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(get_isovar_field());

  const double avgEdgeLength = compute_average_edge_length();

  if (buildAdaptiveFacets)
    build_initial_adaptive_facets_after_nodal_distance_is_initialized_from_initial_surfaces(mesh(), activeFieldSelector, time, get_coordinates_field(), get_isovar_field(), avgEdgeLength, initSurfaces, *facets);
  else
    build_nonadaptive_facets(mesh(), activeFieldSelector, get_coordinates_field(), get_isovar_field(), avgEdgeLength, *facets);
}

//-----------------------------------------------------------------------------------
void
LevelSet::clear_initialization_data(stk::mesh::MetaData & meta)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::clear_initialization_data(stk::mesh::MetaData & meta)"); /* %TRACE% */

  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  for (auto&& ls : surfaceManager.get_levelsets())
    ls->clear_initialization_data();
}

void
LevelSet::clear_initialization_data()
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::clear_initialization_data(stk::mesh::MetaData & meta)"); /* %TRACE% */

  if (my_IC_alg && !get_keep_IC_surfaces())
    my_IC_alg->clear();
}

static void accumulate_levelset_integrals_on_element(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity elem,
    double & area,
    double & negVol,
    double & posVol,
    const double avgElemSize,
    const FieldRef coordsField,
    const FieldRef isoField,
    const double isoVal)
{
  ContourElement contourElement(mesh, elem, coordsField, isoField, isoVal);
  contourElement.compute_subelement_decomposition(avgElemSize);

  area += contourElement.compute_area_of_interface();
  negVol += contourElement.compute_signed_volume(-1);
  posVol += contourElement.compute_signed_volume(1);
}


void LevelSet::locally_conserved_redistance()
{
  const double avgElemSize = compute_average_edge_length();
  const double edgeTolForCalculatingCorrection = 1.e-6;

  sierra::ArrayContainer<double,DIM,NINT> intgPtLocations;
  sierra::ArrayContainer<double,NINT> intgWeights;
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double, NPE_VAR, NINT> weightFnsOnInterface;

  const FieldRef coordsField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();
  const FieldRef dField = get_isovar_field();

  stk::mesh::Selector active_field_selector =
      stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector<stk::mesh::Entity> elements;
  stk::mesh::get_selected_entities(
      active_field_selector, mesh().buckets(stk::topology::ELEMENT_RANK), elements);

  const double isoVal = 0.0;

  std::vector<double> negVolBefore;
  negVolBefore.reserve(elements.size());

  for (auto && elem : elements)
  {
    ContourElement contourElement(mesh(), elem, coordsField, isoField, isoVal);
    contourElement.compute_subelement_decomposition(avgElemSize, edgeTolForCalculatingCorrection);
    const double negVol = contourElement.compute_signed_volume(-1);
    negVolBefore.push_back(negVol);
  }

  redistance();

  const stk::mesh::Selector active_not_ghost_field_selector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dField);

  const int dim = mesh().mesh_meta_data().spatial_dimension();
  const double volNorm = std::pow(avgElemSize, dim);
  const double convergenceTol = 1.e-6;

  const bool useWeightedCorrection = true;
  const bool useWeightedErrorAndArea = !useWeightedCorrection;

  const int maxCorrectionSteps = 30;
  bool done = false;
  int iter = 0;
  while (!done)
  {
    stk::mesh::field_fill(0., myDistanceCorrectionNumerator);
    stk::mesh::field_fill(0., myDistanceCorrectionDenominator);

    double sumVolErrorSquared = 0.;
    double count = 0.;
    for (unsigned iElem = 0; iElem < elements.size(); iElem++)
    {
      const auto elem = elements[iElem];

      ContourElement contourElement(mesh(), elem, coordsField, isoField, isoVal);
      contourElement.compute_subelement_decomposition(avgElemSize, edgeTolForCalculatingCorrection);

      const double negVolCurrent = contourElement.compute_signed_volume(-1);

      const int numAreaIntgPts = contourElement.gather_intg_pts(0, intgPtLocations, intgWeights, determinants);
      const double area = ContourElement::compute_domain_integral(intgPtLocations, intgWeights, determinants);

      if (area == 0.)
        continue;

      weightFnsOnInterface.resize(contourElement.dist_topology().num_nodes(), numAreaIntgPts);
      contourElement.dist_master_elem().shape_fcn(numAreaIntgPts, intgPtLocations.ptr(), weightFnsOnInterface.ptr());

      const double volError = negVolCurrent - negVolBefore[iElem];

      const StkMeshEntities elemNodes{mesh().begin_nodes(elem), mesh().end_nodes(elem)};
      for (size_t i=0; i<elemNodes.size(); ++i)
      {
        double area_wi = 0;
        for (int ip = 0; ip < numAreaIntgPts; ++ip)
        {
          area_wi += weightFnsOnInterface(i, ip) * intgWeights(ip) * determinants(ip);
        }

        double & correctionNum = *(field_data<double>(myDistanceCorrectionNumerator, elemNodes[i]));
        double & correctionDenom = *(field_data<double>(myDistanceCorrectionDenominator, elemNodes[i]));

        if (useWeightedCorrection)
        {
          const double elementCorrection = volError/area;
          correctionNum += elementCorrection * area_wi;
          correctionDenom += area_wi;
        }
        else if (useWeightedErrorAndArea)
        {
          correctionNum += volError * area_wi;
          correctionDenom += area * area_wi;
        }
      }

      sumVolErrorSquared += (negVolCurrent - negVolBefore[iElem])*(negVolCurrent - negVolBefore[iElem]);
      count += 1.;
    }

    stk::mesh::parallel_sum(mesh(), {&myDistanceCorrectionNumerator.field(), &myDistanceCorrectionDenominator.field()});

    double sumCorrectionSquared = 0.;
    size_t countCorrection = 0;
    for ( auto && bucket : mesh().get_buckets( stk::topology::NODE_RANK, active_not_ghost_field_selector) )
    {
      const stk::mesh::Bucket & b = *bucket;

      double *d = field_data<double>(dField , b);
      double * correctionNum = field_data<double>(myDistanceCorrectionNumerator, b);
      double * correctionDenom = field_data<double>(myDistanceCorrectionDenominator, b);
      countCorrection += b.size();

      for (size_t i = 0; i < b.size(); ++i)
      {
        if (correctionDenom[i] != 0)
        {
          const double correction = correctionNum[i]/correctionDenom[i];
          d[i] += correction;
          sumCorrectionSquared += correction*correction;
        }
      }
    }  // end bucket loop

    all_reduce_sum(mesh().parallel(), sumCorrectionSquared);
    all_reduce_sum(mesh().parallel(), countCorrection);

    const double correctionNorm = std::sqrt(sumCorrectionSquared/countCorrection/avgElemSize);
    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "Iteration " << iter << ", error = " << std::sqrt(sumVolErrorSquared/count/volNorm) << ", correction = " << correctionNorm << "\n";
    done = !(++iter < maxCorrectionSteps) || (correctionNorm < convergenceTol);
  }
}

std::pair<double,double> LevelSet::get_conserved_negative_volume_and_time() const
{
  return std::make_pair(myConservedNegVolume, myConservedNegVolumeTime);
}

void LevelSet::set_conserved_negative_volume_and_time(const double vol, const double time)
{
  myConservedNegVolume = vol;
  myConservedNegVolumeTime = time;
}

double LevelSet::constrained_redistance(const bool use_initial_vol, const double & signChangePurturbationTol)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::constrained_redistance(const bool use_initial_vol)"); /* %TRACE% */

  // Steps:
  // 1. measure current volume
  // 2. perform regular redistance
  // 3. find offset needed to conserve volume
  // 4. increment nodal distance by this amount

  sync_all_fields_to_host();

  // measure current area and volumes
  double start_area, start_neg_vol, start_pos_vol;
  compute_sizes( start_area,start_neg_vol, start_pos_vol, 0. );
  if ( 0. == start_neg_vol || 0. == start_area )
  {
    krinolog << "Skipping redistancing operation.  Volume and/or area have zero value.  Area: "
        << start_area << " Volume: " << start_neg_vol << stk::diag::dendl;
    return 0.;
  }

  if (use_initial_vol)
  {
    krinolog << "Performing conserved redistancing..." << stk::diag::dendl;

    if (myConservedNegVolume <= 0.0)
    {
      myConservedNegVolume = start_neg_vol;
    }
    start_pos_vol = start_neg_vol + start_pos_vol - myConservedNegVolume;
    start_neg_vol = myConservedNegVolume;
  }
  else
  {
    krinolog << "Performing constrained redistancing..." << stk::diag::dendl;
  }

  // perform regular redistance
  redistance();

  krinolog << "Correcting for volume change:" << stk::diag::dendl;

  // find correction needed to conserve volume

  const double correction = find_redistance_correction( start_area, start_neg_vol, start_pos_vol );

  double area_i = 0;
  double neg_vol_i = 0;
  double pos_vol_i = 0;

  compute_sizes(area_i, neg_vol_i, pos_vol_i, 0);
  // update nodal distance field
  increment_distance( -correction, true, signChangePurturbationTol);

  compute_sizes(area_i, neg_vol_i, pos_vol_i, 0);

  return myConservedNegVolume;
}

//--------------------------------------------------------------------------------

static std::function<std::pair<double,double>(const double)> build_volume_error_function_with_derivative(LevelSet & ls, const double startingNegVol, const double startingPosVol)
{
  auto volume_error_function_with_derivative = [&ls, startingNegVol, startingPosVol](const double x)
    {
      double area, negVol, posVol;
      ls.compute_sizes( area, negVol, posVol, x );
      const double totVol = startingNegVol+startingPosVol;
      const double relativeError = (negVol - startingNegVol) / totVol;
      krinolog << "  Correction = " << x
            << ", Current volume = " << negVol
            << ", Target volume = " << startingNegVol
            << ", Relative Error = " << std::abs(relativeError)
            << stk::diag::dendl;
      const double derivative = area/totVol;
      return std::make_pair(relativeError, derivative);
    };
  return volume_error_function_with_derivative;
}

double LevelSet::find_redistance_correction(const double /*start_area*/,
    const double start_neg_vol,
    const double start_pos_vol,
    const int max_iterations,
    const double tol)
{ /* %TRACE% */ /* %TRACE% */
  auto volume_error_function_with_derivative =
      build_volume_error_function_with_derivative(*this, start_neg_vol, start_pos_vol);
  const auto result =
      find_root_newton_raphson(volume_error_function_with_derivative, 0., max_iterations, tol);

  if (!result.first)
  {
    stk::RuntimeWarningAdHoc() << "\nConstrained renormalization failed to converge to root within "
                               << max_iterations << " iterations. Continuing with correction "
                               << result.second << "\n";
  }
  return result.second;
}

void LevelSet::sync_all_fields_to_host()
{
  if (get_isovar_field().valid()) get_isovar_field().field().sync_to_host();
  if (get_time_of_arrival_element_speed_field().valid()) get_time_of_arrival_element_speed_field().field().sync_to_host();

  if (get_isovar_field().valid()) get_isovar_field().field().modify_on_host();
  if (get_time_of_arrival_element_speed_field().valid()) get_time_of_arrival_element_speed_field().field().modify_on_host();
}

void LevelSet::redistance_using_existing_facets(const stk::mesh::Selector & volumeSelector, const double elementLengthScale)
{
  // debugging
  if (krinolog.shouldPrint(LOG_FACETS))
    write_facets();

  const stk::mesh::Selector nodeSelector = volumeSelector & aux_meta().active_part();
  const ClosestPointRedistance closestPtRedist(mesh(), get_coordinates_field(), get_isovar_field(), my_redistance_timer);
  const bool doEnforceSignAwayFromInterface = true;
  closestPtRedist.redistance_using_facets(nodeSelector, *facets, my_narrow_band_size, doEnforceSignAwayFromInterface, elementLengthScale);
}

void LevelSet::redistance() { redistance(my_meta.universal_part()); }

void
LevelSet::redistance(const stk::mesh::Selector & volumeSelector)
{
  stk::diag::TimeBlock timer__(my_redistance_timer);

  STK_ThrowErrorMsgIf(!my_time_of_arrival_element_speed_field_name.empty(), "Redistancing a time-of-arrival field will corrupt it.");

  sync_all_fields_to_host();

  if (FAST_MARCHING == my_redistance_method || FAST_ITERATIVE == my_redistance_method)
  {
    fast_methods_redistance(volumeSelector);
    return;
  }
  STK_ThrowRequire(CLOSEST_POINT == my_redistance_method);

  krinolog << "Redistancing the level set field..." << stk::diag::dendl;

  const stk::mesh::Selector activeElementSelector = volumeSelector & aux_meta().active_part();
  const double avgEdgeLength = compute_average_edge_length();
  const ClosestPointRedistance closestPtRedist(mesh(), get_coordinates_field(), get_isovar_field(), my_redistance_timer);
  closestPtRedist.build_isosurface_facets(activeElementSelector, avgEdgeLength, *facets);
  redistance_using_existing_facets(volumeSelector, avgEdgeLength);
}

static std::vector<stk::mesh::Entity> get_owned_and_shared_interface_and_child_element_nodes(const stk::mesh::BulkData & mesh,
    const AuxMetaData & auxMeta,
    const CDFEM_Support & cdfemSupport,
    const stk::mesh::Selector interfaceSelector)
{
  std::vector<stk::mesh::Entity> initialNodes;
  const stk::mesh::Selector ownedOrSharedInterfaceSelector = interfaceSelector & auxMeta.active_part() & (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part());
  stk::mesh::get_selected_entities( ownedOrSharedInterfaceSelector, mesh.buckets(stk::topology::NODE_RANK), initialNodes, false );

  stk::mesh::Selector ownedChildSelector = cdfemSupport.get_child_part() & mesh.mesh_meta_data().locally_owned_part();
  for (auto * bucketPtr :  mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedChildSelector))
    for (auto elem : *bucketPtr)
      for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
        initialNodes.push_back(node);
  stk::util::sort_and_unique(initialNodes, stk::mesh::EntityLess(mesh));
  return initialNodes;
}

static void build_levelset_conforming_facets_with_interface_velocity(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const FieldRef interfaceVelocity,
    const unsigned numVelocityStates,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets)
{
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsIdentifier);
  const stk::mesh::Selector negativeSideBlockSelector = phaseSupport.get_negative_levelset_block_selector(lsIdentifier);

  build_interface_conforming_facets_with_interface_velocity(mesh, interfaceSelector, negativeSideBlockSelector, activePart, coordsField, interfaceVelocity, numVelocityStates, lsIdentifier, facets);
}

static void build_levelset_conforming_facets(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const stk::mesh::Part & activePart,
    const FieldRef coordsField,
    const Surface_Identifier lsIdentifier,
    FacetedSurfaceBase & facets)
{
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(lsIdentifier);
  const stk::mesh::Selector negativeSideBlockSelector = phaseSupport.get_negative_levelset_block_selector(lsIdentifier);

  build_interface_conforming_facets(mesh, interfaceSelector, negativeSideBlockSelector, activePart, coordsField, lsIdentifier, facets);
}

void
LevelSet::interface_conforming_redistance()
{
  if (FAST_MARCHING == my_redistance_method)
    krinolog << "Redistancing the level set field using CDFEM fast marching method..." << stk::diag::dendl;
  else
    krinolog << "Redistancing the level set field using CDFEM method..." << stk::diag::dendl;

  sync_all_fields_to_host();

  const auto & phaseSupport = Phase_Support::get(meta());
  build_levelset_conforming_facets(mesh(), phaseSupport, aux_meta().active_part(), get_coordinates_field(), my_identifier, *facets);

  if (FAST_MARCHING == my_redistance_method)
  {
    fast_marching_interface_conforming_redistance_using_existing_facets();
  }
  else
  {
    const double avgEdgeLength = compute_average_edge_length();
    redistance_using_existing_facets(meta().universal_part(), avgEdgeLength);
  }
}

void
LevelSet::fast_marching_interface_conforming_redistance_using_existing_facets()
{
  const auto & phaseSupport = Phase_Support::get(meta());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(my_identifier);
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta());
  const std::vector<stk::mesh::Entity> initialNodes = get_owned_and_shared_interface_and_child_element_nodes(mesh(), aux_meta(), cdfemSupport, interfaceSelector);

  const ClosestPointRedistance closestPtRedist(mesh(), get_coordinates_field(), get_isovar_field(), my_redistance_timer);
  const bool doEnforceSignAwayFromInterface = false;  // These nodes are just the ones along the interface and thus do not need to enforce the sign
  const double unusedElementLengthScale = 0.;
  closestPtRedist.redistance_given_nodes_using_facets(*facets, initialNodes, my_narrow_band_size, doEnforceSignAwayFromInterface, unusedElementLengthScale);

  const stk::mesh::Selector activeVolumeSelector = aux_meta().active_part();

  std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> get_interface_speed;
  Fast_Marching fm(mesh(),
        activeVolumeSelector,
        get_coordinates_field(),
        get_isovar_field(),
        get_interface_speed,
        my_redistance_timer);
  fm.redistance(initialNodes);
}

template <typename FACET>
void extend_velocity_to_selected_nodes(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & nodeSelector,
    const FieldRef coordsField,
    const FieldRef extendedVelocity,
    const Faceted_Surface<FACET> & interfaceFacets)
{
  for ( auto && bucketPtr : mesh.get_buckets( stk::topology::NODE_RANK, nodeSelector) )
  {
    const size_t length = bucketPtr->size();
    const double *coordsData = field_data<double>(coordsField, *bucketPtr);
    double *extVelData = field_data<double>(extendedVelocity, *bucketPtr);

    for (size_t i = 0; i < length; ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*FACET::DIM, FACET::DIM);
      const auto * nearest = interfaceFacets.get_closest_facet(nodeCoords);
      const stk::math::Vector3d nearestVel = nearest->velocity_at_closest_point(nodeCoords);
      for (int d=0; d<FACET::DIM; ++d)
        extVelData[i*FACET::DIM+d] = nearestVel[d];
    }
  }
}

void LevelSet::extend_interface_velocity_using_closest_point_projection(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef interfaceVelocity, const FieldRef extendedVelocity, const Surface_Identifier lsIdentifier)
{
  krinolog << "Extending velocity for " << Surface_Manager::get(mesh.mesh_meta_data()).get_name(lsIdentifier) << " using closest point projection..." << stk::diag::dendl;
  coordsField.field().sync_to_host();
  interfaceVelocity.field().sync_to_host();
  extendedVelocity.field().sync_to_host();

  const int nDim = mesh.mesh_meta_data().spatial_dimension();
  const unsigned numVelocityStates = 1; // I think we only use NP1?
  std::unique_ptr<FacetedSurfaceBase> interfaceFacets = FacetedSurfaceBase::build_with_velocity(nDim);
  const auto & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const auto & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  build_levelset_conforming_facets_with_interface_velocity(mesh, phaseSupport, auxMeta.active_part(), coordsField, interfaceVelocity, numVelocityStates, lsIdentifier, *interfaceFacets);

  const stk::mesh::Selector nodeSelector = auxMeta.active_not_ghost_selector() & stk::mesh::selectField(extendedVelocity);

  BoundingBox nodeBBox = krino::compute_nodal_bbox(mesh, nodeSelector, coordsField);
  interfaceFacets->prepare_to_compute(nodeBBox, 0.);

  if (3 == nDim)
    extend_velocity_to_selected_nodes(mesh, nodeSelector, coordsField, extendedVelocity, interfaceFacets->as_derived_type<FacetWithVelocity3d>());
  else
    extend_velocity_to_selected_nodes(mesh, nodeSelector, coordsField, extendedVelocity, interfaceFacets->as_derived_type<FacetWithVelocity2d>());
}

void LevelSet::advance_semilagrangian_using_interface_velocity(const FieldRef interfaceCoordsField, const FieldRef evaluationCoordsField, const double timeN, const double timeNp1, const FieldRef interfaceVelocity)
{
  interfaceCoordsField.field().sync_to_host();
  evaluationCoordsField.field().sync_to_host();
  interfaceVelocity.field().sync_to_host();

  const int nDim = mesh().mesh_meta_data().spatial_dimension();
  const unsigned numVelocityStates = (interfaceVelocity.number_of_states() == 1) ? 1 : 2; // If we have multiple states, average N & NP1
  std::unique_ptr<FacetedSurfaceBase> interfaceFacets = FacetedSurfaceBase::build_with_velocity(nDim);
  const auto & phaseSupport = Phase_Support::get(mesh().mesh_meta_data());
  build_levelset_conforming_facets_with_interface_velocity(mesh(), phaseSupport, aux_meta().active_part(), interfaceCoordsField, interfaceVelocity, numVelocityStates, get_identifier(), *interfaceFacets);

  mySemiLagrangianAlg = ADAPTIVE_SINGLE_STEP;

  auto velAtClosestPt = build_extension_velocity_at_closest_point_using_facets_with_velocity(nDim, *interfaceFacets);

  const stk::mesh::Selector activeFieldSelector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(get_isovar_field());
  const BoundingBox paddedNodeBBox = compute_padded_node_bounding_box_for_semilagrangian_using_facets_with_velocity(mesh(), activeFieldSelector, timeNp1-timeN, evaluationCoordsField, *interfaceFacets);

  interfaceFacets->prepare_to_compute(paddedNodeBBox, 0.);

  advance_semilagrangian(evaluationCoordsField, timeN, timeNp1, paddedNodeBBox, velAtClosestPt);
}

void LevelSet::set_interface_velocity( const std::vector<std::string> & interfaceVelocity )
{
  initialize_expression_vector(interfaceVelocity, myInterfaceVelocity);
}

double
LevelSet::get_time_of_arrival_speed(stk::mesh::Entity elem, ParallelErrorMessage& err) const
{
  double speed = 1.0;
  if (myTimeOfArrivalBlockSpeeds.empty())
  {
    if (myTimeOfArrivalElementSpeedField.valid())
    {
      const double * speedFieldData = field_data<double>(myTimeOfArrivalElementSpeedField, elem);
      if (!speedFieldData) err << "Missing element speed on element " << mesh().identifier(elem) << "\n";
      else speed = *speedFieldData;

      if (speed <= 0.0)
      {
        err << "Non positive-definite speed " << speed << " found on element " << mesh().identifier(elem) << "\n";
      }
    }
  }
  else
  {
    const stk::mesh::Part & elemPart = find_element_part(mesh(), elem);
    speed = myTimeOfArrivalBlockSpeeds[elemPart.mesh_meta_data_ordinal()];
    STK_ThrowAssert(speed >= 0.0); // Negative speeds should have already been caught and generated error.
    if (speed == 0.0)
      err << "Speed not specified for block " << elemPart.name() << "\n";
  }

  return speed;
}

void
LevelSet::fast_methods_redistance(const stk::mesh::Selector & volumeSelector, const bool compute_time_of_arrival)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::fast_marching_redistance(const stk::mesh::Selector & selector)"); /* %TRACE% */

  // Unlike redistance() this method provides an approximate (not exact) distance to the isosurface.
  // On the other hand, this method provide solve the Eikonal equation on the mesh.  This is slightly
  // different than a pure distance function because it provides the distance through domain.  It can't see
  // through walls like the redistance() method does.

  std::function<double(ParallelErrorMessage& err, stk::mesh::Entity)> get_interface_speed;
  if (compute_time_of_arrival)
    get_interface_speed = [&](ParallelErrorMessage& err, stk::mesh::Entity elem) { return get_time_of_arrival_speed(elem, err); };

  const stk::mesh::Selector activeVolumeSelector = aux_meta().active_part() & volumeSelector;

  if (my_redistance_method == FAST_MARCHING)
  {
    if (compute_time_of_arrival)
      krinolog << "Initializing the level set field to be the time-of-arrival using a fast marching method..." << stk::diag::dendl;
    else
      krinolog << "Redistancing the level set field using a fast marching method..." << stk::diag::dendl;

    Fast_Marching fm(mesh(),
        activeVolumeSelector,
        get_coordinates_field(),
        get_isovar_field(),
        get_interface_speed,
        my_redistance_timer);
    fm.redistance();
    if (get_isovar_field().number_of_states() > 1)
      stk::mesh::field_copy(get_isovar_field(), get_isovar_field().field_state(stk::mesh::StateOld));
  }
  else
  {
    STK_ThrowRequire(my_redistance_method == FAST_ITERATIVE);
    if (compute_time_of_arrival)
      krinolog << "Initializing the level set field to be the time-of-arrival using a fast iterative method..." << stk::diag::dendl;
    else
      krinolog << "Redistancing the level set field using a fast iterative method..." << stk::diag::dendl;

    FastIterativeMethod fim(mesh(),
        activeVolumeSelector,
        get_coordinates_field(),
        get_isovar_field(),
        get_interface_speed,
        my_redistance_timer);
    fim.redistance();
    STK_ThrowAssertMsg(fim.check_converged_solution(), "Fast iterative method did not fully converge.");
  }
}

//--------------------------------------------------------------------------------

void
LevelSet::set_distance( const double & in_distance ) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::set_distance( const double & in_distance ) const"); /* %TRACE% */

  stk::mesh::field_fill(in_distance, get_isovar_field());
}

void
LevelSet::increment_distance( const double increment, const bool enforce_sign, const double & signChangePurtubationTol )
{
  //
  // increment the distance everywhere by the given value
  //

  const FieldRef dField = get_isovar_field();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double *d = field_data<double>(dField , b);

    for (size_t i = 0; i < length; ++i)
    {
      const double change = (enforce_sign && sign_change(d[i],d[i]+increment)) ? -(1-signChangePurtubationTol)*d[i] : increment;
      d[i] += change;
    }
  }  // end bucket loop
}

void
LevelSet::scale_distance( const double scale) const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::scale_distance( const double & scale ) const"); /* %TRACE% */
  //
  // increment the distance everywhere by the given value
  //

  const FieldRef dField = get_isovar_field();

  const stk::mesh::Selector active_field_selector = aux_meta().active_not_ghost_selector() & stk::mesh::selectField(dField);
  stk::mesh::BucketVector const& buckets = mesh().get_buckets( stk::topology::NODE_RANK, active_field_selector);

  stk::mesh::BucketVector::const_iterator ib = buckets.begin();
  stk::mesh::BucketVector::const_iterator ib_end = buckets.end();

  // iterate nodes, by buckets, and set distance
  for ( ; ib != ib_end ; ++ib )
  {
    const stk::mesh::Bucket & b = **ib;

    const size_t length = b.size();

    double *d = field_data<double>(dField , b);

    for (size_t i = 0; i < length; ++i)
    {
      d[i] *= scale;
    }
  }  // end bucket loop
}

void
LevelSet::negate_distance() const
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::increment_distance( const double & increment ) const"); /* %TRACE% */

  const FieldRef dRef = get_isovar_field();
  stk::mesh::field_scale(-1.0, dRef);
}

//-----------------------------------------------------------------------------------

void
LevelSet::compute_continuous_gradient() const
{ /* %TRACE[ON]% */ /* %TRACE% */
  //
  // Compute a mass-lumped continuous distance gradient
  //

  const int dim = mesh().mesh_meta_data().spatial_dimension();

  const FieldRef contGradRef = aux_meta().get_field(stk::topology::NODE_RANK, "CONT_GRAD");
  const FieldRef nodeAreaRef = aux_meta().get_field(stk::topology::NODE_RANK, "NODE_AREA");


  // initialize
  stk::mesh::field_fill(0.0, contGradRef);
  stk::mesh::field_fill(0.0, nodeAreaRef);

  // intg_pt_locations and intg_weights are just wrappers for framework data
  // determinants and grad_distance are arrays containing data that are resized
  // only as needed
  sierra::Array<const double,DIM,NINT> intg_pt_locations;
  sierra::Array<const double,NINT> intg_weights;
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double,DIM,NINT> grad_distance;

  // **************************
  // Not-covered-in-nightly RWH
  // **************************

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(contGradRef) & aux_meta().active_not_ghost_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    // create element
    ContourElement ls_elem( mesh(), elem, xField, isoField );

    ls_elem.std_intg_pts( intg_pt_locations, intg_weights, determinants,
                          ls_elem.dist_master_elem() );
    ls_elem.compute_distance_gradient( intg_pt_locations, grad_distance );

    const unsigned num_intg_pts = intg_pt_locations.dimension<NINT>();
    const int npe = ls_elem.dist_topology().num_nodes();

    const double * shape_fcn_ptr = ls_elem.dist_master_elem().shape_fcn();
    const sierra::Array<const double,NPE_VAR,NINT> shape_fcn(shape_fcn_ptr,npe,num_intg_pts);

    const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);

    for ( int i = 0; i < npe; ++i )
    {
      double * cont_grad = field_data<double>(contGradRef, elem_nodes[i]);
      double * area_ptr = field_data<double>(nodeAreaRef, elem_nodes[i]);
      double & area = *area_ptr;

      for (unsigned ip = 0; ip < num_intg_pts; ++ip )
      {
        const double NdV = shape_fcn(i,ip) * intg_weights(ip) * determinants(ip);
        area += NdV;

        for ( int d = 0; d < dim; d++ )
        {
          cont_grad[d] += grad_distance(d,ip) * NdV;
          //krinolog << "grad_distance(" << d << "," << ip << ") = " << grad_distance(d,ip) << stk::diag::dendl;
        }
      }
    }
  }

  // iterate nodes, by buckets, and calculate average gradient

  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double * cont_grad = field_data<double>(contGradRef, b);
    double * area = field_data<double>(nodeAreaRef, b);

    for (size_t i = 0; i < length; ++i)
    {
      for ( int d = 0; d < dim; d++ )
      {
        cont_grad[dim*i+d] = cont_grad[dim*i+d] / area[i];
      }
    }
  }
}

//-----------------------------------------------------------------------------------

bool
LevelSet::remove_wall_features() const
{ /* %TRACE[ON]% */ Trace trace__("LevelSet::remove_wall_features()"); /* %TRACE% */

  std::vector<double> dist;
  std::vector<double*> coords;
  int small_feature_removed = false;

  const FieldRef dField = get_isovar_field();
  const FieldRef coordinates_field = get_coordinates_field();

  stk::mesh::Selector surface_selector = selectUnion(my_surface_parts);

  const stk::mesh::Selector active_field_selector = stk::mesh::selectField(dField) & aux_meta().active_locally_owned_selector() & surface_selector;
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(meta().side_rank(), active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const stk::topology dist_topology = MasterElementDeterminer::get_field_topology(b, dField);
    const int npe_dist = dist_topology.num_nodes();
    dist.resize( npe_dist );
    coords.resize( npe_dist );

    const size_t length = b.size();
    for (size_t i = 0; i < length; ++i)
    {
      stk::mesh::Entity side = b[i];

      const stk::mesh::Entity* side_nodes = mesh().begin_nodes(side);

      for ( int n = 0; n < npe_dist; ++n )
      {
        dist[n] = *field_data<double>(dField, side_nodes[n]);

        double *coord_data = field_data<double>(coordinates_field,side_nodes[n]);
        coords[n] = coord_data;
      }

      for ( int n = 0; n < npe_dist; ++n )
      {
        const stk::mesh::Entity *elems = mesh().begin_elements(side); //get the element
        const int num_elems = mesh().num_elements(side);

        for (int l = 0; l < num_elems; ++l )
        {
          const stk::mesh::Entity elem = elems[l];

          //make sure we have the right element (active) and has the distance function
          //at every point in the element (to take gradient)

          if(!aux_meta().active_locally_owned_selector()(mesh().bucket(elem))) continue;
          if(!elem_has_field_data(dField, elem)) continue;

          //const double elem_len = calc_elem_len_normal(elem, side, coordinates_field);

          if(std::fabs(dist[n]) > my_max_feature_size) continue;

          ContourElement ls_elem( mesh(), elem, coordinates_field, dField );
          const stk::math::Vector3d p_coords(1/3., 1/3., 1/3.);
          const stk::math::Vector3d grad_dist_vec = ls_elem.distance_gradient(p_coords);

          stk::math::Vector3d face_normal;

          //assume linear tet or tri elements!
          if(spatial_dimension == 2)
          {
            face_normal = {-coords[1][1]+coords[0][1], coords[1][0]-coords[0][0], 0}; //simple 90 deg rotation
            const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);

            int num_elem_nodes = mesh().num_nodes(elem);

            for (int j = 0; j < num_elem_nodes; ++j) //find node in elem not part of side part, confirm normal points into elem
            {
              if(elem_nodes[j] != side_nodes[0] && elem_nodes[j] != side_nodes[spatial_dimension-1])
              {
                double *coord = field_data<double>(coordinates_field,elem_nodes[j]);
                const stk::math::Vector3d vec_check(coord[0]-coords[0][0], coord[1]-coords[0][1], 0);

                if(Dot(face_normal, vec_check) < 0)
                {
                  face_normal *= -1;
                  break;
                }
              }
            }
          }
          else
          {
            const stk::math::Vector3d x1 (coords[1][0]-coords[0][0], coords[1][1]-coords[0][1],coords[1][2]-coords[0][2]);
            const stk::math::Vector3d x2 (coords[2][0]-coords[0][0], coords[2][1]-coords[0][1],coords[2][2]-coords[0][2]);
            face_normal = -1.0*Cross(x1,x2);
          }

          face_normal.unitize();
          const double eta = -dist[n]/(Dot(grad_dist_vec, face_normal));

          if(eta < my_max_feature_size && eta > 0)
          {
            *field_data<double>(dField, side_nodes[n]) = -1.0 * dist[n];
            small_feature_removed = true;
          }
        }
      }
    }
  }

  if (stk::is_true_on_any_proc(mesh().parallel(), small_feature_removed))
  {
    stk::mesh::communicate_field_data(mesh(), {&dField.field()});
    return true;
  }

  return false;
}
//--------------------------------------------------------------------------------

bool
LevelSet::simple_remove_wall_features() const
{ /* %TRACE[ON]% */ Trace trace__("LevelSet::remove_wall_features()"); /* %TRACE% */

  int small_feature_removed = false;
  const FieldRef dField = get_isovar_field();

  stk::mesh::Selector surface_selector = selectUnion(my_surface_parts);

  const stk::mesh::Selector active_field_selector = stk::mesh::selectField(dField) & aux_meta().active_part() & surface_selector;
  stk::mesh::BucketVector const& buckets = mesh().get_buckets(stk::topology::NODE_RANK, active_field_selector);

  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;
    const size_t length = b.size();
    for (size_t i = 0; i < length; ++i)
    {
      stk::mesh::Entity node = b[i];
      double & dist = *field_data<double>(dField, node);

      if(std::fabs(dist) < my_max_feature_size)
      {
        dist*=-1;
        small_feature_removed = true;
      }
    }
  }

  return small_feature_removed;
}
//--------------------------------------------------------------------------------

void
LevelSet::set_surface_parts_vector()
{
  Phase_Support & my_phase_support = Phase_Support::get(meta());

  const bool includeInterfaceParts = false;
  my_surface_parts =  my_phase_support.get_conformal_parts_of_rank(meta().side_rank(), includeInterfaceParts);

  krinolog << "Removing small features less than " << my_max_feature_size << " on surfaces: ";
  for (auto * part : my_surface_parts)
    krinolog << part->name() << " ";
  krinolog << stk::diag::dendl;
}

//--------------------------------------------------------------------------------

bool
LevelSet::elem_has_field_data(const FieldRef &myField, const stk::mesh::Entity &elem) const
{
  const stk::mesh::Entity* elem_nodes = mesh().begin_nodes(elem);
  int num_nodes = mesh().num_nodes(elem);

  for (int k = 0; k < num_nodes; ++k)
  {
    if(!has_field_data(myField, elem_nodes[k])) return false;
  }

 return true;
}

//--------------------------------------------------------------------------------
void
LevelSet::compute_levelset_sizes( double & area, double & negVol, double & posVol, const FieldRef isovar, const double isoval ) const
{ 
  area = 0.0;
  negVol = 0.0;
  posVol = 0.0;

  const double h_avg = compute_average_edge_length();

  const FieldRef xField = get_coordinates_field();

  stk::mesh::Selector active_field_selector =
      stk::mesh::selectField(isovar) & aux_meta().active_locally_owned_selector();
  std::vector<stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities(
      active_field_selector, mesh().buckets(stk::topology::ELEMENT_RANK), objs);

  for (auto && elem : objs)
    accumulate_levelset_integrals_on_element(mesh(), elem, area, negVol, posVol, h_avg, xField, isovar, isoval);

  all_reduce_sum(mesh().parallel(), area, negVol, posVol);
}
//--------------------------------------------------------------------------------
void
LevelSet::compute_sizes( double & area, double & negVol, double & posVol, const double isoval ) const
{
  compute_levelset_sizes(area, negVol, posVol, get_isovar_field(), isoval);
}

static double get_gradient_magnitude_at_ip(const sierra::ArrayContainer<double,DIM,NINT> & gradDist, const int ip)
{
  double mag2GradPhi = 0.;
  for ( unsigned dim = 0; dim < gradDist.dimension(0); ++dim )
    mag2GradPhi += gradDist(dim,ip) * gradDist(dim,ip);
  return sqrt(mag2GradPhi);
}

//--------------------------------------------------------------------------------
void accumulate_side_area_and_gradient_magnitude_error_integrals(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distField,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Entity side,
    sierra::ArrayContainer<double,NINT> & determinants,
    sierra::ArrayContainer<double,DIM,NINT> & gradDist,
    double & area,
    double & sumL2Error)
{
  sierra::Array<const double,DIM,NINT> intgPtLocations;
  sierra::Array<const double,NINT> intgWeights;

  const StkMeshEntities sideElems{mesh.begin_elements(side), mesh.end_elements(side)};
  const stk::mesh::ConnectivityOrdinal * sideOrdinals = mesh.begin_element_ordinals(side);
  for (size_t i=0; i<sideElems.size(); ++i)
  {
    stk::mesh::Entity sideElem = sideElems[i];
    if (elementSelector(mesh.bucket(sideElem)))
    {
      ContourElement contourElem( mesh, sideElem, coordsField, distField );

      const int numIntgPts = contourElem.std_side_intg_pts( sideOrdinals[i], intgPtLocations, intgWeights, determinants, contourElem.dist_master_elem() );

      contourElem.compute_distance_gradient( intgPtLocations, gradDist );

      for ( int ip = 0; ip < numIntgPts; ++ip )
      {
        const double magGradPhi = get_gradient_magnitude_at_ip(gradDist, ip);

        sumL2Error += (magGradPhi - 1.) * (magGradPhi - 1.) * intgWeights(ip) * determinants(ip);
        area += intgWeights(ip) * determinants(ip);
      }
    }
  }
}

//--------------------------------------------------------------------------------
double
LevelSet::CDFEM_gradient_magnitude_error()
{
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double,DIM,NINT> gradDist;

  const FieldRef coordsField = get_coordinates_field();
  const FieldRef distField = get_isovar_field();
  coordsField.field().sync_to_host();
  distField.field().sync_to_host();

  const auto & phaseSupport = Phase_Support::get(meta());
  const stk::mesh::Selector interfaceSelector = phaseSupport.get_negative_levelset_interface_selector(my_identifier);

  stk::mesh::Selector sideSelector = interfaceSelector & aux_meta().active_locally_owned_selector();
  stk::mesh::Selector elementSelector = aux_meta().active_part();

  double area = 0.;
  double sumL2 = 0.;

  for ( auto && bucket : mesh().get_buckets(meta().side_rank(), sideSelector) )
    for ( auto side : *bucket)
      accumulate_side_area_and_gradient_magnitude_error_integrals(mesh(), coordsField, distField, elementSelector, side, determinants, gradDist, area, sumL2);

  all_reduce_sum(mesh().parallel(), area, sumL2);

  double L2Norm = 0.;
  if (area > 0.)
    L2Norm = std::sqrt(sumL2/area);

  krinolog << "Gradient jump norm for " << name() << ": L2 = " << L2Norm << stk::diag::dendl;

  return L2Norm;
}
//--------------------------------------------------------------------------------
double
LevelSet::gradient_magnitude_error(void)
{ /* %TRACE[ON]% */ /* %TRACE% */

  double area = 0., sum_L2 = 0., errorLoo = 0.;

  sierra::ArrayContainer<double,DIM,NINT> intg_pt_locations;
  sierra::ArrayContainer<double,NINT> intg_weights;
  sierra::ArrayContainer<double,NINT> determinants;
  sierra::ArrayContainer<double,DIM,NINT> grad_dist;

  const double h_avg = compute_average_edge_length();

  const FieldRef xField = get_coordinates_field();
  const FieldRef isoField = get_isovar_field();
  xField.field().sync_to_host();
  isoField.field().sync_to_host();

  stk::mesh::Selector active_field_selector = stk::mesh::selectField(isoField) & aux_meta().active_locally_owned_selector();
  std::vector< stk::mesh::Entity> objs;
  stk::mesh::get_selected_entities( active_field_selector, mesh().buckets( stk::topology::ELEMENT_RANK ), objs );

  for ( auto && elem : objs )
  {
    // create element that is decomposed into subelements
    ContourElement ls_elem( mesh(), elem, xField, isoField );

    ls_elem.compute_subelement_decomposition(h_avg);

    // get integration point locations, weights, and determinants for surface
    int num_intg_pts = ls_elem.gather_intg_pts( 0,		// interface
                                                intg_pt_locations,// gauss pt locations
                                                intg_weights,	// integration wts at gauss pts
                                                determinants ); // determinant at gauss pts

    ls_elem.compute_distance_gradient( intg_pt_locations, grad_dist );
    for ( int ip = 0; ip < num_intg_pts; ++ip )
    {
      double mag_grad_phi = 0.;
      for ( unsigned dim = 0; dim < spatial_dimension; ++dim ) mag_grad_phi += grad_dist(dim,ip) * grad_dist(dim,ip);
      mag_grad_phi = sqrt(mag_grad_phi);

      sum_L2 += (mag_grad_phi - 1.) * (mag_grad_phi - 1.) * intg_weights(ip) * determinants(ip);
      area += intg_weights(ip) * determinants(ip);

      if ( fabs(mag_grad_phi - 1.) > errorLoo ) errorLoo = fabs(mag_grad_phi - 1.);
    }
  }

  all_reduce_sum(mesh().parallel(), area, sum_L2);
  all_reduce_max(mesh().parallel(), errorLoo);

  double errorL2 = 0.;
  if (area > 0.)
    errorL2 = std::sqrt(sum_L2/area);

  return errorL2;
}
//--------------------------------------------------------------------------------

double
LevelSet::compute_average_edge_length() const
{
  const stk::mesh::Selector activeFieldSelector = stk::mesh::selectField(get_isovar_field()) & aux_meta().active_locally_owned_selector();
  return compute_global_average_edge_length_for_selected_elements(mesh(), get_coordinates_field(), activeFieldSelector);
}

//--------------------------------------------------------------------------------

LevelSet &
LevelSet::build(
    stk::mesh::MetaData & in_meta,
    const std::string & ls_name,
    const stk::diag::Timer & parent_timer )
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::build(stk::mesh::MetaData & in_meta, const std::string & ls_name, stk::diag::Timer & parent_timer)"); /* %TRACE% */
  Surface_Manager & surfaceManager = Surface_Manager::get(in_meta);

  STK_ThrowRequire(!surfaceManager.has_levelset(ls_name));
  LevelSet * ls = new LevelSet(in_meta, ls_name, parent_timer);
  surfaceManager.add_levelset(ls);
  return *ls;
}

//--------------------------------------------------------------------------------
LevelSet::LevelSet(stk::mesh::MetaData & in_meta,
    const std::string & in_name,
    const stk::diag::Timer & parent_timer)
    : my_meta(in_meta),
      my_aux_meta(AuxMetaData::get(in_meta)),
      my_identifier(Surface_Manager::get(in_meta).get_identifier(in_name)),
      my_name(Surface_Manager::get(in_meta).get_name(my_identifier)),
      my_parent_timer(parent_timer),
      my_timer("LevelSet", parent_timer),
      my_redistance_timer("Redistance", my_timer),
      spatial_dimension(in_meta.spatial_dimension()),
      my_narrow_band_multiplier(0.0),
      my_narrow_band_size(0.0),
      my_max_feature_size(-1.0),
      my_ic_offset(0.0),
      my_ic_scale(1.0),
      my_perform_initial_redistance(false),
      my_keep_IC_surfaces(false),
      my_threshold(0.0),
      my_redistance_method(CLOSEST_POINT),
      trackIsoSurface(false),
      my_needs_reinitialize_every_step(false)
{ /* %TRACE[ON]% */ Trace trace__("krino::LevelSet::LevelSet(stk::mesh::MetaData & in_meta, const std::string & ls_name, stk::diag::Timer & parent_timer)"); /* %TRACE% */
  my_coordinates_field = my_aux_meta.get_current_coordinates();

  // default names for distance and velocity
  // line commands are available for overriding these names
  const std::string distName = "D_" + name();
  set_levelset_field_name(distName);
}

LevelSet::~LevelSet()
{
}

//--------------------------------------------------------------------------------

std::string
print_sizes(const LevelSet & ls)
{
  std::ostringstream out;

  out << "P" << ls.mesh().parallel_rank() << ": facets for level set '" << ls.name() << "' " << std::endl;
  out << ls.get_facets().print_sizes() << std::endl;

  return out.str();
}
//--------------------------------------------------------------------------------

} // namespace krino
