#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_PostProcess.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_String_Function_Expression.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const std::function<double(const stk::math::Vector3d &)> & analytic_fn)
{

  const stk::mesh::Selector activeOwnedFieldSelector = AuxMetaData::get(mesh.mesh_meta_data()).active_locally_owned_selector() & stk::mesh::selectField(distField);
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  double errorSum = 0.;
  double solutionSum = 0.;

  for ( auto && bucketPtr : mesh.get_buckets( stk::topology::NODE_RANK, activeOwnedFieldSelector) )
  {
    const size_t length = bucketPtr->size();
    const double *coordsData = field_data<double>(coordsField, *bucketPtr);
    double *dist = field_data<double>(distField, *bucketPtr);

    for (size_t i = 0; i < length; ++i)
    {
      const stk::math::Vector3d nodeCoords(coordsData+i*dim, dim);
      const double analytic = analytic_fn(nodeCoords);
      const double error = dist[i]-analytic;
      errorSum += error*error;
      solutionSum += analytic*analytic;
    }
  }
  all_reduce_sum(mesh.parallel(), errorSum);
  all_reduce_sum(mesh.parallel(), solutionSum);
  return std::sqrt(errorSum/solutionSum);
}

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const String_Function_Expression & analyticDist, const double time)
{
  auto analytic_fn = [&analyticDist, time](const stk::math::Vector3d &x) { return analyticDist.evaluate(time, x); };
  return compute_relative_nodal_RMS_error(mesh, coordsField, distField, analytic_fn);
}

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const String_Function_Expression & analyticDist)
{
  auto analytic_fn = [&analyticDist](const stk::math::Vector3d &x) { return analyticDist.evaluate(x); };
  return compute_relative_nodal_RMS_error(mesh, coordsField, distField, analytic_fn);
}

static void compute_and_print_distance_error(const stk::mesh::BulkData & mesh, const double time, const FieldRef coordsField, const FieldRef distField, const String_Function_Expression & analyticDist)
{
  krinolog << "Relative Nodal L2 error of " << distField.name() << " at time " << time << " = " << compute_relative_nodal_RMS_error(mesh, coordsField, distField, analyticDist, time) << stk::diag::dendl;
}

void PostProcessors::add_scalar_postprocesor(const std::string fieldName, const std::string & analyticalExpr)
{
  myScalarPostProcessorStrings.emplace_back(fieldName, analyticalExpr);
}

void PostProcessors::commit(const stk::mesh::MetaData & meta)
{
  for (auto & [fieldName, analyticalExpr] : myScalarPostProcessorStrings)
  {
    FieldRef field = AuxMetaData::get(meta).get_field(stk::topology::NODE_RANK, fieldName);
    STK_ThrowRequireMsg(field.valid(), "Cannot field " << fieldName << " for postprocessor.");
    myScalarPostProcessors.emplace_back(field, analyticalExpr);
  }
}

void PostProcessors::postprocess(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double time) const
{
  for (auto & scalarPostProcessor : myScalarPostProcessors)
    compute_and_print_distance_error(mesh, time, coordsField, scalarPostProcessor.first, scalarPostProcessor.second);
}

}
