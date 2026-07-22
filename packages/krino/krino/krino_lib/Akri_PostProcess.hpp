#ifndef KRINO_KRINO_KRINO_LIB_AKRI_POSTPROCESS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_POSTPROCESS_HPP_

#include <list>
#include <stk_mesh/base/Types.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_String_Function_Expression.hpp>

namespace krino {

class PostProcessors
{
public:
  void add_scalar_postprocesor(const std::string fieldName, const std::string & analyticalExpr);
  void postprocess(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double time) const;
  void commit(const stk::mesh::MetaData & meta);
private:
  std::list<std::pair<std::string, std::string>> myScalarPostProcessorStrings;
  std::list<std::pair<FieldRef,String_Function_Expression>> myScalarPostProcessors;
};

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const std::function<double(const stk::math::Vector3d &)> & analytic_fn);

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const String_Function_Expression & analyticDist, const double time);

double compute_relative_nodal_RMS_error(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distField, const String_Function_Expression & analyticDist);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_POSTPROCESS_HPP_ */
