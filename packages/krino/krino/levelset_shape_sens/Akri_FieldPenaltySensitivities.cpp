#include <Akri_TypeDefs.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Sacado.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace
{
std::vector<double> get_element_nodal_scalar_values(const stk::mesh::BulkData & mesh, stk::mesh::Entity e, const stk::mesh::FieldBase & fld)
{
  unsigned nnodes = mesh.num_nodes(e);
  std::vector<double> result(nnodes);
  auto fldData = krino::FieldRef(fld).data<double>();
  for(unsigned n=0; n<mesh.num_nodes(e); n++)
  {
    result[n] = fldData.entity_values(mesh.begin_nodes(e)[n])(0_comp);
  }
  return result;
}

std::vector<Sacado::Fad::DFad<double>> 
get_element_nodal_scalar_values_with_sensitivities(const stk::mesh::BulkData & mesh, stk::mesh::Entity e, const stk::mesh::FieldBase & fld)
{
  unsigned nnodes = mesh.num_nodes(e);
  std::vector<Sacado::Fad::DFad<double>> result(nnodes);
  auto fldData = krino::FieldRef(fld).data<double>();
  for(unsigned n=0; n<mesh.num_nodes(e); n++)
  {
    result[n] = Sacado::Fad::DFad<double>(nnodes, n, fldData.entity_values(mesh.begin_nodes(e)[n])(0_comp));
  }
  return result;
}

class ElemShapeFcnOperations
{
public:
  ElemShapeFcnOperations(const stk::mesh::BulkData & mesh_, stk::mesh::Entity element, krino::FieldRef coords_field) :
    parametricElem(krino::AuxMetaData::get(mesh_.mesh_meta_data()).parametric_element_library().get_parametric_element(mesh_.bucket(element), coords_field)),
    mesh(mesh_), elem(element)
  {
    using namespace krino;
    const auto dim = mesh.mesh_meta_data().spatial_dimension();

    gradOp.resize(dim, parametricElem.num_nodes(), parametricElem.num_intg_pts());
    detJ.resize(parametricElem.num_intg_pts());
    double gradopError;

    unsigned nnodes = mesh.num_nodes(element);
    STK_ThrowRequire(nnodes == parametricElem.num_nodes());
    std::vector<double> coord(dim * nnodes);
    auto coordData = krino::FieldRef(*mesh.mesh_meta_data().coordinate_field()).data<double>();
    for(unsigned n=0; n<mesh.num_nodes(element); n++)
    {
      auto coordVals = coordData.entity_values(mesh.begin_nodes(element)[n]);
      for(unsigned d=0; d<dim; d++)
      {
        coord[n*dim + d] = coordVals(stk::mesh::ComponentIdx(d));
      }
    }

    parametricElem.gradient_operator(
      dim,                    // Number of coordinate dimensions
      parametricElem.num_intg_pts(),           // Number of target points
      parametricElem.num_nodes(),         // Number of coord shape functions
      parametricElem.shape_fcn_deriv(),  // Mesh shape function derivatives
      parametricElem.num_nodes(),          // Number of dof shape functions
      parametricElem.shape_fcn_deriv(),    // Dof shape function derivatives
      1,                      // Number of elements
      coord.data(),        // Mesh coordinate values
      gradOp.ptr(),          // Gradient operator values (output)
      detJ.ptr(),            // Determinant of the transformation Jacobian for each element (output)
      &gradopError);                // Gradop error (output)
  }

  template <typename DoubleType>
  std::vector<stk::math::Vec<DoubleType, 3>> compute_field_gradient_at_ips(const std::vector<DoubleType> & fldVals)
  {
    std::vector<stk::math::Vec<DoubleType, 3>> result(parametricElem.num_intg_pts(), stk::math::Vec<DoubleType, 3>::ZERO);
    for(unsigned ip=0; ip<parametricElem.num_intg_pts(); ip++)
    {
      for(unsigned d=0; d<parametricElem.topology_dimension(); d++)
      {
        for(unsigned n=0; n<parametricElem.num_nodes(); n++)
        {
          result[ip][d] += gradOp(d, n, ip) * fldVals[n];
        }
      }
    }

    return result;
  }

  template <typename DoubleType>
  std::vector<DoubleType> compute_field_at_ips(const std::vector<DoubleType> & fldVals)
  {
    std::vector<DoubleType> result(parametricElem.num_intg_pts(), DoubleType(0.));
    for(unsigned ip=0; ip<parametricElem.num_intg_pts(); ip++)
    {
      for(unsigned n=0; n<parametricElem.num_nodes(); n++)
      {
        result[ip] += parametricElem.shape_fcn()[ip * parametricElem.num_nodes() + n] * fldVals[n];
      }
    }
    return result;
  }

  template <typename DoubleType>
  DoubleType integrate_over_element(const std::vector<DoubleType> & ipVals)
  {
    DoubleType result(0.);
    for(unsigned ip=0; ip<parametricElem.num_intg_pts(); ip++)
    {
      result += ipVals[ip] * parametricElem.intg_weights()[ip] * detJ(ip);
    }
    return result;
  }
private:
  const krino::ParametricElement & parametricElem;
  const stk::mesh::BulkData & mesh;
  stk::mesh::Entity elem;
  sierra::ArrayContainer<double,krino::DIM,krino::NPE_COORD,krino::NINT> gradOp;
  sierra::ArrayContainer<double,krino::NINT> detJ;
};
}

namespace krino 
{

double compute_l2_penalty_between_fields(const stk::mesh::BulkData & mesh, const stk::mesh::FieldBase & fldWSens, 
  const stk::mesh::FieldBase & otherFld, double normalization, const stk::mesh::Selector & s, std::map<stk::mesh::EntityId, double> & sens)
{
  STK_ThrowRequire(mesh.is_automatic_aura_on());

  double penalty = 0;
  double vol = 0;
  sens.clear();

  for(auto && b : mesh.get_buckets(stk::topology::ELEM_RANK, s &
    (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().aura_part())))
  {
    const bool locallyOwned = b->owned();
    for(auto && e : (*b))
    {
      ElemShapeFcnOperations elemFcns(mesh, e, krino::FieldRef(mesh.mesh_meta_data().coordinate_field()));
      auto sensVals = get_element_nodal_scalar_values_with_sensitivities(mesh, e, fldWSens);
      auto otherVals = get_element_nodal_scalar_values(mesh, e, otherFld);

      auto ipSensVals = elemFcns.compute_field_at_ips(sensVals);
      auto ipOtherVals = elemFcns.compute_field_at_ips(otherVals);

      std::vector<Sacado::Fad::DFad<double>> penaltyIpVals(ipSensVals.size());
      for(unsigned ip=0; ip<ipSensVals.size(); ip++)
      {
        const Sacado::Fad::DFad<double> pen = (ipSensVals[ip]-ipOtherVals[ip])/normalization;
        penaltyIpVals[ip] = pen * pen;
      }
      auto elemPenalty = elemFcns.integrate_over_element(penaltyIpVals);

      if(locallyOwned) 
      {
        penalty += elemPenalty.val();
        vol += elemFcns.integrate_over_element(std::vector<double>(ipSensVals.size(), 1.));
      }

      for(unsigned n = 0; n < mesh.num_nodes(e); n++)
      {
        auto node = mesh.begin_nodes(e)[n];
        if(!mesh.bucket(node).owned()) continue;
        auto nId = mesh.identifier(node);
        auto sensId = sens.find(nId);
        if(sensId == sens.end()) sens[nId] = elemPenalty.dx(n);
        else sens[nId] += elemPenalty.dx(n);
      }
    }
  }

  double globalPenalty;
  double globalVol;
  stk::all_reduce_sum(mesh.parallel(), &penalty, &globalPenalty, 1);
  stk::all_reduce_sum(mesh.parallel(), &vol, &globalVol, 1);
  globalPenalty /= globalVol;
  for(auto && key : sens) sens[key.first] /= globalVol;
  return globalPenalty;
}

double compute_gradient_penalty_between_fields(const stk::mesh::BulkData & mesh, const stk::mesh::FieldBase & fldWSens, 
  const stk::mesh::FieldBase & otherFld, const stk::mesh::Selector & s, std::map<stk::mesh::EntityId, double> & sens)
{
  STK_ThrowRequire(mesh.is_automatic_aura_on());

  double penalty = 0;
  double vol = 0;
  sens.clear();

  for(auto && b : mesh.get_buckets(stk::topology::ELEM_RANK, s &
    (mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().aura_part())))
  {
    const bool locallyOwned = b->owned();
    for(auto && e : (*b))
    {
      ElemShapeFcnOperations elemFcns(mesh, e, krino::FieldRef(mesh.mesh_meta_data().coordinate_field()));
      auto sensVals = get_element_nodal_scalar_values_with_sensitivities(mesh, e, fldWSens);
      auto otherVals = get_element_nodal_scalar_values(mesh, e, otherFld);

      auto ipSensVals = elemFcns.compute_field_gradient_at_ips(sensVals);
      auto ipOtherVals = elemFcns.compute_field_gradient_at_ips(otherVals);

      std::vector<Sacado::Fad::DFad<double>> penaltyIpVals(ipSensVals.size());
      for(unsigned ip=0; ip<ipSensVals.size(); ip++)
      {
        penaltyIpVals[ip] = Sacado::Fad::DFad<double>(0.);
        for(unsigned d=0; d<mesh.mesh_meta_data().spatial_dimension(); d++)
        {
          const Sacado::Fad::DFad<double> diff = (ipSensVals[ip][d]-ipOtherVals[ip][d]);
          penaltyIpVals[ip] += diff * diff;
        }
      }
      auto elemPenalty = elemFcns.integrate_over_element(penaltyIpVals);

      if(locallyOwned) 
      {
        penalty += elemPenalty.val();
        vol += elemFcns.integrate_over_element(std::vector<double>(ipSensVals.size(), 1.));
      }

      for(unsigned n = 0; n < mesh.num_nodes(e); n++)
      {
        auto node = mesh.begin_nodes(e)[n];
        if(!mesh.bucket(node).owned()) continue;
        auto nId = mesh.identifier(node);
        auto sensId = sens.find(nId);
        if(sensId == sens.end()) sens[nId] = elemPenalty.dx(n);
        else sens[nId] += elemPenalty.dx(n);
      }
    }
  }

  double globalPenalty;
  double globalVol;
  stk::all_reduce_sum(mesh.parallel(), &penalty, &globalPenalty, 1);
  stk::all_reduce_sum(mesh.parallel(), &vol, &globalVol, 1);
  globalPenalty /= globalVol;
  for(auto && key : sens) sens[key.first] /= globalVol;
  return globalPenalty;
}
}