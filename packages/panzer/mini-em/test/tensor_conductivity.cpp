// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Sacado_Fad_Kokkos.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

#include "MiniEM_TensorConductivity.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(tensor_conductivity,test2d,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(tensor_conductivity,test3d,TYPE)

namespace mini_em {

namespace {

// Levi-Civita symbol
int levi(const int i, const int j, const int k)
{ return (i-j)*(j-k)*(k-i)/2; }

/** Evaluates the basis coordinates field that TensorConductivity depends
  * on. TensorConductivity does not read it, so it is set to zero.
  */
template<typename EvalT, typename Traits>
class DummyCoordinates : public PHX::EvaluatorWithBaseImpl<Traits>,
                         public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  DummyCoordinates(const panzer::PureBasis& basis)
  {
    const std::string name = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis.name());
    coords = PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim>(name,basis.coordinates);
    this->addEvaluatedField(coords);
    this->setName("DummyCoordinates");
  }

  void evaluateFields(typename Traits::EvalData /* workset */)
  { Kokkos::deep_copy(coords.get_static_view(),0.0); }

private:
  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> coords;
};

/** Evaluates TensorConductivity on a quad (2D) or hex (3D) integration
  * rule with sigma and all three components of beta nonzero, and checks
  * every entry against the closed form
  *
  *   K_ij = sigma/(1+|beta|^2) (delta_ij + beta_i beta_j + eps_ikj beta_k),
  *
  * computed with the Levi-Civita symbol rather than the evaluator's
  * formula. In 2D only the in-plane block (i,j < 2) is checked, which is
  * what the 2D tensor should be. The beta_z terms there are the Hall terms
  * of the in-plane field.
  */
template<typename EvalType>
void testTensorConductivity(const int dim, Teuchos::FancyOStream& out, bool& success)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ScalarT = typename EvalType::ScalarT;

  const int numCells = 2;
  const double sigma = 3.0;
  const double beta[3] = {0.5, -0.7, 1.3};

  RCP<shards::CellTopology> topo = (dim == 2) ?
    rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())) :
    rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
  panzer::CellData cellData(numCells,topo);
  RCP<panzer::IntegrationRule> ir = rcp(new panzer::IntegrationRule(2,cellData));
  RCP<panzer::PureBasis> basis = rcp(new panzer::PureBasis("HCurl",1,cellData));
  panzer::FieldLayoutLibrary fl;
  fl.addFieldAndLayout("E_edge",rcp(new panzer::BasisIRLayout(basis,*ir)));

  panzer::Workset workset;
  workset.num_cells = numCells;
  workset.block_id = "eblock-0_0";

  PHX::FieldManager<panzer::Traits> fm;
  fm.registerEvaluator<EvalType>(rcp(new DummyCoordinates<EvalType,panzer::Traits>(*basis)));
  {
    RCP<PHX::Evaluator<panzer::Traits>> eval =
      rcp(new TensorConductivity<EvalType,panzer::Traits>("sigma",*ir,fl,sigma,beta[0],beta[1],beta[2],"E_edge"));
    fm.registerEvaluator<EvalType>(eval);
    for (const auto& ft : eval->evaluatedFields())
      fm.requireField<EvalType>(*ft);
  }

  panzer::Traits::SD setupData;
  {
    auto worksets = rcp(new std::vector<panzer::Workset>);
    worksets->push_back(workset);
    setupData.worksets_ = worksets;
  }
  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(4);
  fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);
  fm.postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;
  fm.preEvaluate<EvalType>(preEvalData);
  fm.evaluateFields<EvalType>(workset);
  fm.postEvaluate<EvalType>(0);

  PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim,panzer::Dim> K("sigma",ir->dl_tensor);
  fm.getFieldData<EvalType>(K);
  auto K_h = Kokkos::create_mirror_view(K.get_static_view());
  Kokkos::deep_copy(K_h, K.get_static_view());

  // K = sigma/(1+|beta|^2) (I + beta beta^T + [beta]_x), where
  // ([beta]_x)_ij = eps_ikj beta_k. In 2D the field lies in the x-y plane,
  // so K is the in-plane block of the 3D tensor with the same scaling.
  const double scale = sigma/(1.0 + beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2]);
  const double tol = 1e-14;
  for (int c = 0; c < numCells; ++c) {
    for (int p = 0; p < ir->num_points; ++p) {
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          double expected = ((i == j) ? 1.0 : 0.0) + beta[i]*beta[j];
          for (int k = 0; k < 3; ++k)
            expected += levi(i,k,j)*beta[k];
          expected *= scale;
          TEST_FLOATING_EQUALITY(Sacado::scalarValue(K_h(c,p,i,j)),expected,tol);
        }
      }
    }
  }
}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(tensor_conductivity,test2d,EvalType)
{
  testTensorConductivity<EvalType>(2,out,success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(tensor_conductivity,test3d,EvalType)
{
  testTensorConductivity<EvalType>(3,out,success);
}

typedef panzer::Traits::Residual ResidualType;
typedef panzer::Traits::Jacobian JacobianType;

UNIT_TEST_GROUP(ResidualType)
UNIT_TEST_GROUP(JacobianType)

}
