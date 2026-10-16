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
#include "PanzerDiscFE_config.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_CrossProduct.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(cross_product,test2d,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(cross_product,test3d,TYPE)

namespace panzer {

namespace {

// Input vectors at (cell, point, dim). They are not parallel at any point.
KOKKOS_INLINE_FUNCTION
double aValue(const int c, const int p, const int d)
{ return 1.0 + c + 0.1*p + 0.5*d; }

KOKKOS_INLINE_FUNCTION
double bValue(const int c, const int p, const int d)
{ return 2.0 - 0.3*c + 0.2*p*d + d*d; }

// Levi-Civita symbol
int levi(const int i, const int j, const int k)
{ return (i-j)*(j-k)*(k-i)/2; }

/** Evaluates the vectors A and B. For AD types, component d of A is seeded
  * with derivative index d and component d of B with derivative index
  * dim+d, so the derivatives of A x B are its Jacobian with respect to A
  * and B.
  */
template<typename EvalT, typename Traits>
class VectorPair : public PHX::EvaluatorWithBaseImpl<Traits>,
                   public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  VectorPair(const panzer::PointRule& pr)
  {
    a = PHX::MDField<ScalarT,Cell,Point,Dim>("A",pr.dl_vector);
    b = PHX::MDField<ScalarT,Cell,Point,Dim>("B",pr.dl_vector);
    this->addEvaluatedField(a);
    this->addEvaluatedField(b);
    this->setName("VectorPair");
  }

  void evaluateFields(typename Traits::EvalData workset)
  {
    auto av = a.get_static_view();
    auto bv = b.get_static_view();
    const int num_pts = av.extent_int(1);
    const int dim = av.extent_int(2);
    Kokkos::parallel_for("VectorPair",workset.num_cells,KOKKOS_LAMBDA (const int c) {
      for (int p = 0; p < num_pts; ++p) {
        for (int d = 0; d < dim; ++d) {
          av(c,p,d) = aValue(c,p,d);
          bv(c,p,d) = bValue(c,p,d);
          if constexpr (Sacado::IsADType<ScalarT>::value) {
            for (int k = 0; k < 2*dim; ++k) {
              av(c,p,d).fastAccessDx(k) = (k == d) ? 1.0 : 0.0;
              bv(c,p,d).fastAccessDx(k) = (k == dim+d) ? 1.0 : 0.0;
            }
          }
        }
      }
    });
  }

private:
  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,Cell,Point,Dim> a, b;
};

template<typename EvalType>
void testCrossProduct(const int dim, Teuchos::FancyOStream& out, bool& success)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ScalarT = typename EvalType::ScalarT;

  const int numCells = 2, numPoints = 3;

  RCP<shards::CellTopology> topo = (dim == 2) ?
    rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >())) :
    rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
  panzer::CellData cellData(numCells,topo);
  RCP<const panzer::PointRule> pr = rcp(new panzer::PointRule("points",numPoints,cellData));

  panzer::Workset workset;
  workset.num_cells = numCells;
  workset.block_id = "eblock-0_0";

  PHX::FieldManager<panzer::Traits> fm;
  fm.registerEvaluator<EvalType>(rcp(new VectorPair<EvalType,panzer::Traits>(*pr)));
  {
    Teuchos::ParameterList p;
    p.set("Result Name","AxB");
    p.set("Vector A Name","A");
    p.set("Vector B Name","B");
    p.set("Point Rule",pr);
    RCP<PHX::Evaluator<panzer::Traits>> eval = rcp(new panzer::CrossProduct<EvalType,panzer::Traits>(p));
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
  derivative_dimensions.push_back(2*dim);
  fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);
  fm.postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;
  fm.preEvaluate<EvalType>(preEvalData);
  fm.evaluateFields<EvalType>(workset);
  fm.postEvaluate<EvalType>(0);

  // In 2D the result is the scalar z-component, i = 2 below
  const int numComponents = (dim == 2) ? 1 : 3;
  PHX::MDField<ScalarT> axb("AxB", (dim == 2) ? pr->dl_scalar : pr->dl_vector);
  fm.getFieldData<EvalType>(axb);
  auto axb_v = axb.get_static_view();
  auto axb_h = Kokkos::create_mirror_view(axb_v);
  Kokkos::deep_copy(axb_h, axb_v);

  const double tol = 1e-14;
  for (int c = 0; c < numCells; ++c) {
    for (int p = 0; p < numPoints; ++p) {
      double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
      for (int d = 0; d < dim; ++d) {
        a[d] = aValue(c,p,d);
        b[d] = bValue(c,p,d);
      }

      for (int n = 0; n < numComponents; ++n) {
        const int i = (dim == 2) ? 2 : n;
        const ScalarT& result = (dim == 2) ? axb_h(c,p) : axb_h(c,p,n);

        // value: (a x b)_i = eps_ijk a_j b_k
        double expected = 0.0;
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k)
            expected += levi(i,j,k)*a[j]*b[k];
        TEST_FLOATING_EQUALITY(Sacado::scalarValue(result),expected,tol);

        // Jacobian: d(a x b)_i/da_j = eps_ijk b_k, d(a x b)_i/db_k = eps_ijk a_j
        if constexpr (Sacado::IsADType<ScalarT>::value) {
          TEST_EQUALITY(result.size(),2*dim);
          for (int j = 0; j < dim; ++j) {
            double expected_da = 0.0, expected_db = 0.0;
            for (int k = 0; k < 3; ++k) {
              expected_da += levi(i,j,k)*b[k];
              expected_db += levi(i,k,j)*a[k];
            }
            // absolute comparison since many entries are zero
            TEST_COMPARE(std::abs(result.dx(j)-expected_da),<,tol);
            TEST_COMPARE(std::abs(result.dx(dim+j)-expected_db),<,tol);
          }
        }
      }
    }
  }
}

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(cross_product,test2d,EvalType)
{
  testCrossProduct<EvalType>(2,out,success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(cross_product,test3d,EvalType)
{
  testCrossProduct<EvalType>(3,out,success);
}

typedef Traits::Residual ResidualType;
typedef Traits::Jacobian JacobianType;

UNIT_TEST_GROUP(ResidualType)
UNIT_TEST_GROUP(JacobianType)

}
