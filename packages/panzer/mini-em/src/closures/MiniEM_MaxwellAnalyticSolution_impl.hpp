#ifndef MINIEM_MAXWELLANALYTICSOLUTION_IMPL_HPP
#define MINIEM_MAXWELLANALYTICSOLUTION_IMPL_HPP

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

#include "Panzer_Traits.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_MathematicalConstants.hpp"
#include "Kokkos_MathematicalFunctions.hpp"


namespace mini_em {

  //**********************************************************************
  template <typename EvalT,typename Traits>
  MaxwellAnalyticSolution<EvalT,Traits>::MaxwellAnalyticSolution(const std::string & name,
                                                                 const panzer::IntegrationRule & ir,
                                                                 const panzer::FieldLayoutLibrary & fl,
                                                                 const double timeScale,
                                                                 const std::string& basisName)
  {
    using Teuchos::RCP;

    Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(basisName);

    Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
    ir_degree = ir.cubature_degree;
    ir_dim = ir.spatial_dimension;

    source = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);
    this->addEvaluatedField(source);

    std::string n = "Maxwell Analytic Solution";
    this->setName(n);

    timeScale_ = 1.0 / timeScale;

  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void MaxwellAnalyticSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                                    PHX::FieldManager<Traits>& /* fm */)
  {
    ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
  }

  //**********************************************************************
  template <typename EvalT,typename Traits>
  void MaxwellAnalyticSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
  {
    const double time = workset.time;
    const auto pi = Kokkos::numbers::pi_v<double>;

    const double timeScale = timeScale_;

    const auto coords = workset.int_rules[ir_index]->ip_coordinates.get_static_view();
    auto tmp_source = source.get_static_view();

    if (ir_dim == 3) {
      Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,source.extent_int(1)});
      Kokkos::parallel_for("panzer:MaxwellAnalyticSolution 3D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
          auto x = coords(cell,point,0);
          auto y = coords(cell,point,1);
          auto z = coords(cell,point,2);
          tmp_source(cell,point, 0) = pi * Kokkos::sin(timeScale*time) * Kokkos::cos(pi*x) * Kokkos::sin(pi*y) * Kokkos::sin(pi*z);
          tmp_source(cell,point, 1) = pi * Kokkos::sin(timeScale*time) * Kokkos::sin(pi*x) * Kokkos::cos(pi*y) * Kokkos::sin(pi*z);
          tmp_source(cell,point, 2) = pi * Kokkos::sin(timeScale*time) * Kokkos::sin(pi*x) * Kokkos::sin(pi*y) * Kokkos::cos(pi*z);
        });
    } else {
      Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,source.extent_int(1)});
      Kokkos::parallel_for("panzer:MaxwellAnalyticSolution 2D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
          auto x = coords(cell,point,0);
          auto y = coords(cell,point,1);
          tmp_source(cell,point, 0) = pi * Kokkos::sin(timeScale*time) * Kokkos::sin(pi*x) * Kokkos::sin(pi*y);
        });
    }
  }

  //**********************************************************************
}

#endif
