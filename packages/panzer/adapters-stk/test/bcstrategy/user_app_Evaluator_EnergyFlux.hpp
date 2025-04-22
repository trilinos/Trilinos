#ifndef USERAPP_EVALAUTOR_ENERGY_FLUX_HPP
#define USERAPP_EVALAUTOR_ENERGY_FLUX_HPP

#include "Phalanx_MDField.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Panzer_Evaluator_Macros.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace user_app {

  using panzer::Cell;
  using panzer::IP;
  using panzer::Dim;

  template <typename EvalT, typename TRAITS>
  class EnergyFluxEvaluator : public panzer::EvaluatorWithBaseImpl<TRAITS>,
                              public PHX::EvaluatorDerived<EvalT, TRAITS>  {
  public:
    using ScalarT = typename EvalT::ScalarT;

    EnergyFluxEvaluator(const std::string& flux_name,
                        const std::string& dof_name,
                        const std::string& basis_name,
                        const std::string& grad_dof_name,
                        const std::string& velocity_name,
                        const std::string& density_name,
                        const std::string& heat_capacity_name,
                        const std::string& thermal_conductivity_name,
                        const Teuchos::RCP<PHX::DataLayout>& vector_layout,
                        const Teuchos::RCP<PHX::DataLayout>& scalar_layout,
                        const bool apply_weak_dirichlet,
                        const double value)
      : flux_(flux_name,vector_layout),
        weak_dbc_simp_penalty_term_("WEAK_DBC_PENALTY_TERM_"+dof_name,scalar_layout),
        weak_dbc_simp_flux_term_("WEAK_DBC_FLUX_TERM_"+dof_name,vector_layout),
        dof_(dof_name,scalar_layout),
        grad_dof_(grad_dof_name,vector_layout),
        velocity_(velocity_name,vector_layout),
        density_(density_name,scalar_layout),
        heat_capacity_(heat_capacity_name,scalar_layout),
        thermal_conductivity_(thermal_conductivity_name,scalar_layout),
        side_normals_("Side Normal",vector_layout),
        apply_weak_dirichlet_(apply_weak_dirichlet),
        value_(value),
        basis_name_(basis_name)
    {
      this->addEvaluatedField(flux_);
      if (apply_weak_dirichlet) {
        this->addEvaluatedField(weak_dbc_simp_penalty_term_);
        this->addEvaluatedField(weak_dbc_simp_flux_term_);
        this->addDependentField(side_normals_);
      }
      this->addDependentField(dof_);
      this->addDependentField(grad_dof_);
      this->addDependentField(velocity_);
      this->addDependentField(density_);
      this->addDependentField(heat_capacity_);
      this->addDependentField(thermal_conductivity_);
      this->setName("user_app::EnergyFluxEvaluator - "+dof_name);
    }

    void postRegistrationSetup(typename TRAITS::SetupData sd,
                               PHX::FieldManager<TRAITS>& fm)
    {
      basis_index_ = getBasisIndex(basis_name_, (*sd.worksets_)[0], this->wda);
    }

    void evaluateFields(typename TRAITS::EvalData workset)
    {
      using team_policy = Kokkos::TeamPolicy<PHX::exec_space>::member_type;
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
      PHX::MDField<double,panzer::Cell,panzer::BASIS,panzer::IP,panzer::Dim> grad_basis = this->wda(workset).bases[basis_index_]->grad_basis;
      const int num_ip = grad_basis.extent(2);
      const int num_basis = grad_basis.extent(1);
      const int num_dim = grad_basis.extent(3);

      if (apply_weak_dirichlet_) {
        Kokkos::parallel_for("user_app::EnergyFluxEvalautor",policy,KOKKOS_CLASS_LAMBDA (const team_policy& team) {
          const int cell = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_ip), [&] (const int ip) {
            double one_over_h = 0.0;
            for (int bp = 0; bp < num_basis; ++bp) {
              double one_over_h2_bp = 0.0;
              for (int dim = 0; dim < num_dim; ++dim) {
                one_over_h2_bp += grad_basis(cell,bp,ip,dim) * grad_basis(cell,bp,ip,dim);
              }
              one_over_h += Kokkos::sqrt(one_over_h2_bp);
            }

            double penalty = 100.0;
            weak_dbc_simp_penalty_term_(cell,ip) = penalty * one_over_h * (dof_(cell,ip) - value_);
            
            for (int dim = 0; dim < flux_.extent_int(2); ++dim) {
              flux_(cell,ip,dim) = density_(cell,ip) * heat_capacity_(cell,ip) * value_ * velocity_(cell,ip,dim)
                - thermal_conductivity_(cell,ip) * grad_dof_(cell,ip,dim);

              weak_dbc_simp_flux_term_(cell,ip,dim) = thermal_conductivity_(cell,ip) * side_normals_(cell,ip,dim) * (dof_(cell,ip) - value_);
            }
          });
        });
      }
      else {
        Kokkos::parallel_for("user_app::EnergyFluxEvalautor",policy,KOKKOS_CLASS_LAMBDA (const team_policy& team) {
          const int cell = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,flux_.extent(1)), [&] (const int ip) {
            for (int dim = 0; dim < flux_.extent_int(2); ++dim) {
              flux_(cell,ip,dim) = density_(cell,ip) * heat_capacity_(cell,ip) * dof_(cell,ip) * velocity_(cell,ip,dim)
                - thermal_conductivity_(cell,ip) * grad_dof_(cell,ip,dim);
            }
          });
        });
      }
    }

  private:
    PHX::MDField<ScalarT,Cell,IP,Dim> flux_;
    PHX::MDField<ScalarT,Cell,IP> weak_dbc_simp_penalty_term_;
    PHX::MDField<ScalarT,Cell,IP,Dim>  weak_dbc_simp_flux_term_;
    PHX::MDField<const ScalarT,Cell,IP> dof_;
    PHX::MDField<const ScalarT,Cell,IP,Dim> grad_dof_;
    PHX::MDField<const ScalarT,Cell,IP,Dim> velocity_;
    PHX::MDField<const ScalarT,Cell,IP> density_;
    PHX::MDField<const ScalarT,Cell,IP> heat_capacity_;
    PHX::MDField<const ScalarT,Cell,IP> thermal_conductivity_;
    PHX::MDField<const ScalarT,Cell,IP,Dim> side_normals_;
    bool apply_weak_dirichlet_;
    double value_;
    size_t basis_index_;
    std::string basis_name_;
  };

}

#endif
