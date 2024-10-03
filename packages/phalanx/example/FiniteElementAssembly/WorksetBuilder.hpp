// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EXAMPLE_WORKSET_BUILDER
#define EXAMPLE_WORKSET_BUILDER

#include "Mesh.hpp"
#include "Workset.hpp"


struct WorksetBuilder {

  using team_t = Kokkos::TeamPolicy<PHX::exec_space>::member_type;
  struct CopyWorksetDetJac_Tag{};
  struct CopyWorksetGradBasisReal_Tag{};

  Kokkos::View<double**,PHX::Device> mesh_det_jac;
  Kokkos::View<double****,PHX::Device> mesh_grad_basis_real;
  Kokkos::View<double**,PHX::Device> workset_det_jac;
  Kokkos::View<double****,PHX::Device> workset_grad_basis_real;
  int first_cell_global_index;
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const CopyWorksetDetJac_Tag& , const team_t& team) const
  {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,mesh_det_jac.extent(1)), [=] (const int& qp) {
        workset_det_jac(cell,qp) = mesh_det_jac(cell+first_cell_global_index,qp);
        //printf("det_jac=%f\n",workset.det_jac_(cell,qp));
    });
  }
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const CopyWorksetGradBasisReal_Tag& , const team_t& team) const
  {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,mesh_det_jac.extent(1)), [=] (const int& qp) {
        for (int basis=0; basis < static_cast<int>(mesh_grad_basis_real.extent(2)); ++basis)
          for (int dim=0; dim < static_cast<int>(mesh_grad_basis_real.extent(3)); ++dim)
            workset_grad_basis_real(cell,qp,basis,dim) =
              mesh_grad_basis_real(cell+first_cell_global_index,qp,basis,dim);
    });
  }
  
  void buildWorksets(const int& workset_size,
                     const phx_example::Mesh& mesh,
                     std::vector<Workset>& worksets) {
    
    const int num_cells = static_cast<int>(mesh.getCoordinates().extent(0));
    const int num_worksets =
      num_cells / workset_size + ((num_cells % workset_size) == 0 ? 0 : 1);
    
    worksets.resize(num_worksets);
    
    for (int w=0; w < num_worksets; ++w) {
      
      int workset_num_cells = workset_size;
      // If on the last workset, it may not be full
      if ( (w == (num_worksets-1)) && ((num_cells % workset_size) != 0) ) {
        workset_num_cells = num_cells % workset_size;
      }
      
      worksets[w].num_cells_ = workset_num_cells;
      worksets[w].first_cell_global_index_ = w * workset_size;
      
      // need all gids
      worksets[w].gids_ = mesh.getGlobalIndices();
      
      // basis and weights are local, no need to copy data
      worksets[w].basis_ = mesh.getBasis();
      worksets[w].weights_ = mesh.getWeights();

      // bind objects for kokkos functors
      mesh_det_jac = mesh.getDetJac();
      mesh_grad_basis_real = mesh.getGradBasisReal();
      first_cell_global_index = worksets[w].first_cell_global_index_;
      
      worksets[w].det_jac_ = Kokkos::View<double**,PHX::Device>("workset_det_jac",
                                                                worksets[w].num_cells_,
                                                                mesh_det_jac.extent(1));
      
      worksets[w].grad_basis_real_ = Kokkos::View<double****,PHX::Device>("workset_grad_basis_real",
                                                                          worksets[w].num_cells_,
                                                                          mesh_grad_basis_real.extent(1),
                                                                          mesh_grad_basis_real.extent(2),
                                                                          mesh_grad_basis_real.extent(3));
      
      workset_det_jac = worksets[w].det_jac_;
      workset_grad_basis_real = worksets[w].grad_basis_real_;
      
      Kokkos::parallel_for(Kokkos::TeamPolicy<CopyWorksetDetJac_Tag,PHX::exec_space>(worksets[w].num_cells_,Kokkos::AUTO()), *this);
      Kokkos::parallel_for(Kokkos::TeamPolicy<CopyWorksetGradBasisReal_Tag,PHX::exec_space>(worksets[w].num_cells_,Kokkos::AUTO()), *this);
      Kokkos::fence();
    }
  }
  
};

#endif
