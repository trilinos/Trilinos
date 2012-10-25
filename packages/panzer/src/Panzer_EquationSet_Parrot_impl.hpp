#ifndef __Panzer_EquationSet_Parrot_impl_hpp__
#define __Panzer_EquationSet_Parrot_impl_hpp__

namespace panzer {

template <typename EvalT>
EquationSet_Parrot<EvalT>::
EquationSet_Parrot(const EquationSetBase & eqSet,
                   const panzer::InputEquationSet& ies,
                   const panzer::CellData& cell_data,
                   const Teuchos::RCP<panzer::GlobalData>& global_data,
                   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(ies,cell_data,global_data,build_transient_support)
{
  // ********************
  // Assemble DOF names
  // ********************
  const std::vector<std::string> dofNames = eqSet.getDOFNames();
  for(std::size_t i=0;i<dofNames.size();i++) {
    this->m_dof_names->push_back(dofNames[i]);
    this->m_dof_gradient_names->push_back(dofNames[i]);
  }

  // this->m_residual_names 
  // this->m_scatter_name 

  // ********************
  // Build Basis Functions and Integration Rules
  // ********************
  
  this->setupDOFs(cell_data.baseCellDimension());
}

}

#endif
