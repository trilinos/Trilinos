#ifndef __Panzer_EquationSet_Parrot_impl_hpp__
#define __Panzer_EquationSet_Parrot_impl_hpp__

namespace panzer {

template <typename EvalT>
EquationSet_Parrot<EvalT>::
EquationSet_Parrot(const Teuchos::RCP<Teuchos::ParameterList>& plist,
		   const int& default_integration_order,
		   const EquationSetBase & eqSet,
                   const panzer::CellData& cell_data,
                   const Teuchos::RCP<panzer::GlobalData>& global_data,
                   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(plist,default_integration_order,cell_data,global_data,build_transient_support)
{
  // ********************
  // Assemble DOF names
  // ********************
  const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > > & dofs = eqSet.getProvidedDOFs();

  //const std::vector<std::string> dofNames = eqSet.getDOFNames();
  for(std::size_t i=0;i<dofs.size();i++) {
    this->addProvidedDOF(dofs[i].first,dofs[i].second);
    this->addDOFGrad(dofs[i].first,dofs[i].first);
  }

  // ********************
  // Build Basis Functions and Integration Rules
  // ********************
  
  this->setupDOFs(cell_data.baseCellDimension());
}

}

#endif
