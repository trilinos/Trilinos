#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"

//===========================================================================
//===========================================================================
template <typename EvalT,typename LO,typename GO>
panzer::AssemblyEngine<EvalT,LO,GO>::
AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> >& fmb) :
  m_field_manager_builder(fmb)
{ 

}

//===========================================================================
//===========================================================================
template <typename EvalT,typename LO,typename GO>
void panzer::AssemblyEngine<EvalT,LO,GO>::
evaluate(const panzer::AssemblyEngineInArgs& in)
{
  // *********************
  // Volumetric fill
  // *********************
  {
    const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& 
      worksets = m_field_manager_builder->getWorksets();
  
    // Loop over element blocks
    for (std::size_t block = 0; block < worksets.size(); ++block) {

      std::vector<panzer::Workset>& w = *worksets[block]; 

      Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = 
	m_field_manager_builder->getVolumeFieldManagers()[block];

      // Loop over worksets in this element block
      for (std::size_t i = 0; i < w.size(); ++i) {
	panzer::Workset& workset = w[i];

	workset.solution_vector = in.x;
	workset.solution_deriv_vector = in.dxdt;
	workset.residual_vector = in.f;
	workset.jacobian_matrix = in.j;

	fm->template evaluateFields<EvalT>(workset);
      }
    }

  }

  // *********************
  // BC fill
  // *********************
  // NOTE: We have to split neumann and dirichlet bcs since dirichlet
  // bcs overwrite equations where neumann sum into equations.  Make
  // sure all neumann are done before dirichlet.

  this->evaluateBCs(panzer::BCT_Neumann, in);
  this->evaluateBCs(panzer::BCT_Dirichlet, in);

  return;
}

//===========================================================================
//===========================================================================
template <typename EvalT,typename LO,typename GO>
void panzer::AssemblyEngine<EvalT,LO,GO>::
evaluateBCs(const panzer::BCType bc_type,
	    const panzer::AssemblyEngineInArgs& in)
{

  {
    const std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC>& bc_field_managers = 
      m_field_manager_builder->getBCFieldManagers();
  
    const std::map<panzer::BC,
      Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
      panzer::LessBC>& bc_worksets = 
      m_field_manager_builder->getBCWorksets();

    // Must do all neumann before all dirichlet so we need a double loop
    // here over all bcs
    typedef typename std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC>::const_iterator bcfm_it_type;

    typedef typename std::map<panzer::BC,
      Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
      panzer::LessBC>::const_iterator bcwkst_it_type;

    // loop over bcs
    for (bcfm_it_type bcfm_it = bc_field_managers.begin(); 
	 bcfm_it != bc_field_managers.end(); ++bcfm_it) {
      
      const panzer::BC& bc = bcfm_it->first;
      const std::map<unsigned,PHX::FieldManager<panzer::Traits> > bc_fm = 
	bcfm_it->second;

      bcwkst_it_type bc_wkst_it = bc_worksets.find(bc);

      TEST_FOR_EXCEPTION(bc_wkst_it == bc_worksets.end(), std::logic_error,
			 "Failed to find corresponding bc workset!");

      const std::map<unsigned,panzer::Workset>& bc_wkst = 
	*(bc_wkst_it->second);  

      // Only process bcs of the appropriate type (neumann or dirichlet)
      if (bc.bcType() == bc_type) {

	// Loop over local faces
	for (std::map<unsigned,PHX::FieldManager<panzer::Traits> >::const_iterator side = bc_fm.begin(); side != bc_fm.end(); ++side) {
	  
	  unsigned local_side_index = side->first;
	  PHX::FieldManager<panzer::Traits>& local_side_fm = 
	    const_cast<PHX::FieldManager<panzer::Traits>& >(side->second);
	  
	  std::map<unsigned,panzer::Workset>::const_iterator wkst_it = 
	    bc_wkst.find(local_side_index);
	  
	  TEST_FOR_EXCEPTION(wkst_it == bc_wkst.end(), std::logic_error,
			     "Failed to find corresponding bc workset side!");
	  
	  panzer::Workset& workset = 
	    const_cast<panzer::Workset&>(wkst_it->second); 
	  
	  
	  
	  // We have one workset per face
	  workset.solution_vector = in.x;
	  workset.solution_deriv_vector = in.dxdt;
	  workset.residual_vector = in.f;
	  workset.jacobian_matrix = in.j;
	  
	  local_side_fm.template evaluateFields<EvalT>(workset);
	  
	}
      }
    } 
  }

}
