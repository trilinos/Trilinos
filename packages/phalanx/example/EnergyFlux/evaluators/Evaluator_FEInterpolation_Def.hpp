#include "Teuchos_TestForException.hpp"

//**********************************************************************
template< typename ScalarT, typename Traits>
FEInterpolation<ScalarT, Traits>::
FEInterpolation(const Teuchos::ParameterList& p)
{ 
  std::vector<PHX::FieldTag>& scalar_node_ = 
    *(p.get< Teuchos::RCP<std::vector<PHX::FieldTag> > >("Scalar Node"));
  std::vector<PHX::FieldTag>& scalar_qp_ = 
    *(p.get< Teuchos::RCP< std::vector<PHX::FieldTag> > >("Scalar QP"));
  std::vector<PHX::FieldTag>& grad_scalar_node_ = 
    *(p.get< Teuchos::RCP< std::vector<PHX::FieldTag> > >("Grad Scalar Node"));
  std::vector<PHX::FieldTag>& grad_vector_qp_ = 
    *(p.get< Teuchos::RCP< std::vector<PHX::FieldTag> > >("Grad Vector QP"));

  static const std::string msg = "Error, the number of fields to be interpolated is not equal to the number of target fields.";

  TEST_FOR_EXCEPTION(scalar_node_.size() != scalar_qp_.size(), 
		     std::logic_error, msg);

  TEST_FOR_EXCEPTION(grad_scalar_node_.size() != grad_vector_qp_.size(),
		     std::logic_error, msg);

  s_n.resize(scalar_node_.size());
  s_qp.resize(scalar_qp_.size());
  g_n.resize(grad_scalar_node_.size());
  g_qp.resize(grad_vector_qp_.size());

  for (std::size_t i =0; i < s_qp.size(); ++i) {
    s_qp[i].setFieldTag(scalar_qp_[i]);
    this->addEvaluatedField(s_qp[i]);
  }
  
  for (std::size_t i =0; i < g_qp.size(); ++i) {
    g_qp[i].setFieldTag(grad_vector_qp_[i]);
    this->addEvaluatedField(g_qp[i]);
  }
  
  for (std::size_t i =0; i < s_n.size(); ++i) {
    s_n[i].setFieldTag(scalar_node_[i]);
    this->addDependentField(s_n[i]);
  }

  for (std::size_t i =0; i < g_n.size(); ++i) {
    g_n[i].setFieldTag(grad_scalar_node_[i]);
    this->addDependentField(g_n[i]);
  }

  this->setName("FEInterpolation");
}

//**********************************************************************
template<typename ScalarT, typename Traits>
FEInterpolation<ScalarT, Traits>::~FEInterpolation()
{ }

//**********************************************************************
template<typename ScalarT, typename Traits> 
void FEInterpolation<ScalarT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{

  for (std::size_t i =0; i < s_qp.size(); ++i)
    vm.setFieldData(s_qp[i]);

  for (std::size_t i =0; i < s_n.size(); ++i)
    vm.setFieldData(s_n[i]);

  for (std::size_t i =0; i < g_n.size(); ++i)
    vm.setFieldData(g_n[i]);

  for (std::size_t i =0; i < g_qp.size(); ++i)
    vm.setFieldData(g_qp[i]);
  
}

//**********************************************************************
template<typename ScalarT, typename Traits>
void FEInterpolation<ScalarT, Traits>::
evaluateFields(typename Traits::EvalData cell_data)
{ 

  //typename Traits::EvalData& cell_data = *d;

  // **********
  // QP values
  // ********** 

  // Loop over fields
  for (std::size_t var = 0; var < s_qp.size(); ++var) {

    // Loop over number of cells
    for (std::size_t cell = 0; cell < cell_data.size(); ++cell) {
      
      std::vector<double>& phi = cell_data[cell].getBasisFunctions();

      std::size_t qp_offset = cell * s_qp[var].fieldTag().dataLayout()->size();
      std::size_t node_offset = cell * s_n[var].fieldTag().dataLayout()->size();

      // Loop over quad points of cell
      for (std::size_t qp = 0; qp < s_qp[var].fieldTag().dataLayout()->size(); ++qp) {
	s_qp[var][qp_offset + qp] = 0.0;
	
	// Sum nodal contributions to qp
	for (std::size_t node = 0; node < s_n[var].fieldTag().dataLayout()->size(); ++node)
	  s_qp[var][qp_offset + qp] += 
	    phi[node] * s_n[var][node_offset + node];
	
      }

    }
    
  }
  
  // **********
  // QP Gradients
  // ********** 

  // Loop over fields
  for (std::size_t var = 0; var < g_qp.size(); ++var) {

    // Loop over number of cells
    for (std::size_t cell = 0; cell < cell_data.size(); ++cell) {
      
      std::vector< MyVector<double> >& grad_phi = 
	cell_data[cell].getBasisFunctionGradients();

      std::size_t qp_offset = 
	cell * g_qp[var].fieldTag().dataLayout()->size();
      std::size_t node_offset = 
	cell * g_n[var].fieldTag().dataLayout()->size();

      // Loop over quad points of cell
      for (std::size_t qp = 0; qp < g_qp[var].fieldTag().dataLayout()->size(); ++qp) {
	g_qp[var][qp_offset + qp] = MyVector<ScalarT>(0.0, 0.0, 0.0);
	
	// Sum nodal contributions to qp
	for (std::size_t node = 0; node < g_n[var].fieldTag().dataLayout()->size(); ++node)
	  g_qp[var][qp_offset + qp] += 
	    grad_phi[node] * g_n[var][node_offset + node];
	
      }

    }
    
  }
  
}

//**********************************************************************
