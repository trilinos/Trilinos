// @HEADER
// @HEADER

#ifndef PANZER_INTREPID_BASIS_FACTORY_H
#define PANZER_INTREPID_BASIS_FACTORY_H

#include <sstream>
#include <string>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Intrepid_Basis.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"

namespace panzer {

  template <typename ScalarT, typename ArrayT>
    Teuchos::RCP<Intrepid::Basis<ScalarT,ArrayT> >
    createIntrepidBasis(const std::string type, int cell_dimension) {
    
    Teuchos::RCP<Intrepid::Basis<ScalarT,ArrayT> > basis;
    
    if (cell_dimension == 3) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C1_FEM<ScalarT,ArrayT> );
      else if (type == "Q2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C2_FEM<ScalarT,ArrayT> );
      else if (type == "T1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C1_FEM<ScalarT,ArrayT> );
      else if (type == "T2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C2_FEM<ScalarT,ArrayT> );
    }
    else if (cell_dimension == 2) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<ScalarT,ArrayT> );
      else if (type == "Q2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<ScalarT,ArrayT> );
      if (type == "T1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C1_FEM<ScalarT,ArrayT> );
      else if (type == "T2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C2_FEM<ScalarT,ArrayT> );
    }
    else if (cell_dimension == 1) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_LINE_C1_FEM<ScalarT,ArrayT> );
    }

    if (Teuchos::is_null(basis)) {
      std::ostringstream s;
      s << "Failed to create basis.  Either the basis type in the input file, \"" << type << "\" is unsupported or we are out of memory!" << std::endl;
      TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
			 s.str());
    }

    return basis;
  }

  inline 
  std::string getBasisName(const std::string type) {
    const std::string name = "Basis " + type;
    return  name;
  }
  
  inline
  std::string getD1BasisName(const std::string type) {
    const std::string name = "Grad_Basis " + type;
    return  name;
  }

  inline
  std::string getD2BasisName(const std::string type) {
    const std::string name = "D2_Basis " + type;
    return  name;
  }
  
}


#endif
