
#ifndef PANZER_WORKSET_HPP
#define PANZER_WORKSET_HPP

#include <cstddef>
#include <vector>
#include <map>
#include "Teuchos_ArrayRCP.hpp"
#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"
#include "Panzer_BasisValues.hpp"
#include "Panzer_IntegrationValues.hpp"
#include "Intrepid_FieldContainer.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;

namespace panzer {

  struct Workset {
    
    std::size_t num_cells;
    std::vector<std::size_t> cell_local_ids;
    Intrepid::FieldContainer<double> cell_vertex_coordinates;
    std::string block_id;

    //! Value correspondes to integration order.  Use the offest for indexing.
    Teuchos::RCP< std::vector<int> > ir_degrees;
    
    std::vector<Teuchos::RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > > int_rules;
    
    //! Value corresponds to basis type.  Use the offest for indexing.
    Teuchos::RCP< std::vector<std::string> > basis_names;

    //! Static basis function data, key is basis name, value is index in the static_bases vector
    std::vector<Teuchos::RCP< panzer::BasisValues<double,Intrepid::FieldContainer<double> > > > bases;
    
    Teuchos::RCP<Epetra_Vector> solution_vector;
    Teuchos::RCP<Epetra_Vector> solution_deriv_vector;
    Teuchos::RCP<Epetra_Vector> residual_vector;
    Teuchos::RCP<Epetra_CrsMatrix> jacobian_matrix;

  };

} // namespace panzer

#endif
