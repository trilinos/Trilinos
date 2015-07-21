#ifndef MDARRAY_KOKKOS_FUNCTORS_HPP
#define MDARRAY_KOKKOS_FUNCTORS_HPP

#include "Kokkos_View.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Phalanx_MDField.hpp"
#include "Cell.hpp"




// ************* FEInterpolation *************

template < class MDFieldType2_1, class MDFieldType2_2, class MDFieldType3>
class Kokkos_FEInterpolation {

 //input
 MDFieldType2_1 val_qp_;
 MDFieldType2_2 val_node_;
 int num_qp_;
 int num_dim_;
 int num_nodes_;

 //output
 MDFieldType3 val_grad_qp_;

 //temporary
 Kokkos::View <double**, PHX::Device> phi_; 
 Kokkos::View <double***, PHX::Device> grad_phi_;
 
 public:

 Kokkos_FEInterpolation   ( MDFieldType2_1 &val_qp,
                            MDFieldType2_2 &val_node,
                            MDFieldType3 &val_grad_qp,
                            int num_qp,
                            int num_dim,
                            int num_nodes)
                        : val_qp_(val_qp)
                        , val_node_(val_node)
                        , val_grad_qp_(val_grad_qp)
                        , num_qp_ (num_qp)
 			, num_dim_(num_dim)
			, num_nodes_(num_nodes)
 {
   MyCell cell_data_;
   phi_=cell_data_.getBasisFunctions();
   grad_phi_=cell_data_.getBasisFunctionGradients();
 }

 KOKKOS_INLINE_FUNCTION
 void operator () (const int i) const
 {
    for (int qp = 0; qp < num_qp_; ++qp) {
      val_qp_(i,qp) = 0.0;
      for (int dim = 0; dim < num_dim_; ++dim)
        val_grad_qp_(i,qp,dim) = 0.0;
      // Sum nodal contributions to qp
      for (int node = 0; node < num_nodes_; ++node) {
         val_qp_(i,qp) += phi_(qp, node) * val_node_(i,node);
         for (int dim = 0; dim < num_dim_; ++dim)
           val_grad_qp_(i,qp,dim) += grad_phi_(qp, node, dim) * val_node_(i,node);
         }  
     }
  }
};


#endif
