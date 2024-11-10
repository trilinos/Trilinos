// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PRODUCT_IMPL_HPP
#define PANZER_PRODUCT_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Product<EvalT, Traits>::
Product(
  const Teuchos::ParameterList& p)
  : scaling(1.0)
{
  std::string product_name = p.get<std::string>("Product Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  if(p.isType<double>("Scaling"))
    scaling =  p.get<double>("Scaling");
  
  product = PHX::MDField<ScalarT>(product_name, data_layout);
  
  this->addEvaluatedField(product);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<const ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "Product Evaluator";
  this->setName(n);
}


template<int RANK, typename Scalar>
struct V_MultiplyFunctor {

  const PHX::MDField<Scalar> base_;
  const PHX::MDField<const Scalar> source_;

  V_MultiplyFunctor(PHX::MDField<Scalar>& base, const PHX::MDField<const Scalar>& source)
    : base_(base),source_(source) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const PHX::index_t & ind1) const
  {
    using idx_t = PHX::index_t;

    if (RANK == 1){
      base_(ind1) = base_(ind1)*source_(ind1);
    }
    else if (RANK == 2){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        base_(ind1,ind2) = base_(ind1,ind2)*source_(ind1,ind2);
    }
    else if (RANK == 3){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        for (idx_t ind3=0; ind3 < static_cast<idx_t>(base_.extent(2)); ind3++)
          base_(ind1,ind2,ind3) = base_(ind1,ind2,ind3)*source_(ind1,ind2,ind3);
    }
    else if (RANK == 4){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        for (idx_t ind3=0; ind3 < static_cast<idx_t>(base_.extent(2)); ind3++)
          for (idx_t ind4=0; ind4 < static_cast<idx_t>(base_.extent(3)); ind4++)
            base_(ind1,ind2,ind3,ind4) = base_(ind1,ind2,ind3,ind4)*source_(ind1,ind2,ind3,ind4);
    }
    else if (RANK == 5){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        for (idx_t ind3=0; ind3 < static_cast<idx_t>(base_.extent(2)); ind3++)
          for (idx_t ind4=0; ind4 < static_cast<idx_t>(base_.extent(3)); ind4++)
            for (idx_t ind5=0; ind5 < static_cast<idx_t>(base_.extent(4)); ind5++)
              base_(ind1,ind2,ind3,ind4,ind5) = base_(ind1,ind2,ind3,ind4,ind5)*source_(ind1,ind2,ind3,ind4,ind5);
    }
    else if (RANK == 6){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        for (idx_t ind3=0; ind3 < static_cast<idx_t>(base_.extent(2)); ind3++)
          for (idx_t ind4=0; ind4 < static_cast<idx_t>(base_.extent(3)); ind4++)
            for (idx_t ind5=0; ind5 < static_cast<idx_t>(base_.extent(4)); ind5++)
              for (idx_t ind6=0; ind6 < static_cast<idx_t>(base_.extent(5)); ind6++)
                base_(ind1,ind2,ind3,ind4,ind5,ind6) = base_(ind1,ind2,ind3,ind4,ind5,ind6)*source_(ind1,ind2,ind3,ind4,ind5,ind6);
    }
    else if (RANK == 7){
      for (idx_t ind2=0; ind2 < static_cast<idx_t>(base_.extent(1)); ind2++)
        for (idx_t ind3=0; ind3 < static_cast<idx_t>(base_.extent(2)); ind3++)
          for (idx_t ind4=0; ind4 < static_cast<idx_t>(base_.extent(3)); ind4++)
            for (idx_t ind5=0; ind5 < static_cast<idx_t>(base_.extent(4)); ind5++)
              for (idx_t ind6=0; ind6 < static_cast<idx_t>(base_.extent(5)); ind6++)
                for (idx_t ind7=0; ind7 < static_cast<idx_t>(base_.extent(6)); ind7++)
                  base_(ind1,ind2,ind3,ind4,ind5,ind6,ind7) = base_(ind1,ind2,ind3,ind4,ind5,ind6,ind7)*source_(ind1,ind2,ind3,ind4,ind5,ind6,ind7);
    }
  }
};

//**********************************************************************
template<typename EvalT, typename Traits>
void
Product<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{ 
  product.deep_copy(ScalarT(scaling));

  for (std::size_t i = 0; i < values.size(); ++i) {
    const auto length = product.extent(0);
    if (product.rank() == 1){
      Kokkos::parallel_for( length, V_MultiplyFunctor<1,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 2){
      Kokkos::parallel_for( length, V_MultiplyFunctor<2,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 3){
      Kokkos::parallel_for( length, V_MultiplyFunctor<3,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 4){
      Kokkos::parallel_for( length, V_MultiplyFunctor<4,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 5){
      Kokkos::parallel_for( length, V_MultiplyFunctor<5,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 6){
      Kokkos::parallel_for( length, V_MultiplyFunctor<6,ScalarT>(product,values[i]) );
    }
    else if (product.rank() == 7){
      Kokkos::parallel_for( length, V_MultiplyFunctor<7,ScalarT>(product,values[i]) );
    }
  }
}

//**********************************************************************

}

#endif
