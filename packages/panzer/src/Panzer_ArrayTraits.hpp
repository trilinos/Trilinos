
#ifndef PANZER_ARRAY_TRAITS_HPP
#define PANZER_ARRAY_TRAITS_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

  template<typename Scalar, typename Array> struct ArrayTraits;

  // Specialization for Intrepid::FieldContainer
  template<typename Scalar>
  struct ArrayTraits<Scalar,Intrepid::FieldContainer<Scalar> >
  {
    typedef int size_type;
  };

  // Specialization for Intrepid::FieldContainer
  template<typename Scalar>
  struct ArrayTraits<Scalar,PHX::MDField<Scalar> >
  {
    typedef typename PHX::MDField<Scalar>::size_type size_type;
  };

}

#endif
