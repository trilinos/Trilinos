// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_OPERATOR_TRAITS_HPP
#define ANASAZI_OPERATOR_TRAITS_HPP

/*!     \file AnasaziOperatorTraits.hpp
        \brief Virtual base class which defines basic traits for the operator type
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"


namespace Anasazi {


  //! \brief Exceptions thrown to signal error in operator application.
  class OperatorError : public AnasaziError
  {public: OperatorError(const std::string& what_arg) : AnasaziError(what_arg) {}};


  /*! \brief This is the default struct used by OperatorTraits<ScalarType, MV, OP> class to produce a
      compile time error when the specialization does not exist for operator type <tt>OP</tt>.
  */
  template< class ScalarType, class MV, class OP >
  struct UndefinedOperatorTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    /*! \note Any attempt to compile this function results in a compile time error.  This means
      that the template specialization of Anasazi::OperatorTraits class does not exist for type
      <tt>OP</tt>, or is not complete.
    */
    static inline void notDefined() { return OP::this_type_is_missing_a_specialization(); };
  };


  /*!  \brief Virtual base class which defines basic traits for the operator type.

       An adapter for this traits class must exist for the <tt>MV</tt> and <tt>OP</tt> types.
       If not, this class will produce a compile-time error.

       \ingroup anasazi_opvec_interfaces
  */
  template <class ScalarType, class MV, class OP>
  class OperatorTraits 
  {
  public:
    
    //! @name Operator application method.
    //@{ 
    
    //! Application method which performs operation <b>y = Op*x</b>. An OperatorError exception is thrown if there is an error.
    static void Apply ( const OP& Op, 
                        const MV& x, 
                        MV& y )
    { UndefinedOperatorTraits<ScalarType, MV, OP>::notDefined(); };
    
    //@}
    
  };
  
} // end Anasazi namespace

#endif // ANASAZI_OPERATOR_TRAITS_HPP

// end of file AnasaziOperatorTraits.hpp
