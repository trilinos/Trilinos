// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
