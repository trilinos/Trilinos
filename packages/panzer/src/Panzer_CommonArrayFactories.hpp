#ifndef PANZER_COMMON_ARRAY_FACTORIES_HPP
#define PANZER_COMMON_ARRAY_FACTORIES_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Phalanx_MDField.hpp"

#include <string>

//
// This file contains several common array factories
// useful for build arrays through the BasisValues and
// IntegrationValues classes. In particular these objects
// are used in the <code>setupArrays</code> functions.
// Because these class are used as a template argument the
// types and names used are very specific to the BasisValues
// interface.
//

namespace panzer {
  
  /** Implementation for intrepid field container factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects. Notice in this case the string
    * argument is not used.
    */
  class IntrepidFieldContainerFactory {
  public:
     template <typename Scalar,typename T0>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;
  };

  /** Implementation for MDField array factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects.
    */
  class MDFieldArrayFactory {
  public:
     /** Build fields with no prefix, will simply use the string
       * passed into <code>buildArray</code> to name the fields.
       */
     MDFieldArrayFactory() : prefix_("") {}

     /** Build fields with a prefix, will use the string
       * passed into <code>buildArray</code> prefixed with the
       * argument to this constructor to name the fields.
       */
     MDFieldArrayFactory(const std::string & prefix) : prefix_(prefix) {}

 
     template <typename Scalar,typename T0>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;

  private:
     std::string prefix_;     
  };

} // namespace panzer

#include "Panzer_CommonArrayFactories_impl.hpp"

#endif
