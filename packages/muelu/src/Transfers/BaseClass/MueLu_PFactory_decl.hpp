#ifndef MUELU_PFACTORY_DECL_HPP
#define MUELU_PFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

  class Level;

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

class PFactory : public TwoLevelFactoryBase {


  protected:

     bool restrictionMode_;  //< true, if PFactory is used for generating the restriction operator

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PFactory() :  restrictionMode_(false)
    ;

    //! Destructor.
    virtual ~PFactory() ;
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    virtual void BuildP(Level &fineLevel, Level &coarseLevel) const = 0;
    //@}


    //! @name Restriction mode
    //@{

    /// switch prolongator factory to restriction mode
    /// if set to true, the prolongation factory generates a restriciton operator instead of a prolongation operator
    void setRestrictionMode(bool bRestrictionMode = false)    ;

    /// returns restrictionMode flag
    bool isRestrictionModeSet() ;

    //@}
}; //class PFactory

} //namespace MueLu

#define MUELU_PFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_PFACTORY_DECL_HPP
