#ifndef MUELU_PFACTORY_HPP
#define MUELU_PFACTORY_HPP

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
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "PFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~PFactory() {}
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
    void setRestrictionMode(bool bRestrictionMode = false)    {
        restrictionMode_ = bRestrictionMode;
    }

    /// returns restrictionMode flag
    bool isRestrictionModeSet() { return restrictionMode_; }

    //@}
}; //class PFactory

} //namespace MueLu

#define MUELU_PFACTORY_SHORT

#endif //ifndef MUELU_PFACTORY_HPP
