#ifndef MUELU_FACTORYBASE_HPP
#define MUELU_FACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {
  class Level;

  static int generateUniqueFactoryId()
  {
    static int i = 0;
    ++i;
    return i;
  }

  //! Base class for factories (e.g., R, P, and A_coarse).
  class FactoryBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    FactoryBase()
    : id_(MueLu::generateUniqueFactoryId())
    { /*std::cout << "Factory id: " << id_ << std::endl;*/ }

    //! Destructor.
    virtual ~FactoryBase() { /*std::cout << "Deconstructor factory id: " << id_ << std::endl;*/}
    //@}

    //@{
    //! @name Build methods.

    virtual bool NewBuild(Level & requestedLevel) const = 0;

    virtual void callDeclareInput(Level & requestedLevel) const = 0;
    //@}

    //@{
    //! @name Access factory properties

    /// return unique factory id
    const int getID() const { return id_; };

    //@}

  private:
    const int id_;

  }; //class FactoryBase

} //namespace MueLu

#define MUELU_FACTORYBASE_SHORT
#endif //ifndef MUELU_FACTORYBASE_HPP
