#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>
#include <sstream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_DefaultFactoryHandlerBase.hpp"

namespace MueLu {

  /*!
    @class Level
    @brief Class that holds all level-specific information.

    This class stores <tt>A</tt>, <tt>R</tt>, <tt>P</tt>, the presmother and the postsmoother
    explicitly.  All other data is stored in an associative list.
    See the Needs class for more information.
  */
  class Level : public Needs {

  private: 
    mutable int levelID_; // id number associated with level
    RCP<DefaultFactoryHandlerBase> defaultFactoryHandler_;

  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:

    //@{
    //! @name Constructors / Destructors
    Level() : levelID_(-1), out_(this->getOStream()) {
      //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating new uninitialized Level" << std::endl;
    }

    //! Constructor
    Level(RCP<DefaultFactoryHandlerBase>& defaultFactoryHandler) : levelID_(-1), defaultFactoryHandler_(defaultFactoryHandler), out_(this->getOStream()) { }

    //! Copy constructor.
    explicit Level(const Level& source) {
      levelID_ = source.levelID_;
      defaultFactoryHandler_ = source.defaultFactoryHandler_;
    }

    //@}

    //@{
    //! @name Build methods //TODO: merge with copy constructor?
    //! Builds a new Level object.
    RCP<Level> Build(std::ostream &os) { //TODO: why ostream in argument?
      return rcp( new Level(defaultFactoryHandler_) );
    }
    //@}

    virtual ~Level() {}

    void Print(std::ostream &os) {
      os << this << std::endl;
    }

    //@{
    //! @name Set methods.

    //! @brief Set level number.
    void SetLevelID(int i) const {
      levelID_ = i;
    }

    //! Set default factories (used internally by Hierarchy::SetLevel()).
    // Users should not use this method.
    void SetDefaultFactoryHandler(RCP<DefaultFactoryHandlerBase>& defaultFactoryHandler) {
      defaultFactoryHandler_ = defaultFactoryHandler;
    }

    //@}

    //@{
    //! @name Get methods.

    //! @brief Return level number.
    int GetLevelID() const { return levelID_; }

    //! Get default factory.
    const RCP<FactoryBase> & GetDefaultFactory(const std::string& varname) {
      TEST_FOR_EXCEPTION(defaultFactoryHandler_ == null, Exceptions::RuntimeError, "MueLu::Level::GetDefaultFactory(): no DefaultFactoryHandler.");
      return defaultFactoryHandler_->GetDefaultFactory(varname);
    }

    //@}

  }; //class Level

  std::ostream& operator<<(std::ostream& os, Level const &level);

} //namespace MueLu

#define MUELU_LEVEL_SHORT
#endif //ifndef MUELU_LEVEL_HPP
