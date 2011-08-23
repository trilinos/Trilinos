#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>
#include <sstream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Needs.hpp"

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
    mutable int levelID_;                  // id number associated with level

  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:

    //@{
    //! @name Constructors / Destructors
    Level() : levelID_(-1), out_(this->getOStream()) {
      //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating new uninitialized Level" << std::endl;
    }

    //! Copy constructor.
    Level(Level const &Source) : out_(this->getOStream()) {
      //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Copy constructing existing Level" << std::endl;
      levelID_ = Source.levelID_;
    }
    //@}

    //@{
    //! @name Build methods
    //! Builds a new Level object.
    static RCP<Level> Build(std::ostream &os) {
      return rcp( new Level() );
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

    //@}

    //@{
    //! @name Get methods.

    //! @brief Return level number.
    int GetLevelID() const { return levelID_; }

    //@}

  }; //class Level

  std::ostream& operator<<(std::ostream& os, Level const &level);

} //namespace MueLu

#define MUELU_LEVEL_SHORT

#endif //ifndef MUELU_LEVEL_HPP
