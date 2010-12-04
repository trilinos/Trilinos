#ifndef MUELU_NEEDS_HPP
#define MUELU_NEEDS_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

namespace MueLu {
/*!
  @class NeedsObject
  @brief Base class for factories.
  Maintains just 2 things:
   - Ouput level status
   - A list of 'Needs' for the factory. For example, a restriction factory that transposes the tentative
     prolongator 'Needs' the prolongator factory to save this.

  Derived from Teuchos::VerboseObject.
*/
class Needs : public Teuchos::VerboseObject<Needs> {

/*
  //! Friendly print.  FIXME pretty print doesn't work ;(
  friend std::ostream& operator<<(std::ostream &os, Needs const &foo);
*/

  private:

    // TODO JG to JHU: why we cannot use direclty Teuchos::VerboseObject::getVerbLevel() and setVerbLevel() ?
    // FIXME JHU why doesn't  needs->setVerbLevel(Teuchos::VERB_NONE) change output?!
    //! Current output level
    Teuchos::EVerbosityLevel outputLevel_; //FIXME this state should already be stored in VerboseObject
    //! Prior output level
    Teuchos::EVerbosityLevel priorOutputLevel_;
    //! Stores number of references to a Need.
    Teuchos::ParameterList countTable_;
    //! Stores values associated with a Need.
    Teuchos::ParameterList dataTable_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:
    //@{ Constructors/Destructors.
    Needs() : out_(this->getOStream()) {
      countTable_.setName("countTable");
      dataTable_.setName("dataTable");
    }

    virtual ~Needs() {}
    //@}

    //@{ Access methods.
    //! Store need and its value. (This does not indicate that it's needed.)
    template <typename T>
    void Save(const std::string ename, const T &entry) {
      if ( !countTable_.isParameter(ename) ) {
        countTable_.set(ename,0);
      }
      dataTable_.set(ename,entry);
    } //Save

    //! Indicate that an object is needed.
    void Want(const std::string ename) {
      int currentCount = countTable_.get(ename,0);
      countTable_.set(ename,++currentCount);
    } //Want


    //! Get a need and decrement the storage counter.
    template <typename T>
    void Release(const std::string ename, T &value) {
      if (!countTable_.isParameter(ename)) {
        std::string msg =  "Release: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        value = dataTable_.get<T>(ename);
        int currentCount = countTable_.get(ename,0);
        if (currentCount == 1) {
          countTable_.remove(ename);
          dataTable_.remove(ename);
        } else {
          countTable_.set(ename,--currentCount);
        }
      }
    } //Release

    //! Get a need without decrementing the storage counter.
    //! This should be used only for local needs or for debugging purposes.
    template <typename T>
    void Get(const std::string ename, T &value) {
      if (!countTable_.isParameter(ename)) {
        Teuchos::OSTab tab(out_);
        std::string msg =  "Get: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        value = dataTable_.get<T>(ename);
      }
    } //Get
    //@}

    //@{ Utilities.
    //! NeedsTable raw print
    void RawPrint(std::ostream &os, Needs const &foo) {
      std::cout << "name | #requests" << std::endl;
      std::cout << "=================" << std::endl;
      std::cout << foo.countTable_ << std::endl << std::endl;
      std::cout << "name | value" << std::endl;
      std::cout << "============" << std::endl;
      std::cout << foo.dataTable_ << std::endl;
    } //RawPrint

    //! Return count for a need.
    int GetCount(const std::string ename) {
      if (!countTable_.isParameter(ename)) {
        std::string msg =  "GetCount: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        return countTable_.get(ename,0);
      }
    } //GetCount

    //! Test whether a need has been registered.
    //! Note: this tells nothing about whether the value exists.
    bool IsRegistered(const std::string ename) {
      if (countTable_.isParameter(ename)) return true;
      else                                return false;
    }
    //@}

    void SetOutputLevel(Teuchos::EVerbosityLevel outputLevel) {
      outputLevel_ = outputLevel;
    }

    int GetOutputLevel() {
      return outputLevel_;
    }

}; //class Needs

/*
//! NeedsTable pretty print
std::ostream& operator<<(std::ostream &os, Needs const &foo) {
  std::cout << "name  |  #requests  |  value" << std::endl;
  std::cout << "============================" << std::endl;
  for (Teuchos::ParameterList::ConstIterator param=foo.countTable_.begin(); param!=foo.countTable_.end() ; param++) {
    const std::string pname = foo.countTable_.name(param);
    const int& numRequests = foo.countTable_.get<int>(pname);
    const Teuchos::ParameterEntry &someEntry = foo.dataTable_.getEntry(pname);
    std::cout << pname << " | " << numRequests << " | " << someEntry << std::endl;
  }
  return os;
}
*/

} //namespace MueLu

#endif //ifndef MUELU_NEEDS_HPP
