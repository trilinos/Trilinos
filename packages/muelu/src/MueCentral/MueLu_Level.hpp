#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>

#include "Teuchos_RefCountPtr.hpp"
#include "MueLu_Needs.hpp"
#include "Cthulhu_Operator.hpp"
#include "Cthulhu_Vector.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include <sstream>

// JG TODO: remove template parameters

namespace MueLu {

  template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
  class SmootherPrototype;

  /*!
    @class Level
    @brief Class that holds all level-specific information.

    This class stores <tt>A</tt>, <tt>R</tt>, <tt>P</tt>, the presmother and the postsmoother
    explicitly.  All other data is stored in an associative list.
    See the Needs class for more information.
  */
  //template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  //class Level : public Teuchos::VerboseObject<Level<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> > {
  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Level : public Needs {

#include "MueLu_UseShortNames.hpp"

    //friend inline std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &level);

  private: 

    Teuchos::RCP<Operator> A_;             // discretization operator
    Teuchos::RCP<Operator> R_;             // restriction operator
    Teuchos::RCP<Operator> P_;             // prolongator operator
    Teuchos::RCP<SmootherPrototype> PreSmoother_;   // smoother operator
    Teuchos::RCP<SmootherPrototype> PostSmoother_;  // smoother operator
    mutable int levelID_;                  // id number associated with level

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //@{
    //! @name Constructors / Destructors
    Level() : A_(0), R_(0), P_(0), levelID_(-1), out_(this->getOStream())
    {
      //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating new uninitialized Level" << std::endl;
    }

    //! Copy constructor.
    Level(Level const &Source) : out_(this->getOStream()) {
      //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Copy constructing existing Level" << std::endl;
      A_ = Source.A_;
      R_ = Source.R_;
      P_ = Source.P_;
      levelID_ = Source.levelID_;
    }
    //@}

    //@{
    //! @name Build methods
    //! Builds a new Level object.
    static Teuchos::RCP<Level> Build(std::ostream &os) {
      return Teuchos::rcp( new Level() );
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

    //! @brief Set matrix.
    void SetA(Teuchos::RCP<Operator> const A) {
      A_ = A;
    }

    //! @brief Set restriction matrix.
    void SetR(Teuchos::RCP<Operator> const R) {
      R_ = R;
    }

    //! @brief Set prolongation matrix.
    void SetP(Teuchos::RCP<Operator> const P) {
      P_ = P;
    }

    //! @brief Set presmoother.
    void SetPreSmoother(Teuchos::RCP<SmootherPrototype> &preSmoo) {
      PreSmoother_ = preSmoo;
    }

    //! @brief Set postsmoother.
    void SetPostSmoother(Teuchos::RCP<SmootherPrototype> &postSmoo) {
      PostSmoother_ = postSmoo;
    }

    //@}

    //@{
    //! @name Get methods.

    //! @brief Return level number.
    int GetLevelID() const {
      return levelID_;
    }

    //! @brief Return matrix.
    Teuchos::RCP<Operator> GetA() const {
      if (A_ != Teuchos::null)
        return A_;
      else {
        std::ostringstream buf;
        buf << levelID_;
        std::string msg = "Level " + buf.str() + ": A is not set";
        throw(std::logic_error(msg));
      }
    }

    //! @brief Return restriction matrix.
    Teuchos::RCP<Operator> GetR() const {
      if (R_ != Teuchos::null)
        return R_;
      else {
        std::ostringstream buf;
        buf << levelID_;
        std::string msg = "Level " + buf.str() + ": R is not set";
        throw(std::logic_error(msg));
      }
    }

    //! @brief Return prolongation matrix.
    Teuchos::RCP<Operator> GetP() const {
      if (P_ != Teuchos::null)
        return P_;
      else {
        std::ostringstream buf;
        buf << levelID_;
        std::string msg = "Level " + buf.str() + ": P is not set";
        throw(std::logic_error(msg));
      }
    }

    /*! @brief Return presmoother.

        Does not throw exception if smoother is null.
    */
    Teuchos::RCP<SmootherPrototype> GetPreSmoother() const {
      return PreSmoother_;
    }

    /*! @brief Return postsmoother.

        Does not throw exception if smoother is null.
    */
    Teuchos::RCP<SmootherPrototype> GetPostSmoother() const {
      return PostSmoother_;
    }

    /*! @brief Return coarsest solver.

        Does not throw exception if solver is null.
    */

    //@}

  }; //class Level

  //Print function.  Not a friend b/c it uses only public interfaces for data access.
  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> const &level) {
    os << "Printing Level object " << level.GetLevelID() << std::endl;
    if (level.GetA() != Teuchos::null) os << *level.GetA() << std::endl;
    if (level.GetR() != Teuchos::null) os << *level.GetR() << std::endl;
    if (level.GetP() != Teuchos::null) os << *level.GetP() << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_LEVEL_SHORT

#endif //ifndef MUELU_LEVEL_HPP
