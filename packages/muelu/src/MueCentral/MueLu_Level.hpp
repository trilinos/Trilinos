#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>

#include "Teuchos_RefCountPtr.hpp"
#include "Tpetra_CrsMatrix.hpp"    //FIXME replace with Cthulhu Operator
#include "Tpetra_Vector.hpp"       //FIXME replace with Cthulhu Vector
#include "MueLu_Smoother.hpp"

/*!
  @class Level
  @brief Multigrid level object.
*/

namespace MueLu {

template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
class Level {

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Operator;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Vector;

  //friend inline std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> &level);

  private: 

    Teuchos::RCP<Operator> A_;             // discretization operator
    Teuchos::RCP<Operator> R_;             // restriction operator
    Teuchos::RCP<Operator> P_;             // prolongator operator
    //TODO enable this data
    //Smoother PreSmoother_;   // smoother operator
    //Smoother PostSmoother_;  // smoother operator
    mutable int levelID_;                  // id number associated with level
    // auxilliary (optional) data
    //TODO enable this data
    Teuchos::RCP<Operator> AuxMatrix_;     // can be used instead of A_ to steer coarsening/sparsity pattern algorithms
    //TODO enable this data
    //AuxMatrixFunc_; // function used to regenerate AuxMatrix from A & coords instead of via RAP
    Teuchos::RCP<Operator> AuxMatP_;       // an alternative prolongator for projecting AuxMatrix
    Teuchos::RCP<Operator> AuxMatR_;       // an alternative restrictor for projecting AuxMatrix
    Teuchos::RCP<Vector> xCoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    Teuchos::RCP<Vector> yCoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    Teuchos::RCP<Vector> zCoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    Teuchos::RCP<Vector> coordP_;        // an alternative prolongator for projecting coordinates
    Teuchos::RCP<Vector> coordR_;        // an alternative restrictor for projecting coordinates
    //TODO enable this data
    //Graph Graph_;         // saving a graph for SA

  public:

    //@{
    //! @name Constructors / Destructors
    Level() : A_(0), R_(0), P_(0), levelID_(-1) {
      std::cout << "Constructing new unitialized Level" << std::endl;
    }

    //! Copy constructor.
    Level(Level const &Source) {
      std::cout << "Copy constructing existing Level" << std::endl;
      A_ = Source.A_;
      R_ = Source.R_;
      P_ = Source.P_;
      levelID_ = Source.levelID_;
    }
    //@}

    //@{
    //! @name Build methods
    static Teuchos::RCP< Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Build() {
      std::cout << "Building a Level" << std::endl;
      return Teuchos::rcp( new Level<Scalar,LocalOrdinal,GlobalOrdinal,Node>() );
    }
    //@}

    virtual ~Level() {}

    void Print(std::ostream &os) {
      os << this << std::endl;
    }

    //@{
    //! @name Set methods.

    void SetLevelID(int i) const {
      levelID_ = i;
    }

    void SetA(Teuchos::RCP<Operator> const A) {
      A_ = A;
    }

    void SetR(Teuchos::RCP<Operator> const R) {
      R_ = R;
    }

    void SetP(Teuchos::RCP<Operator> const P) {
      P_ = P;
    }

    void SetAuxMatrix(Teuchos::RCP<Operator> const AuxMatrix) {
      return AuxMatrix_;
    }

    void SetXCoords(Teuchos::RCP<Vector> x) {
      xCoords_ = x;
    }

    void SetYCoords(Teuchos::RCP<Vector> y) {
      yCoords_ = y;
    }

    void SetZCoords(Teuchos::RCP<Vector> z) {
      zCoords_ = z;
    }

    void SetAuxMatTransfers(Teuchos::RCP<Operator> &altP, Teuchos::RCP<Operator> &altR) {
      AuxMatP_ = altP;
      AuxMatR_ = altR;
    }

    void SetCoordTransfers(Teuchos::RCP<Operator> &altP, Teuchos::RCP<Operator> &altR) {
      coordP_ = altP;
      coordR_ = altR;
    }

    //@}

    //@{
    //! @name Get methods.

    int GetLevelID() const {
      return levelID_;
    }

    Teuchos::RCP<Operator> GetA() const {
      return A_;
    }

    Teuchos::RCP<Operator> GetR() const {
      return R_;
    }

    Teuchos::RCP<Operator> GetP() const {
      return P_;
    }

    Teuchos::RCP<Operator> GetAuxMatrix() {
      return AuxMatrix_;
    }

    Teuchos::RCP<Vector> GetXCoords() {
      return xCoords_;
    }

    Teuchos::RCP<Vector> GetYCoords() {
      return yCoords_;
    }

    Teuchos::RCP<Vector> GetZCoords() {
      return zCoords_;
    }

    void GetAuxMatTransfers(Teuchos::RCP<Operator> &altP, Teuchos::RCP<Operator> &altR) {
      altP = AuxMatP_;
      altR = AuxMatR_;
    }

    void GetCoordTransfers(Teuchos::RCP<Operator> &altP, Teuchos::RCP<Operator> &altR) {
      altP = coordP_;
      altR = coordR_;
    }
    void SetPreSmoother(Teuchos::RCP<Smoother> &preSmoo) {std::cout << "Need private data member for presmoother" << std::endl;}
    void SetPostSmoother(Teuchos::RCP<Smoother> &postSmoo) {std::cout << "Need private data member for postsmoother" << std::endl;}
/*
    //TODO ==================================================
    //TODO The following methods still need to be implemented.
    //TODO ==================================================
    SetPreSmoother
    SetPostSmoother
    GetPreSmoother
    GetPostSmoother
    GetGraph
    GetAuxMatrixFunc
    SetGraph
    SetAuxMatrixFunc
    ProjectInterface
*/
    //@}



}; //class Level

//Print function.  Not a friend b/c it uses only public interfaces for data access.
template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> const &level) {
  os << "Printing Level object " << level.GetLevelID() << std::endl;
  if (level.GetA() != Teuchos::null) os << *level.GetA() << std::endl;
  if (level.GetR() != Teuchos::null) os << *level.GetR() << std::endl;
  if (level.GetP() != Teuchos::null) os << *level.GetP() << std::endl;
  return os;
}

} //namespace MueLu
#endif //ifndef MUELU_LEVEL_HPP
