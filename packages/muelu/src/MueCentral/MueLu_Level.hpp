#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Tpetra_CrsMatrix.hpp"    //FIXME replace with Cthulhu Operator
#include <iostream>

/*!
  @class Level
  @brief Multigrid level object.
*/

namespace MueLu {

template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
class Level {

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Operator;

  //friend inline std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> &level);

  private: 

    Teuchos::RCP<Operator> A_;             // discretization operator
    Teuchos::RCP<Operator> R_;             // restriction operator
    Teuchos::RCP<Operator> P_;             // prolongator operator
    //TODO enable this data
    //Smoother PreSmoother_;   // smoother operator
    //Smoother PostSmoother_;  // smoother operator
    int levelID_;                  // id number associated with level
    // auxilliary (optional) data
    //TODO enable this data
    //Operator AuxMatrix_;     // can be used instead of A_ to steer coarsening/sparsity pattern algorithms
    //AuxMatrixFunc_; // function used to regenerate AuxMatrix from A & coords instead of via RAP
    //Operator AuxMatP_;       // an alternative prolongator for projecting AuxMatrix
    //Operator AuxMatR_;       // an alternative restrictor for projecting AuxMatrix
    //Vector xcoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    //Vector ycoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    //Vector zcoords_;       // coordinates can be used to steer coarsening/sparsity or for visualization
    //Vector coordP_;        // an alternative prolongator for projecting coordinates
    //Vector coordR_;        // an alternative restrictor for projecting coordinates
    //Graph Graph_;         // saving a graph for SA

  public:

    //@{
    //! @name Constructors / Destructors
    //@}
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

    virtual ~Level() {}

    void Print(std::ostream &os) {
      os << this << std::endl;
    }

    //@{
    //! @name Get/Set methods.
    void SetLevelID(int i) {
      levelID_ = i;
    }

    int GetLevelID() {
      return levelID_;
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

    Teuchos::RCP<Operator> GetA() {
      return A_;
    }

    Teuchos::RCP<Operator> GetR() {
      return R_;
    }

    Teuchos::RCP<Operator> GetP() {
      return P_;
    }

/*
    //TODO ==================================================
    //TODO The following methods still need to be implemented.
    //TODO ==================================================
    SetPreSmoother
    SetPostSmoother
    GetPreSmoother
    GetPostSmoother
    GetGraph
    GetAuxMatrix
    GetAuxMatTransfers
    GetAuxMatrixFunc
    GetCoordTransfers
    GetXCoords
    GetYCoords
    GetZCoords
    SetGraph
    SetAuxMatrix
    SetAuxMatTransfers
    SetcoordTransfers
    SetAuxMatrixFunc
    SetXCoords
    SetYCoords
    SetZCoords
    ProjectInterface
*/


    //@}

}; //class Level

//Print function.  Not a friend b/c it uses only public interfaces for data access.
template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
std::ostream& operator<<(std::ostream& os, Level<Scalar,LocalOrdinal,GlobalOrdinal,Node> &level) {
  os << "printing a Level object" << std::endl;
  os << "levelID = " << level.GetLevelID() << std::endl;
  if (level.GetA() != Teuchos::null) os << *level.GetA() << std::endl;
  if (level.GetR() != Teuchos::null) os << *level.GetR() << std::endl;
  if (level.GetP() != Teuchos::null) os << *level.GetP() << std::endl;
  return os;
}

} //namespace MueLu
#endif //ifndef MUELU_LEVEL_HPP
