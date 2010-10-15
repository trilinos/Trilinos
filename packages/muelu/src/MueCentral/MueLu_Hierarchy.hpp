#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "MueLu_Level.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"

namespace MueLu {

/*!
  @class Hierarchy
  @brief Provides methods to build a multigrid hierarchy and apply multigrid cycles.

  Allows users to manually populate operators at different levels within 
  a multigrid method and push them into the hierarchy via SetBucket() 
  and/or to supply factories for automatically generating prolongators, 
  restrictors, and coarse level discretizations.  Additionally contains 
  a V-cycle apply method.
*/
template<class Scalar,class LO, class GO, class Node>
class Hierarchy {

  typedef MueLu::Level<Scalar,LO,GO,Node> Level;

  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, Hierarchy<AA,BB,CC,DD> &hierarchy);
  //friend std::ostream& operator<< <>(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy);
  //friend std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy);

  private:

    std::vector<Level> Levels_;

  public:

  //! @name Constructors/Destructors
  //@{

    //! Default constructor.
    Hierarchy() {}

    //! Copy constructor.
    Hierarchy(Hierarchy const &inHierarchy) {
      std::cerr << "Not implemented yet." << std::endl;
    }

    //! Destructor.
    virtual ~Hierarchy() {}

   //@}

   //! Set/Get Methods.
   //@{

     //! Assign a level to hierarchy.
     void SetLevel(Level const& level) {
       Levels_.push_back(level);
       level.SetLevelID(Levels_.size());
     }

     //! Retrieve a certain level from hierarchy.
     Level& GetLevel(int levelID) {
       return Levels_[levelID];
     }

     //FIXME should return status
     void FullPopulate() { std::cout << "Hierarchy::FullPopulate()" << std::endl; }

     //FIXME should return status
     void SetSmoothers() { std::cout << "Hierarchy::SetSmoothers()" << std::endl; }

     //FIXME should return status
     void FillHierarchy(Teuchos::RCP<BaseFactory> PFact,
                        Teuchos::RCP<BaseFactory> RFact=Teuchos::null,
                        Teuchos::RCP<BaseFactory> AcFact=Teuchos::null,
                        int startLevel=1, int numDesiredLevels=10 /*,Needs*/)
     {
       if (PFact == Teuchos::null) {throw("FillHierarchy: must supply at least a Prolongator factory");}
       if (RFact == Teuchos::null) RFact = TransPFactory();
       if (AcFact == Teuchos::null) AcFact = RAPFactory();

       FullPopulate();
     }

     //FIXME should return solution vector
     void Iterate() { std::cout << "Hierarchy::Iterate()" << std::endl; }

   //@}
    
}; //class Hierarchy

template<class Scalar,class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy) {
  os << "Printing Hierarchy object" << std::endl;
  typename std::vector< Level<Scalar,LO,GO,Node> >::const_iterator i;
  for (i=hierarchy.Levels_.begin(); i != hierarchy.Levels_.end(); ++i)
    os << *i << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_HIERARCHY_HPP
