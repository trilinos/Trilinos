#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "Teuchos_RefCountPtr.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

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
  typedef MueLu::OperatorFactory<Scalar,LO,GO,Node> OperatorFactory;
  typedef MueLu::SmootherFactory<Scalar,LO,GO,Node> SmootherFactory;

  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, Hierarchy<AA,BB,CC,DD> &hierarchy);
  //friend std::ostream& operator<< <>(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy);
  //friend std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy);

  private:

    std::vector<Teuchos::RCP<Level> > Levels_;

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
     void SetLevel(Teuchos::RCP<Level> const& level) {
       Levels_.push_back(level);
       level->SetLevelID(Levels_.size());
     }

     //! Retrieve a certain level from hierarchy.
     Teuchos::RCP<Level>& GetLevel(int levelID) {
       return Levels_[levelID];
     }

     //FIXME should return status
     //FIXME also calculate complexity here
     void FullPopulate(Teuchos::RCP<OperatorFactory> PFact,
                       Teuchos::RCP<OperatorFactory> RFact=Teuchos::null,
                       Teuchos::RCP<OperatorFactory> AcFact=Teuchos::null,
                       Teuchos::RCP<SmootherFactory> SmooFact=Teuchos::null,
                       int startLevel=0, int numDesiredLevels=10 /*,Needs*/)
     {
       std::cout << "Hierarchy::FullPopulate()" << std::endl;
       bool goodBuild=true;
       int i = startLevel;
       while (i < startLevel + numDesiredLevels - 1)
       {
         //TODO MueMat has a check for an empty level at position i+1
         //TODO I'm not sure how this can happen, but if it does, we must insert code to handle
         //TODO it similarly to MueMat.
         if (PFact != Teuchos::null /*|| isempty(Levels_[i+1] */) {
           if ((i+1) >= (int) Levels_.size())
             Levels_.push_back( Levels_[i]->Build() );
           Levels_[i+1]->SetLevelID(i+1);
           goodBuild = PFact->Build(*(Levels_[i]),*(Levels_[i+1]) /*,MySpecs*/);
         }
         if ((int)Levels_.size() <= i) goodBuild=false; //TODO is this the right wasy to cast?
         if (!goodBuild) /*TODO make Levels_ be length i*/;
         if (RFact != Teuchos::null)
           if ( !RFact->Build(*(Levels_[i]),*(Levels_[i+1]) /*,MySpecs*/) ) {
             //make Levels_ be length i
             break;
           }
         if (AcFact != Teuchos::null)
           if ( !AcFact->Build(*(Levels_[i]),*(Levels_[i+1]) /*,MySpecs*/) ) {
             //make Levels_ be length i
           break;
           }
         if (SmooFact != Teuchos::null) {
           Teuchos::RCP<Smoother> preSm, postSm;
           SmooFact->Build(preSm,postSm,*(Levels_[i]) /*,MySpecs*/);
           if (preSm != Teuchos::null) Levels_[i]->SetPreSmoother(preSm);
           if (postSm != Teuchos::null) Levels_[i]->SetPostSmoother(postSm);
         }
         ++i;
       } //while
     } //FullPopulate()

     //FIXME should return status
     void SetSmoothers() { std::cout << "Hierarchy::SetSmoothers()" << std::endl; }

     //FIXME should return status
     void FillHierarchy(Teuchos::RCP<OperatorFactory> PFact,
                        Teuchos::RCP<OperatorFactory> RFact=Teuchos::null,
                        Teuchos::RCP<OperatorFactory> AcFact=Teuchos::null,
                        int startLevel=0, int numDesiredLevels=10 /*,Needs*/)
     {
       if (PFact == Teuchos::null) {throw("FillHierarchy: must supply at least a Prolongator factory");}
       if (RFact == Teuchos::null) RFact = Teuchos::rcp(new TransPFactory<Scalar,LO,GO,Node>());
       if (AcFact == Teuchos::null) AcFact = Teuchos::rcp(new RAPFactory<Scalar,LO,GO,Node>());

       Teuchos::RCP<SmootherFactory> SmooFact=Teuchos::null;
       FullPopulate(PFact,RFact,AcFact,SmooFact,startLevel,numDesiredLevels);
     }

     //FIXME should return solution vector
     void Iterate() { std::cout << "Hierarchy::Iterate()" << std::endl; }

   //@}
    
}; //class Hierarchy

template<class Scalar,class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LO,GO,Node> &hierarchy) {
  os << "Printing Hierarchy object" << std::endl;
  typename std::vector< Teuchos::RCP<Level<Scalar,LO,GO,Node> > >::const_iterator i;
  for (i=hierarchy.Levels_.begin(); i != hierarchy.Levels_.end(); ++i)
    os << *(*i) << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_HIERARCHY_HPP
