#ifndef MUELU_HIERARCHY_HPP
#define MUELU_HIERARCHY_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {
/*!
  @class Hierarchy
  @brief Provides methods to build a multigrid hierarchy and apply multigrid cycles.

  Allows users to manually populate operators at different levels within 
  a multigrid method and push them into the hierarchy via SetLevel() 
  and/or to supply factories for automatically generating prolongators, 
  restrictors, and coarse level discretizations.  Additionally contains 
  a V-cycle apply method.
*/
template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class Hierarchy : public Teuchos::VerboseObject<Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > {

#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, Hierarchy<AA,BB,CC,DD,EE> &hierarchy);

  private:

    //! vector of Level objects
    std::vector<Teuchos::RCP<Level> > Levels_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

  //! @name Constructors/Destructors
  //@{

    //! Default constructor.
    Hierarchy() : out_(this->getOStream()) {}

    //! Copy constructor.
    Hierarchy(Hierarchy const &inHierarchy) {
      std::cerr << "Not implemented yet." << std::endl;
    }

    //! Destructor.
    virtual ~Hierarchy() {}

   //@}

   //! @name Set/Get Methods.
   //@{

     //! Assign a level to hierarchy.
     void SetLevel(Teuchos::RCP<Level> const& level) {
       Levels_.push_back(level);
       level->SetLevelID(Levels_.size());
     }

     //! Retrieve a certain level from hierarchy.
     Teuchos::RCP<Level>& GetLevel(int const levelID) {
       return Levels_[levelID];
     }
   //@}

   //! @name Populate Methods.
   //@{

     /*!
       @brief Constructs components of the hierarchy.
       FIXME should return status
       FIXME also calculate complexity here

       The only required factory is for the prolongator.  Factories to build the smoothers,
       restriction matrices, and coarse level matrices are optional.  Default behavior is to
       ignore any empty factories.  FillHierarchy invokes this method.
     */
     void FullPopulate(Teuchos::RCP<OperatorFactory> PFact,
                       Teuchos::RCP<OperatorFactory> RFact=Teuchos::null,
                       Teuchos::RCP<OperatorFactory> AcFact=Teuchos::null,
                       Teuchos::RCP<SmootherFactory> SmooFact=Teuchos::null,
                       int startLevel=0, int numDesiredLevels=10 )
     {
       Teuchos::OSTab tab(out_);
       MueLu_cout(Teuchos::VERB_HIGH) << "Hierarchy::FullPopulate()" << std::endl;
       bool goodBuild=true;
       int i = startLevel;
       while (i < startLevel + numDesiredLevels - 1)
       {
         //TODO MueMat has a check for an empty level at position i+1
         //TODO I'm not sure how this can happen, but if it does, we must insert code to handle
         //TODO it similarly to MueMat.
         if (PFact != Teuchos::null) {
           if ( (i+1) >= (int) Levels_.size() || Levels_[i+1] == Teuchos::null ) {
             Teuchos::OSTab tab(out_);
             Levels_.push_back( Levels_[i]->Build(*out_) );
           }
           Levels_[i+1]->SetLevelID(i+1);
           goodBuild = PFact->Build(*(Levels_[i]),*(Levels_[i+1]) /*,MySpecs*/);
         } //if (PFact != Teuchos::null)
         if ((int)Levels_.size() <= i) goodBuild=false; //TODO is this the right way to cast?
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
           Teuchos::RCP<SmootherPrototype> preSm, postSm;
           SmooFact->Build(Levels_[i], preSm, postSm);
           if (preSm != Teuchos::null) Levels_[i]->SetPreSmoother(preSm);
           if (postSm != Teuchos::null) Levels_[i]->SetPostSmoother(postSm);
         }
         ++i;
       } //while
     } //FullPopulate()

     /*! @brief Construct smoothers on all levels.
       TODO should return status
     */
     void SetSmoothers(RCP<SmootherFactory> smooFact=Teuchos::null, LO startLevel=0, LO numDesiredLevels=-1) {
       Teuchos::OSTab tab(out_);
       MueLu_cout(Teuchos::VERB_HIGH) << "Hierarchy::SetSmoothers()" << std::endl;
       if (smooFact == Teuchos::null) {
         throw(Exceptions::NotImplemented("No default smoother is defined"));
       }
       if (numDesiredLevels == -1)
         numDesiredLevels = Levels_.size()-startLevel+1;
       //TODO FullPopulate should return boolean status
       FullPopulate(Teuchos::null,Teuchos::null,Teuchos::null,smooFact,startLevel,numDesiredLevels);
     } //SetSmoothers()

     /*!
       @brief Constructs components of the hierarchy.
       TODO should return status

       The only required factory is for the prolongator.  Factories to build the restriction and
       coarse level matrices are optional.  Default behavior is to ignore any empty factories.
       Internally, this calls FullPopulate.

     */
     void FillHierarchy(Teuchos::RCP<OperatorFactory> PFact,
                        Teuchos::RCP<OperatorFactory> RFact=Teuchos::null,
                        Teuchos::RCP<OperatorFactory> AcFact=Teuchos::null,
                        int startLevel=0, int numDesiredLevels=10 /*,Needs*/)
     {
       if (PFact == Teuchos::null) {throw(std::logic_error("FillHierarchy: must supply at least a Prolongator factory"));}
       if (RFact == Teuchos::null) RFact = Teuchos::rcp(new TransPFactory());
       if (AcFact == Teuchos::null) AcFact = Teuchos::rcp(new RAPFactory());

       Teuchos::RCP<SmootherFactory> SmooFact=Teuchos::null;
       FullPopulate(PFact,RFact,AcFact,SmooFact,startLevel,numDesiredLevels);
     }

     /*!
       @brief Apply the multigrid preconditioner.
     */

     void Iterate(MultiVector const &rhs, LO nIts, MultiVector const &solution,
                  bool InitialGuessIsZero=false, CycleType const Cycle=VCYCLE, LO const startLevel=1)
     {
       Teuchos::OSTab tab(out_);
       MueLu_cout(Teuchos::VERB_HIGH) << "Hierarchy::Iterate()" << std::endl;

     /*
       for (LO i=0; i<nIts; i++) {
         RCP<Level> Fine = Levels_[startLevel];
         RCP<Smoother> preSmoo = Fine->GetPreSmoother();}
         RCP<Smoother> postSmoo = Fine->GetPostSmoother();

         //If on the coarse level, do either smoothing (if defined) or a direct solve.
         if (// check if on coarse level//) {

         //intermediate level
         } else {
           Coarse = Levels_[startLevel+];
           if (preSmoo != Teuchos::null)
             preSmoo->Apply(solution, rhs, InitialGuessIsZero);

           RCP<MultiVector> residual;
           RCP<MultiVector> workVector;
           //TODO FINISH ME 

         }

       */

     } //Iterate()

   //@}
    
}; //class Hierarchy

template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &hierarchy) {
  os << "Printing Hierarchy object" << std::endl;
  typename std::vector< Teuchos::RCP<Level<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> > >::const_iterator i;
  for (i=hierarchy.Levels_.begin(); i != hierarchy.Levels_.end(); ++i)
    os << *(*i) << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_HIERARCHY_SHORT

#endif //ifndef MUELU_HIERARCHY_HPP
