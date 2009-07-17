#ifndef __PB_InverseLibrary_hpp__
#define __PB_InverseLibrary_hpp__

#include "PB_InverseFactory.hpp"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

namespace PB {

class InverseLibrary {
public:
   InverseLibrary();

   //! add an unspecified inverse to the library
   void addInverse(const std::string & label,const Teuchos::ParameterList & pl);

   //! Add a Stratimikos solver with a label to the library
   void addStratSolver(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   //! Add a Stratimikos preconditioner with a label to the library
   void addStratPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   //! Add a PB preconditioner to the library with a label
   void addBlockPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   /** Get the fully constructed parameter list for a particular label
     *
     * \param[in] label Name used for the desired solver.
     *
     * \returns If the label is found in the library the corresponding parameter list
     *          is returned, otherwise <code>Teuchos::null</code> is returned.
     */
   Teuchos::RCP<const Teuchos::ParameterList> getParameterList(const std::string & label) const;

   //! Get the inverse factory associated with a particular label
   Teuchos::RCP<const InverseFactory> getInverseFactory(const std::string & label) const;

protected:

   // stratimikos type Inverse objects: mapping the label to a parameter list
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > stratSolver_;
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > stratPrecond_;
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > blockPrecond_;

   // vectors showing which string types are in Stratimikos
   std::vector<std::string> stratValidSolver_;
   std::vector<std::string> stratValidPrecond_;
   std::vector<std::string> blockValidPrecond_;
    
public:

   /** \brief Build an inverse library from a parameter list.
     * 
     * Build an inverse library from a parameter list. This will
     * contain all the labeled inverses specified.
     *
     * \param[in] pl Parameter list to build the library from
     * \param[in] useStratDefaults Also add the default parameters from Stratimikos
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromParameterList(const Teuchos::ParameterList & pl,bool useStratDefaults=true);

   /** \brief Build an inverse library from Stratimikos
     * 
     * Build an inverse library from Stratimkos. The labels
     * will just be the names in Stratimikos.
     *
     * \param[in] strat Stratimikos object to use
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromStratimikos(
         const Stratimikos::DefaultLinearSolverBuilder & strat=Stratimikos::DefaultLinearSolverBuilder());
};

} // end namespace PB

#endif
