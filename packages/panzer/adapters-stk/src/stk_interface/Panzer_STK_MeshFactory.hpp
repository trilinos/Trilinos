// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Panzer_STK_MeshFactory_hpp__
#define Panzer_STK_MeshFactory_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include "Panzer_STK_PeriodicBC_Parser.hpp"

namespace panzer_stk {

class STK_Interface;

/** Pure virtual interface that constructs a 
  * STK_Mesh interface object.
  */
class STK_MeshFactory : public Teuchos::ParameterListAcceptorDefaultBase {
public:
   STK_MeshFactory() : enableRebalance_(false) {}

   /** Construct a STK_Inteface object described
     * by this factory.
     *
     * \param[in] parallelMach Descriptor for machine to build this mesh on.
     *
     * \returns Pointer to <code>STK_Interface</code> object with 
     *          <code>isModifiable()==false</code>.
     */ 
   virtual Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const = 0;

   /** This builds all the meta data of the mesh. Does not call metaData->commit.
     * Allows user to add solution fields and other pieces. The mesh can be "completed"
     * by calling <code>completeMeshConstruction</code>.
     */
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const = 0;

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const = 0;

   /** Parse the periodic boundary condition parameter list and build a vector of periodic boundary
     * conditions (a convenience function)
     */
   static void parsePeriodicBCList(const Teuchos::RCP<Teuchos::ParameterList> & pl,
                                   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & periodicBC,
                                   bool & useBBoxSearch)
   {
      panzer_stk::PeriodicBC_Parser parser;
      parser.setParameterList(pl);
      periodicBC = parser.getMatchers();
      useBBoxSearch = parser.useBoundingBoxSearch();
   }

   void enableRebalance(bool enable,const Teuchos::RCP<const Teuchos::ParameterList> & rebalanceList=Teuchos::null) 
   { enableRebalance_ = enable; 
     rebalanceList_ = rebalanceList; }

   void rebalance(STK_Interface & mesh) const
   {
     if(rebalanceList_!=Teuchos::null) {
       // loop over user specified partitioning lists
       for(Teuchos::ParameterList::ConstIterator itr=rebalanceList_->begin();
           itr!=rebalanceList_->end();++itr) {

         const Teuchos::ParameterEntry & entry = rebalanceList_->entry(itr);
         TEUCHOS_TEST_FOR_EXCEPTION(!entry.isList(),std::runtime_error,
                                    "Rebalance list is incorrect:\n" << entry << "\nA Zoltan list formated with strings is expected.");

         // partition according to the list
         mesh.rebalance(Teuchos::getValue<Teuchos::ParameterList>(entry));

         // rebuild mesh internals
         mesh.buildLocalElementIDs();
       }
     }
     else if(enableRebalance_) {
       // do the default thing, once
       Teuchos::ParameterList emptyList;
       mesh.rebalance(emptyList);

       // rebuild mesh internals
       mesh.buildLocalElementIDs();
     }
   }

   double getMeshCoord(const int nx, const double deltaX, const double x0) const {
      double x = static_cast<double>(nx)*deltaX;
      double modX = std::abs(x);
      double modX0 = std::abs(x0);
      double val = x+x0;
      if ((x0*x < 0.0) && (std::abs(modX-modX0) < std::numeric_limits<double>::epsilon()*modX0)) val=0.0;
      return (val);
   }

protected:
   // vector of periodic boundary condition objects
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCVec_; 
   // flag indicating which periodic search algorithm to use (bounding box or direct search)
   bool useBBoxSearch_;

   // for managing rebalance
   bool enableRebalance_;
   Teuchos::RCP<const Teuchos::ParameterList> rebalanceList_;
};

}

#endif
