#ifndef Panzer_STK_ExodusReaderFactory_hpp__
#define Panzer_STK_ExodusReaderFactory_hpp__

#include <string>

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_MeshFactory.hpp"

#ifdef HAVE_IOSS

namespace panzer_stk {

class STK_Interface;

/** Concrete mesh factory instantiation. This reads
  * a mesh from an exodus file and builds a STK_Interface object
  */
class STK_ExodusReaderFactory : public STK_MeshFactory {
public:

   STK_ExodusReaderFactory();

   STK_ExodusReaderFactory(const std::string & fileName);

   /** Construct a STK_Inteface object described
     * by this factory.
     *
     * \param[in] parallelMach Descriptor for machine to build this mesh on.
     *
     * \returns Pointer to <code>STK_Interface</code> object with 
     *          <code>isModifiable()==false</code>.
     */ 
   virtual Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   /** This builds all the meta data of the mesh. Does not call metaData->commit.
     * Allows user to add solution fields and other pieces. The mesh can be "completed"
     * by calling <code>completeMeshConstruction</code>.
     */
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const 
   { return Teuchos::null; }

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const
   { return; }

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! Get file name mean is read from.
   const std::string & getFileName() const
   { return fileName_; }

protected:

   void registerElementBlocks(STK_Interface & mesh) const;
   void registerSidesets(STK_Interface & mesh) const;

   std::string fileName_;
};

}

#endif
#endif
