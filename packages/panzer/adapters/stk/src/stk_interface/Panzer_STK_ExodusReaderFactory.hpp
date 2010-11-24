#ifndef Panzer_STK_ExodusReaderFactory_hpp__
#define Panzer_STK_ExodusReaderFactory_hpp__

#include <string>

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_MeshFactory.hpp"

#include <stk_mesh/fem/TopologicalMetaData.hpp>

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

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! Get file name mean is read from.
   const std::string & getFileName() const
   { return fileName_; }

protected:

   void registerElementBlocks(STK_Interface & mesh,const stk::mesh::TopologicalMetaData & md) const;
   void registerSidesets(STK_Interface & mesh,const stk::mesh::TopologicalMetaData & md) const;

   std::string fileName_;
};

}

#endif
#endif
