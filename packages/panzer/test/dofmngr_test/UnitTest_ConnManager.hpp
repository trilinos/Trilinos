#ifndef __UnitTest_ConnManager_hpp__
#define __UnitTest_ConnManager_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"

#ifdef HAVE_MPI
   #include "Teuchos_DefaultMpiComm.hpp"
   #include "mpi.h"
#else
   #include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace panzer {
namespace unit_test {

// simple debugging call back
class ConnCallback { public: virtual void buildConnectivity(const FieldPattern & fp) = 0; };

/** This a dummy connection manager used for testing. It has three element blocks
  *    "block_0", "block_1", and "block_2"
  * each with two elements.  Elements in each block reside on different processors.
  */
class ConnManager : public virtual panzer::ConnManager<short,int> {
public:
   typedef short LocalOrdinal;
   typedef int GlobalOrdinal;

   ConnManager(int rank,int procCount);

   ~ConnManager() {}

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const FieldPattern & fp);

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const;

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const;

   /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
   virtual std::string getBlockId(LocalOrdinal localElmtId) const;

   /** How many element blocks in this mesh?
     */
   virtual std::size_t numElementBlocks() const;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockID) const;

   void setBuildConnectivityCallback(const Teuchos::RCP<ConnCallback> & callback)
   { callback_ = callback; }

private:
   int procRank_;
   
   Teuchos::RCP<ConnCallback> callback_;
   std::map<std::string,std::vector<short> > elements_; // local element IDs
   std::vector<std::vector<int> > connectivity_;
};

} // end unit test
} // end panzer

#endif
