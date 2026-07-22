// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Panzer_STK_ExodusReaderFactory_hpp__
#define Panzer_STK_ExodusReaderFactory_hpp__

#include <string>

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_MeshFactory.hpp"

#ifdef PANZER_HAVE_IOSS

#include <stk_io/StkMeshIoBroker.hpp>

namespace panzer_stk {

class STK_Interface;

/** External function to return the dimension of an Exodus or Pamgen
 * mesh. This uses a quick temporary read of the meta data from input
 * file/mesh and is intended if you need the mesh dimension before the
 * MeshFactory is used to build the uncommitted mesh. Once
 * buildUncommittedMesh() is called, you can query the dimension from
 * the MetaData in the STK_Interface object (will be faster than
 * creating mesh metadata done in this function).
 *
 * \param[in] meshStr Filename containing the mesh string, or the mesh string itself.
 * \param[in] parallelMach Descriptor for machine to build this mesh on.
 *
 * \returns Integer indicating the spatial dimension of the mesh.
 */
  int getMeshDimension(const std::string & meshStr,stk::ParallelMachine parallelMach, const std::string & typeStr = "Exodus");

  std::string fileTypeToIOSSType(const std::string & fileType);

/** Concrete mesh factory instantiation. This reads
  * a mesh from an exodus file and builds a STK_Interface object.
  *
  * Also, if a nonzero restart index (the Exodus indices are 1 based) is
  * specified then this will set the initial state time in the STK_Interface.
  * However, as prescribed by that interface the currentStateTime is only
  * set by the writeToExodus call, thus the currentStateTime of the created
  * STK_Interface object will be zero. It is up to the user to rectify this
  * when calling writeToExodus.
  */
class STK_ExodusReaderFactory : public STK_MeshFactory {
public:

   STK_ExodusReaderFactory();


  /** \brief Ctor
   *
   * \param[in] fileName Name of the input file.
   * \param[in] restartIndex Index used for restarts.
   */
  STK_ExodusReaderFactory(const std::string & fileName, const int restartIndex=0);

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
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor. Must be called if the empty ctor is used to construct this object.
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! Get file name mean is read from.
   const std::string & getFileName() const
   { return fileName_; }

protected:

   void registerElementBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const;
   void registerSidesets(STK_Interface & mesh) const;
   void registerNodesets(STK_Interface & mesh) const;
   void registerEdgeBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const;
   void registerFaceBlocks(STK_Interface & mesh,stk::io::StkMeshIoBroker & meshData) const;

   void addEdgeBlocks(STK_Interface & mesh) const;
   void addFaceBlocks(STK_Interface & mesh) const;

   std::string mkBlockName(std::string base, std::string topo_name) const;
   void createUniqueEdgeTopologyMap(STK_Interface & mesh, const stk::mesh::Part *elemBlockPart) const;
   void createUniqueFaceTopologyMap(STK_Interface & mesh, const stk::mesh::Part *elemBlockPart) const;

   void buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const;

   bool doPerceptRefinement() const;

   std::string fileName_;
   std::string fileType_;
   int restartIndex_;

   /* The ExodusReaderFactory creates one edge/face block for each
    * unique edge/face topology in the mesh.  There are a few
    * situations where it's desirable to have a list of unique
    * topologies for each element block.  Instead of creating it
    * on the fly, they are created and saved when the element
    * blocks are added to the STK_Interface in
    * registerElementBlocks().
    */
   mutable std::map<std::string,std::vector<stk::topology>> elemBlockUniqueEdgeTopologies_;
   mutable std::map<std::string,std::vector<stk::topology>> elemBlockUniqueFaceTopologies_;

private:

  //! Did the user request mesh scaling
  bool userMeshScaling_;

  //! Did the user request to keep the Percept data
  bool keepPerceptData_;

  //! Did the user request to keep the Percept parent element data
  bool keepPerceptParentElements_;

  //! The type of mesh rebalancing to be performed after creation
  std::string rebalancing_;

  //! If requested, scale the input mesh by this factor
  double meshScaleFactor_;

  //! Number of levels of inline uniform mesh refinement to be applied to exodus mesh
  int levelsOfRefinement_;

  //! Did the user request to create missing edge blocks
  bool createEdgeBlocks_;

  //! Did the user request to create missing face blocks
  bool createFaceBlocks_;

  std::string geometryName_;
};

}

#endif
#endif
