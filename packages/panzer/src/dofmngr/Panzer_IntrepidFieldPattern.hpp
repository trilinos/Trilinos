#ifndef __Panzer_IntrepidFieldPattern_hpp__
#define __Panzer_IntrepidFieldPattern_hpp__

#include "Panzer_FieldPattern.hpp"

// Trilinos includes
#include "Intrepid_Basis.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

/** This is a derived class that specializes based
  * on a single intrepid basis function.
  */
class IntrepidFieldPattern : public FieldPattern {
public:
   IntrepidFieldPattern(const Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > & intrepidBasis)
      : intrepidBasis_(intrepidBasis) {}

   virtual int getSubcellCount(int dim) const;
   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const;
   virtual int getDimension() const;

   virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const;

   // static functions for examining shards objects

   /** For a given sub cell find the set of sub cells at all dimensions contained
     * internally.  This is inclusive, so that (dim,subCell) will be in the set.
     * 
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] closure Set of sub cells associated with specified sub cell.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>closure</code>.
     *       Previous information will not be removed.
     */
   static void buildSubcellClosure(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                                   std::set<std::pair<unsigned,unsigned> > & closure);

   /** Search a cell topology for sub cells containing a specfic set of nodes.
     * This is a downward search (inclusive) from a user specified dimension.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] nodes    Nodes forming the super set
     * \param[in,out] subCells Specific sub cells containing the nodes.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>subCells</code>.
     *       Previous information will not be removed.
     */
   static void findContainedSubcells(const shards::CellTopology & cellTopo,unsigned dim,
                                     const std::vector<unsigned> & nodes,
                                     std::set<std::pair<unsigned,unsigned> > & subCells);

   /** Get the set of nodes making up the user specified sub cells.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] nodes Nodes associated with sub cell.
     */
   static void getSubcellNodes(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                               std::vector<unsigned> & nodes);

protected:
   Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepidBasis_;
   std::vector<int> empty_;
};

}

#endif 
