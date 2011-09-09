
#include <stk_adapt/RefinementInfoByType.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_percept/stk_mesh.hpp>

namespace stk {
  namespace adapt {

    /** iRefinePass is used to predict number of new elements in query_only mode of UniformRefiner.
     *  Pass in 0 for iRefinePass unless you are using UniformRefiner in query_only mode, in which case
     *  pass in the current refinement pass
     */
    void RefinementInfoByType::printTable(std::ostream& os, std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass, bool printAll)
    {
      RefinementInfoCount numOrigTot = 0;
      RefinementInfoCount numNewTot = 0;
      for (unsigned irank = 0; irank < refinementInfoByType.size(); irank++)
        {
          RefinementInfoCount numOrig = refinementInfoByType[irank].m_numOrigElems;
          RefinementInfoCount numNew = refinementInfoByType[irank].m_numNewElems;
          if (numOrig)
            {
              double refFactor = ((double)numNew)/((double)numOrig);
              double refFactorNew = std::pow(refFactor, ((double)(iRefinePass+1) ));
              double refFactorOld = std::pow(refFactor, ((double)(iRefinePass) ));
              numNew = (RefinementInfoCount)((double)numOrig * refFactorNew);
              numOrig = (RefinementInfoCount)((double)numOrig * refFactorOld);
            }
          numOrigTot += numOrig;
          numNewTot += numNew;
        }

      RefinementInfoCount numOrigNodes = refinementInfoByType[0].m_numOrigNodes;
      RefinementInfoCount numNewNodes = refinementInfoByType[0].m_numNewNodes;

      //os << "Refinement Info> total original elements = " << numOrigTot << "\n";
      //os << "Refinement Info> total new elements = " << numNewTot << "\n";

      stk::PrintTable table;
      table.setTitle("Refinement Info\n");
      //table.setAutoEndCol(false);

      table << "|" << "                     " <<  "|" << justify(PrintTable::Cell::CENTER) 
        << "Original" << stk::end_col << "     " << "|" << "New     " << stk::end_col << "     " << "|" << stk::end_header;
      //    << stk::span << "Original"  << "|" << stk::span << "New     "  << "|" << stk::end_header;
      table << justify(PrintTable::Cell::LEFT) ;
      table << "|" << "Element Topology Type" <<  "|" << "Elements" << stk::end_col << "Nodes" << "|" << "Elements" << stk::end_col << "Nodes" << "|" << stk::end_header;

      for (unsigned irank = 0; irank < refinementInfoByType.size(); irank++)
        {
          if (!printAll && refinementInfoByType[irank].m_numOrigElems == 0)
            continue;

          RefinementInfoCount numOrig = refinementInfoByType[irank].m_numOrigElems;
          RefinementInfoCount numNew = refinementInfoByType[irank].m_numNewElems;
          if (numOrig)
            {
              double refFactor = ((double)numNew)/((double)numOrig);
              double refFactorNew = std::pow(refFactor, ((double)(iRefinePass+1) ));
              double refFactorOld = std::pow(refFactor, ((double)(iRefinePass) ));
              numNew = (RefinementInfoCount)((double)numOrig * refFactorNew);
              numOrig = (RefinementInfoCount)((double)numOrig * refFactorOld);
            }

          table << "|" << refinementInfoByType[irank].m_topology.getName() << "|"
                << numOrig << stk::end_col << " " << "|"
                << numNew  << stk::end_col << " " << "|"
            //                 << numOrig << stk::end_col << numOrigNodes << "|"
            //                 << numNew  << stk::end_col << numNewNodes << "|"
                << stk::end_row;
        }

      table << "|" << "Totals" << "|"
            << numOrigTot << stk::end_col << numOrigNodes << "|"
            << numNewTot  << stk::end_col << numNewNodes << "|"
            << stk::end_row;

#if 0
      table << "|" << "Node Totals" << "|"
            << numOrigNodes << "|"
            << numNewNodes << "|"
            << stk::end_row;
#endif

      os << "\n" << table;
    }

    void RefinementInfoByType::countCurrentNodes(stk::percept::PerceptMesh& eMesh, std::vector< RefinementInfoByType >& refinementInfoByType)
    {
      mesh::Selector selector(eMesh.getFEM_meta_data()->locally_owned_part());
      std::vector<unsigned> count ;
      stk::mesh::count_entities( selector, *eMesh.getBulkData(), count );
      
      unsigned nnodes = count[0];

      stk::ParallelMachine pm = eMesh.getBulkData()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );

      for (unsigned i = 0; i < refinementInfoByType.size(); i++)
        {
          refinementInfoByType[i].m_numNewNodes = nnodes;
        }
    }

  }
}
