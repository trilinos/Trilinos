
#include <stk_adapt/RefinementInfoByType.hpp>
#include <stk_util/util/PrintTable.hpp>

namespace stk {
  namespace adapt {

    /** iRefinePass is used to predict number of new elements in query_only mode of UniformRefiner.
     *  Pass in 0 for iRefinePass unless you are using UniformRefiner in query_only mode, in which case
     *  pass in the current refinement pass
     */
    void RefinementInfoByType::printTable(std::ostream& os, std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass, bool printAll)
    {
      unsigned numOrigTot = 0;
      unsigned numNewTot = 0;
      for (unsigned irank = 0; irank < refinementInfoByType.size(); irank++)
        {
          int numOrig = refinementInfoByType[irank].m_numOrigElems;
          int numNew = refinementInfoByType[irank].m_numNewElems;
          if (numOrig)
            {
              double refFactor = ((double)numNew)/((double)numOrig);
              double refFactorNew = std::pow(refFactor, ((double)(iRefinePass+1) ));
              double refFactorOld = std::pow(refFactor, ((double)(iRefinePass) ));
              numNew = (int)((double)numOrig * refFactorNew);
              numOrig = (int)((double)numOrig * refFactorOld);
            }
          numOrigTot += numOrig;
          numNewTot += numNew;
        }

      //os << "Refinement Info> total original elements = " << numOrigTot << "\n";
      //os << "Refinement Info> total new elements = " << numNewTot << "\n";

      stk::PrintTable table;
      table.setTitle("Refinement Info\n");

      table << "|" << "Element Topology Type" << "|"
            << "Original Number of Elements" << "|" << "New Number of Elements" << "|"
            << stk::end_header;

      for (unsigned irank = 0; irank < refinementInfoByType.size(); irank++)
        {
          if (!printAll && refinementInfoByType[irank].m_numOrigElems == 0)
            continue;

          int numOrig = refinementInfoByType[irank].m_numOrigElems;
          int numNew = refinementInfoByType[irank].m_numNewElems;
          if (numOrig)
            {
              double refFactor = ((double)numNew)/((double)numOrig);
              double refFactorNew = std::pow(refFactor, ((double)(iRefinePass+1) ));
              double refFactorOld = std::pow(refFactor, ((double)(iRefinePass) ));
              numNew = (int)((double)numOrig * refFactorNew);
              numOrig = (int)((double)numOrig * refFactorOld);
            }

          table << "|" << refinementInfoByType[irank].m_topology.getName() << "|"
                << numOrig << "|"
                << numNew << "|"
                << stk::end_row;
        }

      table << "|" << "Totals" << "|"
            << numOrigTot << "|"
            << numNewTot << "|"
            << stk::end_row;


      os << "\n" << table;
    }

  }
}
