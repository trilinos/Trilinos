// -*- c-file-style: "krlong" -*-
void Assembler::getGraph(int br, int bc,
                         Array<int>& graphData,
                         Array<int>& rowPtrs,
                         Array<int>& nnzPerRow) const
{
  int i;
  {
    for (unsigned int d=0; d<eqn_->numRegions(); d++)
      {
        while (iter != cells.end())
          {
            if (pairs.get() != 0)
              {
                for (int c=0; c<nCells; c++)
                  {
                    for (int t=0; t<nt; t++)
                      {
                        for (unsigned int uit=0; uit<unksForTests[t].size(); uit++)
                          { 
                            for (int n=0; n<nTestNodes; n++)
                              {
                                int row 
                                  = testDOFs[(c*nTestFuncs + testFuncIndex)*nTestNodes + n];
                                if ( row < lowestRow_[br] || row >= highestRow
                                     || (*(isBCRow_[br]))[row-lowestRow_[br]]
                                     )
                                  continue;
                                Set<int>& colSet = tmpGraph[row-lowestRow_[br]];
                                for (int m=0; m<nUnkNodes; m++)
                                  {
                                    int col
                                      = unkDOFs[(c*nUnkFuncs + unkFuncIndex)*nUnkNodes + m];
                                    colSet.put(col);
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }
  }
}
