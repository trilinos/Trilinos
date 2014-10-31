#ifdef BASKER_DECL_HPP
#define BASKER_DECL_HPP


namespace basker{

  template<class Int, class Entry>
  Int basker
  (
   Int Ap [],    //column pointer
   Int Ai [],    //row index
   Entry Ax[],   //values
   Int anrow,    //number of rows in A
   Int ancol,    //number of columns in A
   Int ws[],     //Int workspace
   Entry X[],    //Entry workspace
   Int *LP,      //Column pointer
   Int **Li,     //row index array
   Entry **Lx,   //values 
   Int *Up, 
   Int **Ui,
   Entry **Ux, 
   Int *lnnz,
   Int *unnz,
   Int *pinv
   );


  template<class Int, class Entry>
  void basker_dfs
  (
   Int n,
   Int j,
   Int Li[],
   Int Lp[],
   Int color[],
   Int pattern[],
   Int *top,
   Int k,
   Int pinv[],
   Int stack[]
   );

}
#endif
