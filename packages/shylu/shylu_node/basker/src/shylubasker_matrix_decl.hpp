#ifndef SHYLUBASKER_MATRIX_DECL_HPP
#define SHYLUBASKER_MATRIX_DECL_HPP

/*Basker Includes*/
#include "shylubasker_types.hpp"

/*System Includes*/
#include <limits>
#include <string>

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#else
#include <omp.h>
#endif

using std::string;

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  class BaskerMatrix
  {
 
  public:
    
    //Constructors and deconstructors
    BASKER_INLINE
    BaskerMatrix();
    BASKER_INLINE
    BaskerMatrix(string _label);
    BASKER_INLINE
    BaskerMatrix(Int _m, Int _n, Int _nnz, 
                 Int *col_ptr, Int *row_idx, Entry *val);
    BASKER_INLINE
    BaskerMatrix(string _label, Int _m, Int _n, Int _nnz, 
                 Int *col_ptr, Int *row_idx, Entry *val);
    BASKER_INLINE
    ~BaskerMatrix();

    /*  We might want to define this in the future.
    BASKER_INLINE
    BaskerMatrix<Int,Entry,Exe_Space>& operator= (const BaskerMatrix<Int,Entry,Exe_Space>&);
    */

    //init_matrix (want to change these to malloc_matrix)
    BASKER_INLINE
    void init_matrix(string _label, Int _m, Int _n, Int _nnz);

    BASKER_INLINE
    void init_matrix(string _label, Int _m, Int _n, Int _nnz,
                    Int *_col_ptr, Int *_row_idx, Entry *_val);
    BASKER_INLINE
    void init_matrix(string _label, Int _sr, Int _m, 
                    Int _sc, Int _n, Int _nnz);

    //finalize, used to delete any array structure that where created
    BASKER_INLINE
    void Finalize();

    //
    BASKER_INLINE
    int copy_values(Int _sr, Int _m, Int _sc, Int _n, Int _nnz,
		    Int *_col_ptr, Int *_row_idx, Entry *_val);
    BASKER_INLINE
    int copy_values(Int _m, Int _n, Int _nnz,
		    Int *_col_ptr, Int *_row_idx, Entry *_val);

    BASKER_INLINE
    void init_col();
    BASKER_INLINE
    void clean_col();
    BASKER_INLINE
    void convert2D(BASKER_MATRIX &M, 
		   BASKER_BOOL alloc, 
		   Int kid);

    
    //just set shape, do not init
    void set_shape(Int _sr, Int _m, 
		  Int _sc, Int _n);


    BASKER_INLINE
    int fill();

    BASKER_INLINE
    void init_inc_lvl();


    //****Deprecated*******
    BASKER_INLINE
    void malloc_perm(Int n);
    BASKER_INLINE
    void init_perm();

    //****Deprecated*****
    //malloc union_bit
    BASKER_INLINE
    void malloc_union_bit();
    BASKER_INLINE
    void init_union_bit();
    BASKER_INLINE
    void init_union_bit(Int kid);

    //helper functions
    void copy_vec(Int* ptr, Int size, INT_1DARRAY a);
    void copy_vec(Entry *ptr, Int size,  ENTRY_1DARRAY a);
    BASKER_INLINE
    void init_vectors(Int _m, Int _n, Int _nnz);

    //information
    BASKER_INLINE
    void info();
    BASKER_INLINE
    void level_info();
    BASKER_INLINE
    void print();


    //Note: These need to be reordered to make better use of 
    //Class size.

    string label;


    BASKER_BOOL v_fill;

    Int srow, scol; //start col (wrt global matrix, if a block)
    Int erow, ecol; //end col (wrt global matrix, if a block)
    Int ncol, nrow, nnz;
    Int mnnz; //malloc nnz
    
    INT_1DARRAY   col_ptr;
    INT_1DARRAY   row_idx;
    ENTRY_1DARRAY val;

    //**Deprecated***
    INT_1DARRAY   lpinv;

    //***Deprecated***
    BOOL_1DARRAY union_bit;
   
   
    //#ifdef BASKER_INC_LVL
    BASKER_BOOL   inc_lvl_flg;
    INT_1DARRAY   inc_lvl;
    //#endif

    #ifdef BASKER_2DL
    BASKER_BOOL   w_fill;
    ENTRY_1DARRAY ews;
    INT_1DARRAY   iws;
    Int           iws_size;
    Int           ews_size;
    Int           iws_mult;
    Int           ews_mult;
    Int           p_size;
    #endif
    
    Entry tpivot;

    //***Deprecated***
    //Remove..... will not be used in future ver
    static const Int max_idx = (Int) -1;
    
    //Used for end
    INT_1DARRAY pend;
    void init_pend();
    void clear_pend();

  };//end class BaskerMatrix

}//End namespace BaskerNS

#endif //End ifndef basker_matrix_decl_hpp
