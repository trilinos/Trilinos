                     Basker (Version 0.1)
                    =====================

<Copyright> 
Basker: Threaded Sparse LU package for matrices with low fill-in
Copyright 2015 Sandia Corporation
Contact Joshua Dennis Booth (jdbooth@sandia.gov)

<Introduction to Basker>
Basker is a multithreaded LU factorization.  
It depends on SuiteSparse ordering, Scotch 6.0,  Kokkos, and C++11.
Currently, it is unknown how Basker will behave on a random system and with HWLOC.
Tested systems have been: Intel Sandy, Intel Phi, and IBM Power8.
This version does not contain incomplete factorization interfaces

<To Use>
As a package WITHOUT Amesos2.  
Enable ShyLUBasker with Kokkos OpeMP:
-DTrilinos_ENABLE_OpenMP:BOOL=ON
-DTrilinos_ENABLE_Kokkos:BOOL=ON
-DTrilinos_ENABLE_ShyLU:BOOL=ON
-DTrilinos_ENABLE_ShyLUBasker:BOOL=ON
See test for example.

As a package WITH Amesos2.
Do instruction with Amesos2, PLUS
-DCMAKE_CXX_FLAGS:STRING="-DSHYLU_NODEBASKER"

<Current Limitations>
Due to converting my experiemental code into a package that can easily be compiled in the Trilinos FrameWork, the following options where set.

Options Currently On:
BTF_LOWER, Multiple Thread Schur Lower Update, Multiple Threads Upper Update, Repeat Factor
Options Currently Off:
In Block Pivoting, BTF_UPPER


Basker File Structure:
**CNU = Currently not used.

Basker/README.txt          This document.  Overview of files and Basker
Basker/doc/                Personal working notes at this time
Basker/example/            Example drivers for basker native interface
Basker/src/                Source files
   basker.cpp              Main CPP file declaring the basker class
   basker_decl.hpp         Declaration of the basker class
   basker_def.hpp          Definitions of the public interfaces of the basker class
   basker_order.hpp        Definitions of ordering interfaces
   basker_order_amd.hpp    Definitions of AMD based ordering interfaces
   basker_order_btf.hpp    Definitions of BTF based ordering interfaces
   basker_order_hund.hpp   CNU
   basker_order_match.hpp  Definition of interface to mathcing ordering
   basker_order_scotch.hpp Definition of interface to ND from Scotch
   basker_sfactor.hpp      Definition of interface to symbolic factor
   basker_nfactor_blk.hpp  Definitions of private numerical factorization (ND Blks)
   basker_nfactor_col.hpp  Definitions of private numerical factorization (Sep Blks)
   basker_nfactor_col2.hpp Definitions of private numerical factorization (Sep Blks)
   basker_solve_rhs.hpp    Definitions of private solve phase
   basker_spbls.hpp        CNU
   basker_blas.hpp         CNU
   basker_stats.hpp        Statical measurement class (Not fully finished)
   basker_structs.hpp      Declaration of structures used by basker
   basker_thread.hpp       Definition of functions used to handle P2P syncs
   basker_tree.hpp         Definition of function for e/sep tree
   basker_types.hpp        Definition of MACROS for Kokkos
   basker_util.hpp         Random functions for i/o, initalization, and typde conversions
   basker_vector.hpp       Definitions of function for inital vector
   mwm2.hpp                A bottle-neck weighted card. matching


Compiler Option Flags:
NA

    


    
