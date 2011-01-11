#ifndef stk_percept_ParallelUtilDef_hpp
#define stk_percept_ParallelUtilDef_hpp


namespace stk { 
  namespace percept { 


    //========================================================================================================================
    // implementation
    namespace {

#if defined( STK_HAS_MPI )

      template <typename T>
      struct stk_percept_Datatype;

      template <>
      struct stk_percept_Datatype<char>
      {
        static MPI_Datatype type() {
          return MPI_CHAR;
        }
      };

      template <>
      struct stk_percept_Datatype<signed char>
      {
        static MPI_Datatype type() {
          return MPI_CHAR;
        }
      };

      template <>
      struct stk_percept_Datatype<unsigned char>
      {
        static MPI_Datatype type() {
          return MPI_BYTE;
        }
      };

      template <>
      struct stk_percept_Datatype<int>
      {
        static MPI_Datatype type() {
          return MPI_INT;
        }
      };

      template <>
      struct stk_percept_Datatype<unsigned int>
      {
        static MPI_Datatype type() {
          return MPI_UNSIGNED;
        }
      };

      template <>
      struct stk_percept_Datatype<short>
      {
        static MPI_Datatype type() {
          return MPI_SHORT;
        }
      };

      template <>
      struct stk_percept_Datatype<unsigned short>
      {
        static MPI_Datatype type() {
          return MPI_UNSIGNED_SHORT;
        }
      };

      template <>
      struct stk_percept_Datatype<long>
      {
        static MPI_Datatype type() {
          return MPI_LONG;
        }
      };

      template <>
      struct stk_percept_Datatype<unsigned long>
      {
        static MPI_Datatype type() {
          return MPI_UNSIGNED_LONG;
        }
      };

      template <>
      struct stk_percept_Datatype<float>
      {
        static MPI_Datatype type() {
          return MPI_FLOAT;
        }
      };

      template <>
      struct stk_percept_Datatype<double>
      {
        static MPI_Datatype type() {
          return MPI_DOUBLE;
        }
      };


      //========================================================================================================================
      extern "C" {
        typedef void (*stk_percept_ParallelReduceOp)
          (void * inv, void * outv, int *, MPI_Datatype *);
      }

      template < typename T >
      struct stk_percept_reduce_min_lex
      {
        static void void_op(void * inv, void * outv, int *len, MPI_Datatype *)
        {
          T* inv_T = reinterpret_cast<T*>(inv);
          T* outv_T = reinterpret_cast<T*>(outv);

          bool is_less = false;
          for (int iv = 0; iv < *len; iv++)
            {
              if (inv_T[iv] < outv_T[iv])
                {
                  is_less = true;
                  break;
                }
              else if (inv_T[iv] == outv_T[iv])
                {
                  // continue
                }
              else
                {
                  is_less = false;
                  break;
                }
            }
          if (is_less)
            {
              std::copy(inv_T, inv_T+(*len), outv_T);
            }
        }
      };

#endif

      template<class T>
      inline
      void stk_percept_global_lex_min(stk::ParallelMachine comm,  int n , T local_min[] , T global_min[] )
      {
#if defined( STK_HAS_MPI )
        ParallelReduceOp p_op = reinterpret_cast<ParallelReduceOp>(&    stk_percept_reduce_min_lex<double>::void_op );

        if ( n < 1 ) return;
        MPI_Op mpi_op = MPI_OP_NULL ;

        MPI_Op_create(p_op, 0, & mpi_op);

        if ( MPI_SUCCESS != MPI_Allreduce( local_min, global_min, (int) n,
                                           stk_percept_Datatype<T>::type(),
                                           mpi_op,
                                           comm ) )  
          {
            throw std::runtime_error("MPI_Allreduce failed");
          }

        MPI_Op_free(& mpi_op);

#else
        for ( register int i = 0 ; i < n ; ++i )
          global_min[i] = local_min[i] ;
#endif
      }


    } // empty namespace


  } // percept
}// stk

#endif
