
#ifndef stk_util_util_IndexList_h
#define stk_util_util_IndexList_h

namespace stk {

/** \brief Compile-time list of indices.
 *  \ingroup util_module
 */
template< unsigned  I0 = 0 , unsigned  I1 = 0 ,
          unsigned  I2 = 0 , unsigned  I3 = 0 ,
          unsigned  I4 = 0 , unsigned  I5 = 0 ,
          unsigned  I6 = 0 , unsigned  I7 = 0 ,
          unsigned  I8 = 0 , unsigned  I9 = 0 ,
          unsigned I10 = 0 , unsigned I11 = 0 ,
          unsigned I12 = 0 , unsigned I13 = 0 ,
          unsigned I14 = 0 , unsigned I15 = 0 ,
          unsigned I16 = 0 , unsigned I17 = 0 ,
          unsigned I18 = 0 , unsigned I19 = 0 ,
          unsigned I20 = 0 , unsigned I21 = 0 ,
          unsigned I22 = 0 , unsigned I23 = 0 ,
          unsigned I24 = 0 , unsigned I25 = 0 ,
          unsigned I26 = 0 , unsigned I27 = 0 ,
          unsigned I28 = 0 , unsigned I29 = 0 ,
          unsigned I30 = 0 , unsigned I31 = 0 >
struct IndexList {};

/** \brief Access member of compile-time list of indices. <br>
 *         Defines <b> enum { value = index_at_J }; </b> 
 *  \ingroup util_module
 */
template< class List , unsigned J > struct IndexListAt {};

#define INDEX_LIST_AT_SPECIALIZATION( J , K )	\
  template< unsigned  I0 , unsigned  I1 ,	\
            unsigned  I2 , unsigned  I3 ,	\
            unsigned  I4 , unsigned  I5 ,	\
            unsigned  I6 , unsigned  I7 ,	\
            unsigned  I8 , unsigned  I9 ,	\
            unsigned I10 , unsigned I11 ,	\
            unsigned I12 , unsigned I13 ,	\
            unsigned I14 , unsigned I15 ,	\
            unsigned I16 , unsigned I17 ,	\
            unsigned I18 , unsigned I19 ,	\
            unsigned I20 , unsigned I21 ,	\
            unsigned I22 , unsigned I23 ,	\
            unsigned I24 , unsigned I25 ,	\
            unsigned I26 , unsigned I27 ,	\
            unsigned I28 , unsigned I29 ,	\
            unsigned I30 , unsigned I31 >	\
struct IndexListAt<	\
  IndexList< I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,	\
             I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,	\
            I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,	\
            I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > , J >	\
{ enum { value = K }; };

INDEX_LIST_AT_SPECIALIZATION(  0 ,  I0 )
INDEX_LIST_AT_SPECIALIZATION(  1 ,  I1 )
INDEX_LIST_AT_SPECIALIZATION(  2 ,  I2 )
INDEX_LIST_AT_SPECIALIZATION(  3 ,  I3 )
INDEX_LIST_AT_SPECIALIZATION(  4 ,  I4 )
INDEX_LIST_AT_SPECIALIZATION(  5 ,  I5 )
INDEX_LIST_AT_SPECIALIZATION(  6 ,  I6 )
INDEX_LIST_AT_SPECIALIZATION(  7 ,  I7 )
INDEX_LIST_AT_SPECIALIZATION(  8 ,  I8 )
INDEX_LIST_AT_SPECIALIZATION(  9 ,  I9 )
INDEX_LIST_AT_SPECIALIZATION( 10 , I10 )
INDEX_LIST_AT_SPECIALIZATION( 11 , I11 )
INDEX_LIST_AT_SPECIALIZATION( 12 , I12 )
INDEX_LIST_AT_SPECIALIZATION( 13 , I13 )
INDEX_LIST_AT_SPECIALIZATION( 14 , I14 )
INDEX_LIST_AT_SPECIALIZATION( 15 , I15 )
INDEX_LIST_AT_SPECIALIZATION( 16 , I16 )
INDEX_LIST_AT_SPECIALIZATION( 17 , I17 )
INDEX_LIST_AT_SPECIALIZATION( 18 , I18 )
INDEX_LIST_AT_SPECIALIZATION( 19 , I19 )
INDEX_LIST_AT_SPECIALIZATION( 20 , I20 )
INDEX_LIST_AT_SPECIALIZATION( 21 , I21 )
INDEX_LIST_AT_SPECIALIZATION( 22 , I22 )
INDEX_LIST_AT_SPECIALIZATION( 23 , I23 )
INDEX_LIST_AT_SPECIALIZATION( 24 , I24 )
INDEX_LIST_AT_SPECIALIZATION( 25 , I25 )
INDEX_LIST_AT_SPECIALIZATION( 26 , I26 )
INDEX_LIST_AT_SPECIALIZATION( 27 , I27 )
INDEX_LIST_AT_SPECIALIZATION( 28 , I28 )
INDEX_LIST_AT_SPECIALIZATION( 29 , I29 )
INDEX_LIST_AT_SPECIALIZATION( 30 , I30 )
INDEX_LIST_AT_SPECIALIZATION( 31 , I31 )

#undef INDEX_LIST_AT_SPECIALIZATION

} // namespace stk

#endif // stk_util_util_IndexList_h

