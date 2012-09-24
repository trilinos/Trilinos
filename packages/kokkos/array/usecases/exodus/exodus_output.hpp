/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <string>
#include <vector>
#include <exodusII.h>

namespace KokkosArray {
namespace Exodus {

/** \brief  Create output mesh
 *
 *  'node_output_values_names' and 'elem_output_values_names'
 *  are the names for each node and element scalar value to be output.
 *  For example, output of a nodal displacement array ( #nodes , 3 )
 *  requires 3 names.
 */
template< class CoordinateArray , class ElementNodeArray >
int create_output_mesh( const std::string & file_path ,
                        const std::string & title ,
                        const CoordinateArray  & node_coord ,
                        const ElementNodeArray & elem_nodes ,
                        const std::vector<std::string> & node_output_values_names ,
                        const std::vector<std::string> & elem_output_values_names );

/* \brief  Output time value for step.  Step = 1..N */

void output_mesh_time( const int exo_id , const int step , const double time );

/* \brief  Output nodal values for step.  Step = 1..N
 *
 *  Each value will be named as follows:
 *
 *  for ( k = 0 ; k < array.dimension(2) ; ++k ) {
 *    for ( j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
 *      node_output_values_names[ var_index++ ]
 *    }
 *  }
 */
template< class NodalArray >
void output_mesh_nodal_variable( const int exo_id ,
                                 const int step ,
                                 int var_index ,
                                 const NodalArray & array );

/* \brief  Output element values for step.  Step = 1..N 
 *
 *  Each value will be named as follows:
 *
 *  for ( k = 0 ; k < array.dimension(2) ; ++k ) {
 *    for ( j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
 *      elem_output_values_names[ var_index++ ]
 *    }
 *  }
 */
template< class ElementArray >
void output_mesh_element_variable( const int exo_id ,
                                   const int step ,
                                   int var_index ,
                                   const ElementArray & array );

void close_output_mesh( const int exo_id );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace {
inline void copy_string( char * dst , const char * src )
{ strncpy( dst , src , strlen( src ) ); }
}

void close_output_mesh( const int exo_id )
{ ex_close( exo_id ); }

template< class CoordinateArray , class ElementNodeArray >
int create_output_mesh( const std::string & file_path ,
                        const std::string & title ,
                        const CoordinateArray  & node_coord ,
                        const ElementNodeArray & elem_nodes ,
                        const std::vector<std::string> & node_output_values_names ,
                        const std::vector<std::string> & elem_output_values_names )
{
  typedef double value_type ;

  const const * const char_ptr_null = 0 ;

  const int input_word_size = sizeof(value_type);
  const int output_word_size = input_word_size ;
  const int num_dim = 3 ;
  const int num_nodes_global = node_coord.dimension(0);
  const int num_elems_global = elem_nodes.dimension(0);
  const int num_elem_blk = 1 ;
  const int num_node_sets = 0 ;
  const int num_side_sets = 0 ;

  const int exo_id = ex_create( file_path.c_str() , EX_CLOBBER , & input_word_size , & output_word_size );

  if ( m_exo_id < 0 ) throw std::runtime_error( std::string("ex_create") );
  
  int exo_error ;

  exo_error = ex_put_init( exo_id , title.c_str() ,
                           num_dim , num_nodes_global , num_elems_global ,
                           num_elem_blk , num_node_sets , num_side_sets );

  if ( exo_error ) throw std::runtime_error( std::string("ex_put_init") );
  
  // Put the model nodal coordinate names:
  const char coord_name_0[ MAX_STR_LENGTH ] = "model_coordinates_x" ;
  const char coord_name_1[ MAX_STR_LENGTH ] = "model_coordinates_y" ;
  const char coord_name_2[ MAX_STR_LENGTH ] = "model_coordinates_z" ;
    
  const char * coord_names[3] = { coord_name_0 , coord_name_1 , coord_name_2 };
    
  exo_error = ex_put_coord_names( exo_id , coord_names );

  if ( exo_error ) throw std::runtime_error( std::string("ex_put_coord_names") );

  // Put element blocks' description:

  char * const char_ptr_null = NULL ;

  std::vector<char> tmp_type( num_elem_blk * MAX_STR_LENGTH );
  std::vector<char> tmp_name( num_elem_blk * MAX_STR_LENGTH );
    
  std::vector<int>    elem_blk_id(   num_elem_blk , 0 );
  std::vector<char *> elem_blk_name( num_elem_blk , char_ptr_null );
  std::vector<char *> elem_blk_type( num_elem_blk , char_ptr_null );
  std::vector<int>    num_nodes_per_elem( num_elem_blk , 0 );
  std::vector<int>    num_attr( num_elem_blk , 0 );

  elem_blk_id[0]        = 1 ;
  num_nodes_per_elem[0] = 8 ;
  num_attr[0]           = 0 ;
  elem_blk_type[0]      = & tmp_type[ 0 * MAX_STR_LENGTH ] ;
  elem_blk_name[0]      = & tmp_name[ 0 * MAX_STR_LENGTH ] ;

  copy_string( elem_blk_name[0] , "ElementBlock_1" );
  copy_string( elem_blk_type[0] , "HEX_8" );

  const int will_define_identifiers = false ; // will not do this

  exo_error = ex_put_concat_elem_block( exo_id ,
                                        & elem_blk_id[0] ,
                                        & elem_blk_type[0] ,
                                          num_elems_global ,
                                        & num_nodes_per_elem[0] ,
                                        & num_attr[0] ,
                                        will_define_identifiers );

  if ( exo_error ) throw std::runtime_error( std::string("ex_put_concat_elem_block") );

  exo_error = ex_put_names( m_exo_id, EX_ELEM_BLOCK, & elem_blk_name[0] );

  if ( exo_error ) throw std::runtime_error( std::string("ex_put_names") );

  // Put results variables descriptions

  const int num_var_node = node_output_scalars_name.size();
  const int num_var_elem = elem_output_scalars_name.size();
  const int num_var_global = 0 ;
  const int num_var_side_set = 0 ;
  const int num_var_node_set = 0 ;

  std::vector<int> exists_elem_var( num_elem_blk * num_output_scalars_per_elem , 1 );

  exo_error = ex_put_all_var_param( exo_id ,
                                    num_var_global ,
                                    num_var_node ,
                                    num_var_elem ,
                                    & exist_elem_var[0] ,
                                    num_var_node_set ,
                                    NULL ,
                                    num_var_side_set ,
                                    NULL );

  if ( exo_error ) throw std::runtime_error( std::string("ex_put_all_var_param") );

  if ( num_var_node ) {
    std::vector<char *> ptr_var_names( num_var_node );

    for ( int i = 0 ; i < num_var_node ; ++i ) {
      ptr_var_names[i] = const_cast<char*>( node_output_values_names[i].c_str() );
    }

    exo_error = ex_put_var_names( exo_id , "n" , num_var_node , & ptr_var_names[0] );

    if ( exo_error ) throw std::runtime_error( std::string("ex_put_var_names") );
  }

  if ( num_var_elem ) { 
    std::vector<char *> ptr_var_names( num_var_elem );

    for ( unsigned i = 0 ; i < name_elem_var.size() ; ++i ) {
      ptr_var_names[i] = const_cast<char*>( elem_output_values_names[i].c_str() );
    }

    exo_error = ex_put_var_names( m_exo_id , "e" , num_var_elem , & ptr_var_names[0] );

    if ( exo_error ) throw std::runtime_error( std::string("ex_put_var_names") );
  }

  // Write nodal coordinates

  {
    std::vector<value_type> x( num_nodes_global ), y( num_nodes_global ), z( num_nodes_global );

    for ( int i = 0 ; i < num_nodes_global ; ++i ) {
      x[i] = node_coord(i,0);
      y[i] = node_coord(i,1);
      z[i] = node_coord(i,2);
    }

    exo_error = ex_put_coord( exo_id , & x[0] , & y[0] , & z[0] );

    if ( exo_error ) throw std::runtime_error( std::string("ex_put_coord") );
  }

  // Write element connectivity for the one and only element block

  {
    std::vector<int> connect( num_elems_global * num_nodes_per_elem[0] );

    int k = 0 ;
    for ( int i = 0 ; i < num_elems_global ; ++i ) {
      for ( int j = 0 ; j < num_nodes_per_elem[0] ; ++j , ++k ) {
        connect[k] = elem_nodes(i,j);
      }
    }

    exo_error = ex_put_elem_conn( exo_id , 0 , & connect[0] );

    if ( exo_error ) throw std::runtime_error( std::string("ex_put_elem_conn") );
  }

  return exo_id ;
}

//----------------------------------------------------------------------------

void output_mesh_time( const int exo_id , const int step , const double time )
{
  int exo_error = ex_put_time( exo_id , step , time );
  if ( exo_error ) throw std::runtime_error( std::string("ex_put_time") );
}

//----------------------------------------------------------------------------

template< class NodalArray >
void output_mesh_nodal_variable( const int exo_id ,
                                 const int step ,
                                 int var_index ,
                                 const NodalArray & array )
{
  const int length = array.dimension(0);

  std::vector<double> values( length );

  int exo_error ;

  switch( array.rank() ) {
  case 1 :
    for ( int i = 0 ; i < length ; ++i ) { values[i] = array(i); }
    exo_error = ex_put_nodal_var( exo_id , step , var_index , length , & values[0] );
    if ( exo_error ) throw std::runtime_error( std::string("ex_put_nodal_var") );
    break ;

  case 2 :
    for ( int j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
      for ( int i = 0 ; i < length ; ++i ) { values[i] = array(i,j); }
      exo_error = ex_put_nodal_var( exo_id , step , var_index , length , & values[0] );
      if ( exo_error ) throw std::runtime_error( std::string("ex_put_nodal_var") );
    }
    break ;

  case 3 :
    for ( int k = 0 ; k < array.dimension(2) ; ++k ) {
      for ( int j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
        for ( int i = 0 ; i < length ; ++i ) { values[i] = array(i,j,k); }
        exo_error = ex_put_nodal_var( exo_id , step , var_index , length , & values[0] );
        if ( exo_error ) throw std::runtime_error( std::string("ex_put_nodal_var") );
      }
    }
    break ;
  }
}

template< class ElementArray >
void output_mesh_element_variable( const int exo_id ,
                                   const int step ,
                                   int var_index ,
                                   const ElementArray & array )
{
  const int elem_blk_id = 1 ;
  const int length = array.dimension(0);

  std::vector<double> values( num_elems );

  int exo_error ;

  switch( array.rank() ) {
  case 1 :
    for ( int i = 0 ; i < num_elems ; ++i ) { values[i] = array(i); }
    exo_error = ex_put_elem_var( exo_id , step , var_index , elem_blk_id , length , & values[0] );
    if ( exo_error ) throw std::runtime_error( std::string("ex_put_elem_var") );
    break ;

  case 2 :
    for ( int j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
      for ( int i = 0 ; i < length ; ++i ) { values[i] = array(i,j); }
      exo_error = ex_put_elem_var( exo_id , step , var_index , elem_blk_id , length , & values[0] );
      if ( exo_error ) throw std::runtime_error( std::string("ex_put_elem_var") );
    }
    break ;

  case 3 :
    for ( int k = 0 ; k < array.dimension(2) ; ++k ) {
      for ( int j = 0 ; j < array.dimension(1) ; ++j , ++var_index ) {
        for ( int i = 0 ; i < length ; ++i ) { values[i] = array(i,j,k); }
        exo_error = ex_put_elem_var( exo_id , step , var_index , elem_blk_id , num_nodes , & values[0] );
        if ( exo_error ) throw std::runtime_error( std::string("ex_put_elem_var") );
      }
    }
    break ;
  }
}

//----------------------------------------------------------------------------

} // namespace Exodus 
} // namespace KokkosArray 





