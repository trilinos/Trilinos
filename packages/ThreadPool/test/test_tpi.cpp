/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <TPI.hpp>

/*------------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

template<unsigned N> class TEST ;

template<unsigned N>
class TEST {
public:
  int m_flag[N] ;
  ~TEST() {}
  TEST();
  void flag( TPI::Work & );
  void verify();
private:
  TEST( const TEST & );
  TEST & operator = ( const TEST & );
};

template<unsigned N>
TEST<N>::TEST()
{
  for ( unsigned i = 0 ; i < N ; ++i ) { m_flag[i] = 0 ; }
}

template<unsigned N>
void TEST<N>::flag( TPI::Work & work )
{
  static const char method[] = "TEST::flag" ;
  if ( work.work_count != (int) N ) {
    std::cerr << method
              << "<" << N << "> work_count(" << work.work_count << ") failed"
              << std::endl ;
    throw std::exception();
  }
  m_flag[ work.work_rank ] = 1 ;
}

template<unsigned N>
void TEST<N>::verify()
{
  static const char method[] = "TEST::verify" ;

  for ( unsigned i = 0 ; i < N ; ++i ) {
    if ( ! m_flag[i] ) {
      std::cerr << method
                << "<" << N << "> m_flag[" << i << "] failed"
                << std::endl ;
      throw std::exception();
    }
    else {
      m_flag[i] = 0 ;
    }
  }
}

void test_tpi_cpp( int np )
{
  TEST<1> test_1 ;
  TEST<2> test_2 ;
  TEST<4> test_4 ;
  TEST<8> test_8 ;
  TEST<16> test_16 ;

  TPI::Init( np );

  TPI::Run( test_1 , & TEST<1>::flag , 1 );
  TPI::Run( test_2 , & TEST<2>::flag , 2 );
  TPI::Run( test_4 , & TEST<4>::flag , 4 );
  TPI::Run( test_8 , & TEST<8>::flag , 8 );
  TPI::Run( test_16 , & TEST<16>::flag , 16 );

  test_1.verify();
  test_2.verify();
  test_4.verify();
  test_8.verify();
  test_16.verify();

  TPI::Finalize();
}

int main( int argc , char ** argv )
{
  if ( argc ) { std::cout << argv[0] ; }
  else        { std::cout << "test" ; }
  test_tpi_cpp(1); std::cout << " 1 " ;
  test_tpi_cpp(2); std::cout << " 2 " ;
  test_tpi_cpp(4); std::cout << " 4 " ;
  test_tpi_cpp(8); std::cout << " 8 " ;
  test_tpi_cpp(16); std::cout << " 16 " ;
  std::cout << " passed" << std::endl ;
  return 0 ;
}

