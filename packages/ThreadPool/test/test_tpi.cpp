/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/* Redistribution and use in source and binary forms, with or without     */
/* modification, are permitted provided that the following conditions are */
/* met:                                                                   */
/*                                                                        */
/* 1. Redistributions of source code must retain the above copyright      */
/* notice, this list of conditions and the following disclaimer.          */
/*                                                                        */
/* 2. Redistributions in binary form must reproduce the above copyright   */
/* notice, this list of conditions and the following disclaimer in the    */
/* documentation and/or other materials provided with the distribution.   */
/*                                                                        */
/* 3. Neither the name of the Corporation nor the names of the            */
/* contributors may be used to endorse or promote products derived from   */
/* this software without specific prior written permission.               */
/*                                                                        */
/* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY        */
/* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR     */
/* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE    */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    */
/* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR     */
/* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
/* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING   */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           */
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
  if ( work.count != (int) N ) {
    std::cerr << method
              << "<" << N << "> count(" << work.count << ") failed"
              << std::endl ;
    throw std::exception();
  }
  m_flag[ work.rank ] = 1 ;
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

