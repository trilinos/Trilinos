
#include <utility>
#include <iostream>

std::pair<int,int> map_a_to_ab( int n , int i , int j )
{
  std::pair<int,int> result ;

  const int nb = ( n + 1 ) / 2 ;
  const int diag = j - i ;

  if ( diag <= - nb || ( 0 <= diag && diag < nb ) ) {
    result.first = i ;
    result.second = ( n + diag ) % n ;
  }
  else {
    result.first = j ;
    result.second = ( n - diag ) % n ;
  }

  return result ;
}

std::pair<int,int> map_ab_to_a( int n , int i , int k )
{
  std::pair<int,int> result ;

  const int nb = ( n + 1 ) / 2 ;

  result.first = i ;
  result.second = ( i + k ) % n ;

  return result ;
}

int main()
{
  const int n = 5 ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < n ; ++i ) {
      std::pair<int,int> y = map_a_to_ab( n , j , i );

      std::cout << " (" << y.first << "," << y.second << ")" ;
    }
    std::cout << std::endl ;
  }

  std::cout << std::endl ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < n ; ++i ) {
      std::pair<int,int> x = map_a_to_ab( n , i , j );

      std::cout << " (" << x.first << "," << x.second << ")" ;
    }
    std::cout << std::endl ;
  }

  std::cout << std::endl ;

  const int nb = ( n + 1 ) / 2 ;

  for ( int j = 0 ; j < n ; ++j ) {
    for ( int i = 0 ; i < nb ; ++i ) {
      std::pair<int,int> x = map_ab_to_a( n , j , i );
      std::cout << " (" << x.first << "," << x.second << ")" ;
    }
    std::cout << std::endl ;
  }
}

