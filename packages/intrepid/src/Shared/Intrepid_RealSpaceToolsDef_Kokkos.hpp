#ifndef INTREPID_REALSPACETOOLSDEF_KOKKOS_HPP
#define INTREPID_REALSPACETOOLSDEF_KOKKOS_HPP

#ifdef INTREPID_OLD_KOKKOS_CODE
#ifdef Kokkos_ENABLE_CXX11
#include "Intrepid_RealSpaceTools.hpp"

namespace Intrepid{
template<class Scalar>	
template<class Scalar1, class Scalar2,class Layout,class MemorySpace>
void RealSpaceTools<Scalar>::inverseTemp(Kokkos::View<Scalar1,Layout,MemorySpace> & inverseMats, const Kokkos::View<Scalar2,Layout,MemorySpace> & inMats){
		
 ArrayWrapper<Scalar1,Kokkos::View<Scalar1,Layout,MemorySpace>, Rank<Kokkos::View<Scalar1,Layout,MemorySpace> >::value, false>inverseMatsWrap(inverseMats);
 ArrayWrapper<Scalar2,Kokkos::View<Scalar2,Layout,MemorySpace>, Rank<Kokkos::View<Scalar2,Layout,MemorySpace> >::value, true>inMatsWrap(inMats);

  int arrayRank = getrank(inMats);

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.dimension(arrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = inMats.dimension(0);
      dim_i1 = inMats.dimension(1);
       switch(dim) {
    case 3: {
     
	Kokkos::parallel_for (dim_i0, [&] (const int i0) {
  
        for (int i1=0; i1<dim_i1; i1++) {
         

          int i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i0,i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i0,i1,i,j) );
              }
            }
          }

          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(i0,i1,rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];


          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i0,i1,i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        } // for i1
      }); // for i0
      break;
    } // case 3

    case 2: {
     
	Kokkos::parallel_for (dim_i0, [&] (const int i0) {


        for (int i1=0; i1<dim_i1; i1++) {
 

          Scalar determinant    = inMatsWrap(i0,i1,0,0)*inMatsWrap(i0,i1,1,1)-inMatsWrap(i0,i1,0,1)*inMatsWrap(i0,i1,1,0);

          inverseMatsWrap(i0,i1,0,0)   = inMatsWrap(i0,i1,1,1) / determinant;
          inverseMatsWrap(i0,i1,0,1) = - inMatsWrap(i0,i1,0,1) / determinant;
          //
          inverseMatsWrap(i0,i1,1,0) = - inMatsWrap(i0,i1,1,0) / determinant;
          inverseMatsWrap(i0,i1,1,1) =   inMatsWrap(i0,i1,0,0) / determinant;
        } // for i1
      }); // for i0
      break;
    } // case 2

    case 1: {
       	Kokkos::parallel_for (dim_i0, [&] (const int i0) {

        for (int i1=0; i1<dim_i1; i1++) {
          inverseMatsWrap(i0,i1,0,0) = (Scalar)1 / inMatsWrap(i0,i1,0,0);
        } // for i1
      }); // for i2
          
      
      break;
    } // case 1

  } // switch (dim)	
      break;
    case 3:
      dim_i1 = inMats.dimension(0);
       switch(dim) {
    case 3: {
	Kokkos::parallel_for (dim_i1, [&] (const int i1) {

          int i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i1,i,j) );
              }
            }
          }

          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(i1,rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];


          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i1,i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        }); // for i1
      
      break;
    } // case 3

    case 2: {

       Kokkos::parallel_for (dim_i1, [&] (const int i1) {

          Scalar determinant    = inMatsWrap(i1,0,0)*inMatsWrap(i1,1,1)-inMatsWrap(i1,0,1)*inMatsWrap(i1,1,0);

          inverseMatsWrap(i1,0,0)   = inMatsWrap(i1,1,1) / determinant;
          inverseMatsWrap(i1,0,1) = - inMatsWrap(i1,0,1) / determinant;
          //
          inverseMatsWrap(i1,1,0) = - inMatsWrap(i1,1,0) / determinant;
          inverseMatsWrap(i1,1,1) =   inMatsWrap(i1,0,0) / determinant;
        }); // for i1
 
      break;
    } // case 2

    case 1: {
   
      	Kokkos::parallel_for (dim_i1, [&] (const int i1) {
        inverseMatsWrap(i1,0,0) = (Scalar)1 / inMatsWrap(i1,0,0); 
        });
    
          
      
       break;
      } // case 1

    } // switch (dim)	
      break;
    case 2:
     switch(dim) {
    case 3: {
          int i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMatsWrap(i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMatsWrap(i,j) );
              }
            }
          }

          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMatsWrap(rowperm[i],colperm[j]);
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];

          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMatsWrap(i,j) = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
   
      break;
    } // case 3

    case 2: {
          Scalar determinant    = inMatsWrap(0,0)*inMatsWrap(1,1)-inMatsWrap(0,1)*inMatsWrap(1,0);

          inverseMatsWrap(0,0)   = inMatsWrap(1,1) / determinant;
          inverseMatsWrap(0,1) = - inMatsWrap(0,1) / determinant;
          //
          inverseMatsWrap(1,0) = - inMatsWrap(1,0) / determinant;
          inverseMatsWrap(1,1) =   inMatsWrap(0,0) / determinant;
      
      break;
    } // case 2

    case 1: {
           inverseMatsWrap(0,0) = (Scalar)1 / inMatsWrap(0,0);       
      break;
    } // case 1

  }
    break;  
  }

	

	
}
	
}


#endif
#endif
#endif
