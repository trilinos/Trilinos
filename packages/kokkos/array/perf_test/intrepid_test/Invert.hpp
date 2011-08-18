
namespace Test {

template< class Scalar , class DeviceType, int dimension>
struct Invert;

template<class Scalar , int dimension>
struct Invert<Scalar, KOKKOS_MACRO_DEVICE, dimension>{

	typedef KOKKOS_MACRO_DEVICE     device_type ;
  	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

	array_type inMatrices;
	array_type outMatrices;
	
	int cubature;

  	Invert(array_type & arg_inMat, array_type & arg_outMat, int cub){

		inMatrices = arg_inMat;    	
		outMatrices = arg_outMat;    	

		cubature = cub;

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()( int ielem )const {

		switch(dimension){

			case 3:
				
				for(int c = 0; c < cubature; c++){

	 				int i, j, rowID = 0, colID = 0;
					int rowperm[3]={0,1,2};
					int colperm[3]={0,1,2}; // Complete pivoting
					Scalar emax(0);

					for(i = 0; i < 3; i++){
						for(j = 0; j < 3; j++){

							Scalar entry_abs = sqrt(pow(inMatrices(ielem, c, i, j), 2));

							if(entry_abs > emax){

								rowID = i;
								colID = j;
								emax = entry_abs;

							} // if

						} // for, j
					} // for, i

					if( rowID ){
						rowperm[0] = rowID;
						rowperm[rowID] = 0;
					}
				  	
					if( colID ){
						colperm[0] = colID;
						colperm[colID] = 0;
					}

					Scalar B[3][3], S[2][2], Bi[3][3]; // B = rowperm inMat colperm, S = Schur complement(Boo)

					for(i=0; i < 3; ++i){
						for(j=0; j < 3; ++j){

							B[i][j] = inMatrices(ielem,c,rowperm[i],colperm[j]);

						}
					}
				
				  	B[1][0] /= B[0][0]; 
					B[2][0] /= B[0][0];// B(:,0)/=pivot
					for(i=0; i < 2; ++i){
						for(j=0; j < 2; ++j){

							S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'

						} // for, j
					} // for, i

					B[1][0] /= B[0][0]; 
					B[2][0] /= B[0][0];// B(:,0)/=pivot
					for(i=0; i < 2; ++i){
						for(j=0; j < 2; ++j){

							S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'

						} // for, j
					} // for, i

			  		Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];

				  	Si[0][0] =  S[1][1]/detS;	Si[0][1] = -S[0][1]/detS;
				  	Si[1][0] = -S[1][0]/detS;	Si[1][1] =  S[0][0]/detS;

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

					  		outMatrices(ielem, c, i, j) = Bi[colperm[i]][rowperm[j]]; // set inverse

						}
				  	}

				} // for, c 	

				break;
		} // switch, dimension

	}

}; // struct

} // namespace

