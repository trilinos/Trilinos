
namespace Test {

template<class Scalar , class DeviceType >
struct Determinant;

template<class Scalar >
struct Determinant<Scalar , KOKKOS_MACRO_DEVICE >
{
	typedef KOKKOS_MACRO_DEVICE 		device_type;
	typedef device_type::size_type 		size_type;
	
	typedef  Kokkos::MDArrayView<Scalar,device_type> device_array;
	typedef  Kokkos::MDArrayView<Scalar,Kokkos::DeviceHost> host_array;
	
  private:
	
	device_array out;
	device_array in;
	int dim_i1;

  public:
  
  	Determinant(device_array & arg_out ,const device_array & arg_in) : out(arg_out) , in(arg_in) 
  	{
  		//Assume rank 4
		dim_i1 = in.dimension(1);
  	}
  	//Assuming 3 dimension
  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(size_type ielem) const {
        for (int i1=0; i1<dim_i1; i1++) {
          int i,j,rowID = 0;
          int colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0), determinant(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( in(ielem,i1,i,j) ) >  emax){
                rowID = i;  colID = j; emax = std::abs( in(ielem,i1,i,j) );
              }
            }
          }
          if( emax > 0 ){
            if( rowID ){
              rowperm[0] = rowID;
              rowperm[rowID] = 0;
            }
            if( colID ){
              colperm[0] = colID;
              colperm[colID] = 0;
            }
            Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
            for(i=0; i < 3; ++i){
              for(j=0; j < 3; ++j){
                B[i][j] = in(ielem,i1,rowperm[i],colperm[j]);
              }
            }
            B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
            for(i=0; i < 2; ++i){
              for(j=0; j < 2; ++j){
                S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
              }
            }
            determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
            if( rowID ) determinant = -determinant;
            if( colID ) determinant = -determinant;
          }
          out(ielem,i1) = determinant;
        } // for i1
  	}
};

} //Test
