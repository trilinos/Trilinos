/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <iostream>
#include <fstream>

// Print out sparse Matrix (A) to a file. Useful for visualizing matrix
// to see correctness.
template<class scalar_vector , class int_vector>
void printSparse(const std::string & filename , scalar_vector & value , 
												int_vector & row ,
												int_vector & col )
{
	std::ofstream outfile(filename.c_str());
	int end = row.length()-1;
	outfile<<end << " " << col.length()<<std::endl;
	for( int i = 0 ; i < end ; i++)
	{
		int stop = row(i+1);
		for(int j = row(i) ; j < stop ; j++) 
		{
			if(value(j) != 0)
				outfile << i+1 <<" " << col(j)+1 << " " << value(j)<< std::endl;
		}
	}
}

// Print out answer (X) to a file in a specialized format for viewing a GLUT visualization
// of answer.
template<class Scalar , class scalar_vector_d, class HostView_scalar, class HostView_int>
void printGLUT(const std::string & filename , scalar_vector_d & X , HostView_scalar & elem_coords_h, 
											HostView_int & elem_nodeIDs_h , int x , int y, int z)
{
	typedef Kokkos::MultiVectorView<Scalar , Kokkos::DeviceHost> scalar_vector_h;

	int nelem = x * y * z;
	int nnodes = X.length();
	std::ofstream outfile(filename.c_str());
	outfile<<x<<" "<<y<<" "<<z<< " " << 1<<std::endl;
	scalar_vector_h X_host = Kokkos::create_labeled_multivector<scalar_vector_h>("X_host", nnodes);
	Kokkos::deep_copy(X_host , X);
	for(int i = 0 ; i < nelem ; i++) 
	{
			for(int j = 0 ; j < 8 ; j++)
			{
				//Coordinates temperature
				outfile << elem_coords_h(i,0,j) << " " 
				<< elem_coords_h(i,1,j)<< " "
				<< elem_coords_h(i,2,j) <<" " 
				<< elem_nodeIDs_h(i,j) << std::endl;
			}
	}

	for(int i = 0 ; i < nnodes ; i++) 
	{
		outfile << i << " " << X_host(i) << std::endl;
	}
	outfile.close();
} 

