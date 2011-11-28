#include <sys/time.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>


	int internalForce(	   const int, 
				           const double,
				           double,
				           double * const,
				           double * const,
				           double * const,
				           double * const, 
				           double * const, 
				           double * const, 
				           double * const,
				           double * const, 
				           double * const, 
				           double * const, 
						   double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const,
				           double * const, 
						   double * const);



/***********************************************************/

	void fill(int n, double * pos, double * vel, double * rot){

		int i;

		for(i = 0; i < n; i++){

			pos[24 * i +  0] = 0.0;
			pos[24 * i +  3] = 1.0;
			pos[24 * i +  6] = 1.0;
			pos[24 * i +  9] = 0.0;
			pos[24 * i + 12] = 0.0;
			pos[24 * i + 15] = 1.0;
			pos[24 * i + 18] = 1.0;
			pos[24 * i + 21] = 0.0;

			pos[24 * i +  1] = 0.0;
			pos[24 * i +  4] = 0.0;
			pos[24 * i +  7] = 1.0;
			pos[24 * i + 10] = 1.0;
			pos[24 * i + 13] = 0.0;
			pos[24 * i + 16] = 0.0;
			pos[24 * i + 19] = 1.0;
			pos[24 * i + 22] = 1.0;

			pos[24 * i +  2] = 0.0;
			pos[24 * i +  5] = 0.0;
			pos[24 * i +  8] = 0.0;
			pos[24 * i + 11] = 0.0;
			pos[24 * i + 14] = 1.0;
			pos[24 * i + 17] = 1.0;
			pos[24 * i + 20] = 1.0;
			pos[24 * i + 23] = 1.0;


			vel[24 * i +  0] = 0.25;
			vel[24 * i +  3] = 0.0;
			vel[24 * i +  6] = 0.0;
			vel[24 * i +  9] = 0.0;
			vel[24 * i + 12] = 0.0;
			vel[24 * i + 15] = 0.0;
			vel[24 * i + 18] = 0.0;
			vel[24 * i + 21] = 0.0;

			vel[24 * i +  1] = 0.25;
			vel[24 * i +  4] = 0.0;
			vel[24 * i +  7] = 0.0;
			vel[24 * i + 10] = 0.0;
			vel[24 * i + 13] = 0.0;
			vel[24 * i + 16] = 0.0;
			vel[24 * i + 19] = 0.0;
			vel[24 * i + 22] = 0.0;

			vel[24 * i +  2] = 0.25;
			vel[24 * i +  5] = 0.0;
			vel[24 * i +  8] = 0.0;
			vel[24 * i + 11] = 0.0;
			vel[24 * i + 14] = 0.0;
			vel[24 * i + 17] = 0.0;
			vel[24 * i + 20] = 0.0;
			vel[24 * i + 23] = 0.0;


			rot[9 * i + 0] = 0.0;
			rot[9 * i + 1] = 0.0;
			rot[9 * i + 2] = 0.0;
			rot[9 * i + 3] = 0.0;
			rot[9 * i + 4] = 0.0;
			rot[9 * i + 5] = 0.0;
			rot[9 * i + 6] = 0.0;
			rot[9 * i + 7] = 0.0;
			rot[9 * i + 8] = 0.0;

		}

	}

/***********************************************************/


	double run_kernel(int n){

		timeval start,stop,result;

		int size_32 = sizeof(double) * n * 32;
		int size_24 = sizeof(double) * n * 24;
		int size_12 = sizeof(double) * n * 12;
		int size_9 = sizeof(double)  * n *  9;
		int size_6 = sizeof(double)  * n *  6;
		int size_3 = sizeof(double)  * n *  3;
		int size_1 = sizeof(double)  * n *  1;

		int nelem = n;
	  	double * position = 		(double *)malloc(size_24);		
		double * velocity = 		(double *)malloc(size_24);
		double * mid_pos = 			(double *)malloc(size_24);
	  	double * rot_old =	 		(double *)malloc(size_9);
	  	double * rot_new = 			(double *)malloc(size_9);			
		double * gradop12 = 		(double *)malloc(size_24);
		double * vel_grad = 		(double *)malloc(size_9);
		double * stretch = 			(double *)malloc(size_6);
		double * vorticity = 		(double *)malloc(size_3);
		double * rot_stret = 		(double *)malloc(size_6);
		double * mid_vol = 			(double *)malloc(size_1);
		double * hgop = 			(double *)malloc(size_32);		
		double * strain_rate = 		(double *)malloc(size_6);
		double * shrmod = 			(double *)malloc(size_1);
		double * dilmod = 			(double *)malloc(size_1);
		double * elem_mass = 		(double *)malloc(size_1);
		double * elem_t_step = 		(double *)malloc(size_1);
		double * force_new = 		(double *)malloc(size_24);
		double * hg_energy = 		(double *)malloc(size_1);
		double * hg_resist_n = 		(double *)malloc(size_12);		//old and new
		double * hg_resist_o = 		(double *)malloc(size_12);
		double * intern_energy = 	(double *)malloc(size_1);
		double * stress_new = 		(double *)malloc(size_6);			
		double * rot_stress = 		(double *)malloc(size_6);
		double * bulk_modulus = 	(double *)malloc(size_1);
		double * two_mu = 			(double *)malloc(size_1);
		double * s_temp = 			(double *)malloc(size_6);

		double dt = 0.25;
		double	stable_time_step = 0;

		fill(n, position, velocity, rot_old);

		gettimeofday(&start, NULL);

		internalForce(	nelem, 
						dt, 
						stable_time_step, 
						elem_t_step, 
						position, 
						velocity, 
						rot_old, 
						rot_new, 
						mid_vol, 
						vorticity, 
						stretch, 
						strain_rate,
						hgop, 
						stress_new,
						rot_stret,
						bulk_modulus,
						two_mu,
						shrmod,
						dilmod,
						elem_mass,
						force_new,
						hg_energy,
						intern_energy,
						hg_resist_o,
						hg_resist_n);

		gettimeofday(&stop, NULL);
		timersub(&stop, &start, &result);
		double time = (result.tv_sec + result.tv_usec/1000000.0);


/***********************************************************/


		free(position);
		free(velocity);
		free(mid_pos);
		free(rot_old);
		free(rot_new);
		free(gradop12);
		free(vel_grad);
		free(stretch);
		free(vorticity);
		free(rot_stret);
		free(mid_vol);
		free(hgop);
		free(strain_rate);
		free(shrmod);
		free(dilmod);
		free(elem_mass);
		free(elem_t_step);
		free(force_new);
		free(hg_energy);
		free(hg_resist_o);
		free(hg_resist_n);
		free(intern_energy);
		free(stress_new);
		free(rot_stress);
		free(bulk_modulus);
		free(two_mu);
		free(s_temp);

		return time;

	}


/************************************************************/

namespace test{

	void test_Original_Host( int beg, int end, int runs){

		std::cout << "Original Host: " << std::endl;

		for(int i = beg; i < end; i++){

			int n = 1 << i;
			double time = 100000, min = 100000;

			for(int j = 0; j < runs; j++){

				time = run_kernel(n);
				if (time < min)
					min = time;

			}

			std::cout << 	std::setw(8) << n << ", " << 
							std::setw(8) << 1000 * time << ", " << 
							std::setw(8) << time / n << std::endl;

		}//for

	}//test_host

}// namespace

