
/*TEST
PATH='tests/trapezoidal.c'
CCFLAGS=""
INPUT=""
OUTPUT="With n = 1024 trapezoids, our estimate\nof the intergral from 0.000000 to 1.000000 = -1.799356"
STATUS=0
SKIP=1
TEST*/

/*===========================================================================================*/
/* Parallel Trapezoidal Rule
 * Peter S. Pacheco
 *
 * Input: None.
 * Output: Estimate of the intergral from a to b of f(x)
 *          using the trapezoidal rule and n trapezoids.
 *
 * Algorithm:
 *   1. Each process calculates "its" interval of integration.
 *   2. Each process estimates the integral of f(x) over its
 *      interval using the trapezoidal rule.
 *   3a. Each process != 0 sends its integral to 0.
 *   3b. Process 0 sums the calculations received from 
 *       the individual processes and prints the result.
 * Note: f(x), a, b, and n are all hardwired.
 */
  #include <stdio.h>
  #include "mpi.h"
/*===========================================================================================*/
  float Trap(float local_a, float local_b, int local_n, float h)
  {
    float integral; 
    float x; 
    int i;

    float f(float x); /*function we are integrating */
    integral = (f(local_a)+f(local_b))/2.0;
    x = local_a;
    for (i = 1; i<=local_n-1; i++) {
      x = x+h;
      integral = integral+ f(x);
    }
    integral = integral*h;
    return integral;
  } /*Trap*/
/*===========================================================================================*/
  float f(float x) {
    float return_val;
    /*Calculate f(x) */
    return_val = x*x*x + (50*x);
    return return_val;
  }
/*===========================================================================================*/
  int main(int argc, char**argv)
  {
    int my_rank;
    int p;
    float a = 0.0;
    float b = 1.0;
    int n = 1024;
    float h;
    float local_a;
    float local_b;
    int local_n;
    float integral;
    float total;
    int source;
    int dest = 0;
    int tag = 0;
    MPI_Status status;

    float Trap(float local_a, float local_b, int local_n, float h); /*Calculate local integral*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    h = (b-a)/n;    /*h is the same for all processes*/
    local_n = n/p;  /*So is the number of trapezoids */
    
    /*Length of each process's interval of integration = local_n*h. So my interval starts at:*/
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral - Trap(local_a, local_b, local_n, h);

    /*Add up the integrals calculated by each process*/
    if (my_rank == 0) {
      total = integral; 
      for (source = 1; source < p; source++) {
        MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
        total = total + integral;
      }
    } else {
      MPI_Send(&integral, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    }

    /*Print the result*/
    if (my_rank == 0) {
      printf("With n = %d trapezoids, our estimate\n",n);
      printf("of the intergral from %f to %f = %f\n",a, b, total);
    }
    
    /*Shut down MPI*/
    MPI_Finalize();
    return 0;
  } /*main*/

    
/*===========================================================================================*/
