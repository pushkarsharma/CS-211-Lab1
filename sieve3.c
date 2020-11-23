/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[])
{
   unsigned long int count; /* Local prime count */
   double elapsed_time;     /* Parallel execution time */
   unsigned long int first; /* Index of first multiple */
   int local_first;
   unsigned long int global_count = 0; /* Global prime count */
   unsigned long long int high_value;  /* Highest value on this proc */
   unsigned long int i;
   int id;                           /* Process ID number */
   unsigned long int index;          /* Index of current prime */
   unsigned long long int low_value; /* Lowest value on this proc */
   char *marked;                     /* Portion of 2,...,'n' */
   char *local_prime_marked;
   unsigned long long int n;     /* Sieving from 2, ..., 'n' */
   int p;                        /* Number of processes */
   unsigned long int proc0_size; /* Size of proc 0's subarray */
   unsigned long int prime;
   unsigned long int local_prime; /* Current prime */
   unsigned long int size;        /* Elements in 'marked' */
   unsigned long int local_prime_size;
   unsigned long int local_size;
   char *local_marked;
   unsigned long int blocked_size = 1048576;
   unsigned long long int blocked_low_value;
   unsigned long long int blocked_high_value;

   MPI_Init(&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2)
   {
      if (!id)
         printf("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit(1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   /* Add you code here  */
   low_value = 2 + ((id) * (n - 1) / (p));
   high_value = 2 + (((id + 1) * (n - 1) / (p)) - 1);
   low_value = low_value + (low_value + 1) % 2;
   high_value = high_value - (high_value + 1) % 2;
   size = (high_value - low_value) / 2 + 1;
   local_size = (int)sqrt(n) - 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */
   proc0_size = (n / 2 - 1) / p;
   if ((2 + proc0_size) < (int)sqrt((double)n / 2))
   {
      if (!id)
         printf("Too many processes.\n");
      MPI_Finalize();
      exit(1);
   }

   /* Allocate this process's share of the array. */
   marked = (char *)malloc(size);
   local_marked = (char *)malloc(local_size);
   if (marked == NULL || local_marked == NULL)
   {
      printf("Cannot allocate enough memory.\n");
      MPI_Finalize();
      exit(1);
   }

   /*Calculating the sieving primes within each machine separately.*/
   local_prime = 2;
   for (i = 0; i < local_size; i++)
      local_marked[i] = 0;
   index = 0;
   do
   {
      local_first = local_prime * local_prime - 2;
      for (i = local_first; i < local_size; i += local_prime)
         local_marked[i] = 1;
      while (local_marked[++index] == 1)
         ;
      local_prime = 2 + index;
   } while (local_prime * local_prime <= n);

   for (i = 0; i < size; i++)
      marked[i] = 0;

   blocked_low_value = low_value;
   blocked_high_value = blocked_low_value + 2 * (blocked_size - 1);
   do
   {
      index = 0;
      prime = 3;
      while (prime * prime <= blocked_high_value)
      {
         if (prime * prime > blocked_low_value)
            first = (prime * prime - blocked_low_value) / 2;
         else
         {
            /*Find the offset in the array as in how far we are from
              the next multiple of the current seiving prime number.*/
            if (!(blocked_low_value % prime))
               first = 0;
            else
               first = (prime - (blocked_low_value % prime) + blocked_low_value / prime % 2 * prime) / 2;
         }
         /*Mark the multiples of the seiving prime.*/
         for (i = first + (blocked_low_value - low_value) / 2; i <= (blocked_high_value - low_value) / 2; i += prime)
            marked[i] = 1;
         while (local_marked[++index] == 1)
            ;
         /*The mapping of the index and value is index+2. So, we
           calculate the next prime value with the next un-marked
           index.*/
         prime = index + 2;
      }
      blocked_low_value = blocked_high_value + 2;
      blocked_high_value = blocked_low_value + 2 * (blocked_size - 1);
      if (blocked_high_value > high_value)
         blocked_high_value = high_value;
   } while (blocked_low_value <= high_value);
   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i])
         count++;
   if (!id)
      count++;
   if (p > 1)
      MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

   /* Stop the timer */
   elapsed_time += MPI_Wtime();

   /* Print the results */

   if (!id)
   {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);
   }
   MPI_Finalize();
   return 0;
}
