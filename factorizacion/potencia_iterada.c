#include <stdio.h>
#include <math.h>

#define NDEBUG 1
#include "dbg.h"
#define ARRAY_SIZE(array) (sizeof (array) / sizeof (array[0]))
#define EPSILON	10e-9
#define MAX_IT	100

struct crs_matrix
{
	double *aa;
	size_t saa;
	
	/* jad must have aa length */
	size_t *jad;
	
	size_t *jr;
	size_t sjr;

	size_t *dia;

	size_t n;
};

struct crs_system
{
	struct crs_matrix *A;
	double *x;
	double *b;
};

void pi(struct crs_matrix *A, double *q0, double *x)
{
	size_t i, j, k, jr_i, jr_i_1;
	double a_ii, x_i, l_ij, l_jk, l_jj;
	double norm_max;
	size_t norm_i;
	double beta, q;
	double d_max, q_max, er;
	for(k = 0; k < MAX_IT; k++)
	{
		for(i = 0; i < A->n; i++)
		{
			debug("\ni = %d\n", i);
			jr_i = A->jr[i];
			jr_i_1 = A->jr[i+1];
			debug("jr[%d] = %d | jr[%d] = %d\n", i, jr_i, i+1, jr_i_1);

			x[i] = 0.0;
			for(j = jr_i; j < jr_i_1; j++)
			{
				debug("aa[%d] = %.2f\n", j, A->aa[j]);
				debug("aa[%d] is at A[%d,%d]\n", j, i, A->jad[j]);
				debug("%d =?= %d\n", i, A->jad[j]);

				x[i] += A->aa[j] * q0[A->jad[j]];
				
			}
			printf("%.2f ", x[i]);
		}
		printf("| ");
		
		norm_max = 0.0;
		for(i = 0; i < A->n; i++)
		{
			if(norm_max < fabs(x[i]))
			{
				norm_max = fabs(x[i]);
				norm_i = i;
			}
		}
		beta = x[norm_i];
		d_max = 0.0;
		q_max = 0.0;
		for(i = 0; i < A->n; i++)
		{
			q = x[i] / beta;

			if(d_max < fabs(q0[i] - q))
				d_max = fabs(q0[i] - q);
			
			if(q_max < fabs(q))
				q_max = q;
			
			printf("%.2f ", q);
			q0[i] = q;
		}
		er = d_max / q_max;
		printf("| err = %e\n", er);

		if(er < EPSILON)
		{
			printf("%d iteraciones. Error %e\n", k, er);
			break;	
		}
	}
}

int main(int argc, char *argv[])
{
/*
	double  aa[] = { 10,  0,  4,
	                 -9,  2,  2,
			  1, -1,  2 };

	size_t jad[] = { 0, 1, 2,
	                 0, 1, 2,
	                 0, 1, 2 };

	size_t dia[] = { 0, 4, 8 };
	size_t  jr[] = { 0, 3, 6, 9 };
	double  q0[] = { 1, 1, 1 };
	double   x[] = { 0, 0, 0 };
*/
	double  aa[] = { -4, 14,  0,
	                 -5, 13,  0,
			 -1,  0,  2 };

	size_t jad[] = { 0, 1, 2,
	                 0, 1, 2,
	                 0, 1, 2 };

	size_t dia[] = { 0, 4, 8 };
	size_t  jr[] = { 0, 3, 6, 9 };
	double  q0[] = { 1, 1, 1 };
	double   x[] = { 0, 0, 0 };

	struct crs_matrix A;
	A.aa = aa;
	A.saa = ARRAY_SIZE(aa);
	A.jad = jad;
	A.jr = jr;
	A.sjr = ARRAY_SIZE(jr);
	A.dia = dia;
	A.n = 3;
	
	pi(&A, q0, x);
	return 0;
}
