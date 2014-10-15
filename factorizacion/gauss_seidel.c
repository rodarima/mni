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

void gs(struct crs_system *s, double *x2)
{
	struct crs_matrix *A = s->A;
	double *x = s->x;
	double *b = s->b;
	double d_max = 0.0;
	double x_1_max = 0.0;
	double er;

	size_t i, j, k, jr_i, jr_i_1;
	double a_ii, x_i, l_ij, l_jk, l_jj;
	for(k = 0; k < MAX_IT; k++)
	{
		d_max = 0.0;
		x_1_max = 0.0;

		for(i = 0; i < A->n; i++)
		{
			debug("\ni = %d\n", i);
			jr_i = A->jr[i];
			jr_i_1 = A->jr[i+1];
			debug("jr[%d] = %d | jr[%d] = %d\n", i, jr_i, i+1, jr_i_1);

			x2[i] = x[i];

			x[i] = b[i];

			for(j = jr_i; j < A->dia[i]; j++)
			{
				debug("1|aa[%d] = %.2f\n", j, A->aa[j]);
				debug("1|aa[%d] is at A[%d,%d]\n", j, i, A->jad[j]);
				debug("1|%d =?= %d\n", i, A->jad[j]);

				x[i] -= A->aa[j] * x[A->jad[j]];
			}
			for(j = A->dia[i] + 1; j < jr_i_1; j++)
			{
				debug("2|aa[%d] = %.2f\n", j, A->aa[j]);
				debug("2|aa[%d] is at A[%d,%d]\n", j, i, A->jad[j]);
				debug("2|%d =?= %d\n", i, A->jad[j]);

				x[i] -= A->aa[j] * x[A->jad[j]];
			}
			x[i] /= A->aa[A->dia[i]];
		}
		for(i = 0; i < A->n; i++)
		{
			if(d_max < fabs(x[i] - x2[i]))
				d_max = fabs(x[i] - x2[i]);
			if(x_1_max < fabs(x[i]))
				x_1_max = fabs(x[i]);
			printf("%.3f ", x[i]);
		}
		er = d_max / x_1_max;
		printf("| err = %.9f\n", er);

		if(er < EPSILON)
		{
			printf("%d iteraciones. Error %.9f\n", k, er);
			break;	
		}
	}
}

int main(int argc, char *argv[])
{
/*
	double  aa[] = { 1, 2, 3, 4, -5, 6.2, 7, 8, -9.1, 10, 11, 12.3 };
	size_t dia[] = { 0, 3, 6, 10, 11 };
	size_t jad[] = { 0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4 };
	size_t  jr[] = { 0, 2, 5, 9, 11, 12 };
	double   b[] = { 1, 2, 3, 4, 5 };
	double   x[] = { 0, 0, 0, 0, 0 };
	double  x2[] = { 0, 0, 0, 0, 0 };
*/
	double  aa[] = { 10, -2, -1, -1,
	                 -2, 10, -1, -1,
			 -1, -1, 10, -2,
			 -1, -1, -2, 10 };

	size_t jad[] = { 0, 1, 2, 3,
	                 0, 1, 2, 3,
	                 0, 1, 2, 3,
	                 0, 1, 2, 3 };

	size_t dia[] = { 0, 5, 10, 15 };
	size_t  jr[] = { 0, 4, 8, 12, 16 };
	double   b[] = { 3,15,27,-9 };
	double   x[] = { 0, 0, 0, 0 };
	double  x2[] = { 0, 0, 0, 0 };

	struct crs_matrix A;
	A.aa = aa;
	A.saa = ARRAY_SIZE(aa);
	A.jad = jad;
	A.jr = jr;
	A.sjr = ARRAY_SIZE(jr);
	A.dia = dia;
	A.n = 4;

	struct crs_system s;
	s.A = &A;
	s.b = b;
	s.x = x;
	
	gs(&s, x2);
	return 0;
}
