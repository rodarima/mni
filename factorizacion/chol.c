#include <stdio.h>
#include <stdlib.h>

//Funciones inline para get y set
#define HAVE_INLINE

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#define NDEBUG 1
#include "dbg.h"

void chol(gsl_matrix *A, gsl_matrix *L, size_t n)
{
	size_t i, j, k;
	double s, l_ii, l_ik, l_ij, l_jk, l_jj;
	for(i = 0; i < n; i++)
	{
		debug("\ni = %d\n", i);
		/* Calcular l_ij */

		for(j = 0; j < i; j++)
		{
			l_ij = gsl_matrix_get(A, i, j);
			debug("1 Ai%dj%d = %f\n", i, j, l_ij);
			for(k = 0; k < j; k++)
			{
				debug("1.1 - Li%dk%d * Lj%dk%d\n", i, k, j, k);
				l_ik = gsl_matrix_get(L, i, k);
				l_jk = gsl_matrix_get(L, j, k);
				l_ij -= l_ik * l_jk;
			}
			l_jj = gsl_matrix_get(L, j, j);
			debug("1 / Lj%dj%d = %f\n", j, j, l_jj);
			l_ij /= l_jj;

			/* XXX Set */
			gsl_matrix_set(L, i, j, l_ij);
			debug("1 Li%dj%d = %f\n", i, j, l_ij);
		}

		/* Calcular l_ii */

		l_ii = gsl_matrix_get(A, i, i);
		debug("2 Ai%di%d = %f\n", i, i, l_ii);
		for(k = 0; k < i; k++)
		{
			l_ik = gsl_matrix_get(L, i, k);
			debug("2.1 - Li%dk%d^2 = %f^2\n", i, k, l_ik);
			l_ii -= l_ik * l_ik;
		}
		l_ii = sqrt(l_ii);
		debug("2 Li%di%d = %f\n", i, i, l_ii);

		/* XXX Set */
		gsl_matrix_set(L, i, i, l_ii);
	}

	if(A != L)
	{
		for(i = 0; i < n; i++)
		{
			for(j = i+1; j < n; j++)
			{
				gsl_matrix_set(L, i, j, 0.0);
			}
		}
	}
}

int main(int argc, char *argv[])
{
	double a_data[] =
	{
		 4,10, 4, 2,
		10,61,28,53,
		 4,28,38,41,
		 2,53,41,78
	};
	gsl_matrix *L = gsl_matrix_calloc(4, 4);
	gsl_matrix *U = gsl_matrix_calloc(4, 4);
	gsl_matrix_view m = gsl_matrix_view_array(a_data, 4, 4);
	
	chol(&m.matrix, L, 4);

	int i,j;
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("% 10f ", gsl_matrix_get(&m.matrix, i, j));
		}
		printf("\n");
	}
	printf("\n");
	
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("% 10f ", gsl_matrix_get(L, i, j));
		}
		printf("\n");
	}
	return 0;
}
