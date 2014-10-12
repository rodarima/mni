#include <stdio.h>
#include <stdlib.h>

//Funciones inline para get y set
#define HAVE_INLINE

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#define NDEBUG 1
#include "dbg.h"

#define ARRAY_SIZE(array) (sizeof (array) / sizeof (array[0]))

double skyline_matrix_get(gsl_vector *aa, gsl_vector *jd, size_t n, size_t i, size_t j)
{
	double jd_i = gsl_vector_get(jd, i);
	debug("jd[%d] = %f\n", i, jd_i);

}

void chol_sky(gsl_vector *aa, gsl_vector *jd, size_t n)
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
		 3, 0,-1, 0, 0,
		 0, 1, 0, 0, 1,
		-1, 0, 3, 0, 0,
		 0, 0, 0, 5, 0,
		 0, 1, 0, 0, 8
	};

	double aa_data[] = { 3, 1, -1, 0, 3, 5, 1, 0, 0, 8 };
	double jd_data[] = { 0, 1, 4, 5, 9 };
	int N = 5;

	gsl_vector_view aa_v, jd_v;

	aa_v = gsl_vector_view_array(aa_data, ARRAY_SIZE(aa_data));
	jd_v = gsl_vector_view_array(jd_data, ARRAY_SIZE(jd_data));

	chol_sky(&aa_v.vector, &jd_v.vector, N);

	return 0;
}
