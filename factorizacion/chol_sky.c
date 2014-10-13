#include <stdio.h>
#include <stdlib.h>

//Funciones inline para get y set
#define HAVE_INLINE

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#define NDEBUG 1
#include "dbg.h"

#define ARRAY_SIZE(array) (sizeof (array) / sizeof (array[0]))
#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })


double skyline_get(gsl_vector *aa, gsl_vector *jd, size_t i, size_t j)
{
	double jd_i;
	size_t x;
	double aa_x;

	jd_i = gsl_vector_get(jd, i);
	//debug("jd[%d] = %f\n", i, jd_i);
	x = jd_i - i + j;
	aa_x = gsl_vector_get(aa, x);
	debug("get(%d,%d) = aa[%d] = %f\n", i, j, x, aa_x);

	return aa_x;
}

double skyline_set(gsl_vector *aa, gsl_vector *jd,
	size_t i, size_t j, double value)
{
	double jd_i;
	size_t x;

	jd_i = gsl_vector_get(jd, i);
	//debug("jd[%d] = %f\n", i, jd_i);
	x = jd_i - i + j;
	debug("set(%d,%d) = aa[%d] = %f\n", i, j, x, value);
	gsl_vector_set(aa, x, value);
}

gsl_vector *skyline_st(gsl_vector *jd)
{
	gsl_vector *st;
	size_t i;

	if(!(st = gsl_vector_alloc(jd->size)))
	{
		perror("gsl_vector_alloc");
		return NULL;
	}
	gsl_vector_set(st, 0, 0);
	debug("st[%d] = %f\n", 0, 0);
	for(i = 1; i < jd->size; i++)
	{
		double jd_i = gsl_vector_get(jd, i);
		double jd_i_1 = gsl_vector_get(jd, i - 1);
		double st_i = i - (jd_i - jd_i_1) + 1;
		debug("st[%d] = %f\n", i, st_i);
		gsl_vector_set(st, i, st_i);
	}

	return st;
}

void chol_sky(gsl_vector *aa, gsl_vector *jd, gsl_vector *st, size_t n)
{
	size_t i, j, k, j_init, j2_init;
	double s, l_ii, l_ik, l_ij, l_jk, l_jj;
	for(i = 0; i < n; i++)
	{
		debug("\ni = %d\n", i);
		/* Calcular l_ij */
		j_init = gsl_vector_get(st, i);
		for(j = j_init; j < i; j++)
		{
			debug("\n## %d,%d\n", i, j);

			l_ij = skyline_get(aa, jd, i, j);
			debug("1 Ai%dj%d = %f\n", i, j, l_ij);
			j2_init = gsl_vector_get(st, j);
			for(k = max(j_init, j2_init); k < j; k++)
			{
				debug("1.1 - Li%dk%d * Lj%dk%d\n", i, k, j, k);
				l_ik = skyline_get(aa, jd, i, k);
				l_jk = skyline_get(aa, jd, j, k);
				l_ij -= l_ik * l_jk;
			}
			l_jj = skyline_get(aa, jd, j, j);
			debug("1 / Lj%dj%d = %f\n", j, j, l_jj);
			l_ij /= l_jj;

			/* XXX Set */
			skyline_set(aa, jd, i, j, l_ij);
			debug("1 Li%dj%d = %f\n", i, j, l_ij);
		}

		/* Calcular l_ii */

		debug("\n## %d,%d\n", i, i);
		l_ii = skyline_get(aa, jd, i, i);
		debug("2 Ai%di%d = %f\n", i, i, l_ii);
		for(k = j_init; k < i; k++)
		{
			l_ik = skyline_get(aa, jd, i, k);
			debug("2.1 - Li%dk%d^2 = %f^2\n", i, k, l_ik);
			l_ii -= l_ik * l_ik;
		}
		l_ii = sqrt(l_ii);
		debug("2 Li%di%d = %f\n", i, i, l_ii);

		/* XXX Set */
		skyline_set(aa, jd, i, i, l_ii);
	}

	/*
	for(i = 0; i < n; i++)
	{
		for(j = i+1; j < n; j++)
		{
			skyline_set(aa, jd, i, j, 0.0);
		}
	}
	*/
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
	/* Comprobar con octave:

		chol([ 3 0 -1 0 0; 0 1 0 0 1; -1 0 3 0 0; 0 0 0 5 0; 0 1 0 0 8 ])'
	*/

	double aa_data[] = { 3, 1, -1, 0, 3, 5, 1, 0, 0, 8 };
	double jd_data[] = { 0, 1, 4, 5, 9 };
	int N = 5;
	size_t i, j;

	gsl_vector_view aa_v, jd_v;
	gsl_vector *aa, *jd;

	aa_v = gsl_vector_view_array(aa_data, ARRAY_SIZE(aa_data));
	jd_v = gsl_vector_view_array(jd_data, ARRAY_SIZE(jd_data));

	aa = &aa_v.vector;
	jd = &jd_v.vector;

	gsl_vector *st = skyline_st(jd);
	/*
	skyline_get(aa, jd, 0, 0);
	skyline_set(aa, jd, 0, 0, 2);
	skyline_get(aa, jd, 0, 0);
	*/
	chol_sky(aa, jd, st, N);

	printf("aa = ");
	for(i=0; i<aa->size; i++)
	{
		printf("% 5f,", gsl_vector_get(aa, i));
	}
	printf("\n");

	return 0;
}
