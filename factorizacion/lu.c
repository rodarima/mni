#include <stdio.h>
#include <stdlib.h>

//Funciones inline para get y set
#define HAVE_INLINE

//#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#define NDEBUG 1
#include "dbg.h"



inline double sum_s(gsl_matrix *L, gsl_matrix *U, int i, int j)
{
	double s = 0.0;
	int k;
	debug(" S=0.0\n"); 
	for(k = 0; k < j; k++)
	{
		double l_ik = gsl_matrix_get(L, i, k);
		double u_kj = gsl_matrix_get(U, k, j);
		s += l_ik * u_kj;
		debug(" S+= l_%d%d * u_%d%d = %1.f * %1.f = %1.f\n", 
		               i,k,     k,j,   l_ik, u_kj,   s);
	}
	return s;
}

inline double sum_u(gsl_matrix *L, gsl_matrix *U, int i, int j)
{
	double s = 0.0;
	int k;
	debug(" S=0.0\n"); 
	for(k = 0; k <= i-1; k++)
	{
		double l_ik = gsl_matrix_get(L, i, k);
		double u_kj = gsl_matrix_get(U, k, j);
		s += l_ik * u_kj;
		debug(" S+= l_%d%d * u_%d%d = %1.f * %1.f = %1.f\n", 
		               i,k,     k,j,   l_ik, u_kj,   s);
	}
	return s;
}

void lu(gsl_matrix *A, gsl_matrix *L, gsl_matrix *U, size_t n)
{
	size_t i, j, k;
	double s, a_ij, a_jj, a_ii, l_ij, u_ij, u_ii, u_jj;
	if(A != U)
	{
		for(j = 0; j < n; j++)
		{
			double a_ij = gsl_matrix_get(A, 0, j);
			gsl_matrix_set(U, 0, j, a_ij);
		}
	}
	if(A != L)
	{
		for(i = 0; i < n; i++)
		{
			gsl_matrix_set(L, i, i, 1);
		}
	}
	for(i = 1; i < n; i++)
	{
		debug("\n-- i=%d\n", i);
		for(j = 0; j < i; j++)
		{
			s = sum_s(L, U, i, j);
			a_ij = gsl_matrix_get(A, i, j);
			u_jj = gsl_matrix_get(U, j, j);
			l_ij = (a_ij - s) / u_jj;
			
			debug("L%d%d\ts=%1.f\ta_%d%d=%1.f  \ta_%d%d=%1.f  \tl_%d%d=%1.f\n", 
			          i,j,     s,       i,j, a_ij,    j,j, a_jj,  i,j, l_ij);
			
			/* XXX SET */
			gsl_matrix_set(L, i, j, l_ij);
		}
		
		s = sum_s(L, U, i, i);
		a_ii = gsl_matrix_get(A, i, i);
		u_ii = a_ii - s;
		
		debug("U%d%d\ts=%1.f\ta_%d%d=%1.f  \tu_%d%d=%1.f\n", 
		          i,i,    s,      i,i, a_ii,     i,i, u_ii);
			

		/* XXX SET */
		gsl_matrix_set(U, i, i, u_ii);

		for(j = i+1; j < n; j++)
		{
			s = sum_u(L, U, i, j);
			a_ij = gsl_matrix_get(A, i, j);
			u_ij = a_ij - s;
			
			debug("U%d%d\ts=%1.f\ta_%d%d=%1.f  \tu_%d%d=%1.f\n", 
		                  i,j,    s,      i,j, a_ij,     i,j, u_ij);

			/* XXX SET */
			gsl_matrix_set(U, i, j, u_ij);
		}
	}
}

int main(int argc, char *argv[])
{
	double a_data[] =
	{
		 2,-1, 2, 3,
		 2, 3, 3, 5,
		-2, 9, 3, 2,
		 8,-8,10,15
	};
	gsl_matrix *L = gsl_matrix_calloc(4, 4);
	gsl_matrix *U = gsl_matrix_calloc(4, 4);
	gsl_matrix_view m = gsl_matrix_view_array(a_data, 4, 4);

	lu(&m.matrix, L, U, 4);
	int i,j;
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("%4.f ", gsl_matrix_get(L, i, j));
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("%4.f ", gsl_matrix_get(U, i, j));
		}
		printf("\n");
	}

	return 0;
}
