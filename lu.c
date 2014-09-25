#include <stdio.h>
#include <stdlib.h>
//#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

//Funciones inline para get y set
#define HAVE_INLINE 1

inline double sum_s(gsl_matrix *A, size_t rows, size_t cols, int i, int j)
{
	double s = 0.0;
	int k;
	printf(" S=0.0\n"); 
	for(k = 0; k < j; k++)
	{
		double l_ik = gsl_matrix_get(A, i, k);
		double u_kj = gsl_matrix_get(A, k, j);
		s += l_ik * u_kj;
		printf(" S+= l_%d%d * u_%d%d = %1.f * %1.f = %1.f\n", 
		               i,k,     k,j,   l_ik, u_kj,   s);
	}
	return s;
}

inline double sum_u(gsl_matrix *A, size_t rows, size_t cols, int i, int j)
{
	double s = 0.0;
	int k;
	printf(" S=0.0\n"); 
	for(k = 0; k <= i-1; k++)
	{
		double l_ik = gsl_matrix_get(A, i, k);
		double u_kj = gsl_matrix_get(A, k, j);
		s += l_ik * u_kj;
		printf(" S+= l_%d%d * u_%d%d = %1.f * %1.f = %1.f\n", 
		               i,k,     k,j,   l_ik, u_kj,   s);
	}
	return s;
}

void lu(gsl_matrix *A, size_t rows, size_t cols)
{
	printf("%dx%d\n", rows, cols);
	size_t i = 1, j = 0, k;
	double s;
	for(i = 1; i < rows; i++)
	{
		printf("\n-- i=%d\n", i);
		for(j = 0; j < i; j++)
		{
			s = sum_s(A, rows, cols, i, j);
			double a_ij = gsl_matrix_get(A, i, j);
			double a_jj = gsl_matrix_get(A, j, j);
			
			double l_ij = gsl_matrix_get(A, i, j) - s;
			l_ij /= gsl_matrix_get(A, j, j);
			
			printf("L%d%d\ts=%1.f\ta_%d%d=%1.f  \ta_%d%d=%1.f  \tl_%d%d=%1.f\n", 
			          i,j,     s,       i,j, a_ij,    j,j, a_jj,  i,j, l_ij);
			
			/* XXX SET */
			gsl_matrix_set(A, i, j, l_ij);
		}
		
		s = sum_s(A, rows, cols, i, i);
		double a_ii = gsl_matrix_get(A, i, i);
		double u_ii = a_ii - s;
		
		printf("U%d%d\ts=%1.f\ta_%d%d=%1.f  \tu_%d%d=%1.f\n", 
		          i,i,    s,      i,i, a_ii,     i,i, u_ii);
			

		/* XXX SET */
		gsl_matrix_set(A, i, i, u_ii);

		for(j = i+1; j < rows; j++)
		{
			s = sum_u(A, rows, cols, i, j);
			double a_ij = gsl_matrix_get(A, i, j);
			double u_ij = a_ij - s;
			
			printf("U%d%d\ts=%1.f\ta_%d%d=%1.f  \tu_%d%d=%1.f\n", 
		                  i,j,    s,      i,j, a_ij,     i,j, u_ij);

			/* XXX SET */
			gsl_matrix_set(A, i, j, u_ij);
		}
	}
}

int main(int argc, char *argv[])
{
	gsl_matrix *A;
	A = gsl_matrix_alloc(4, 4);
	if(A==NULL)
	{
		perror("gsl_matrix_alloc");
		exit(1);
	}

	double m1[4][4] =
	{
		{  2.0, -1.0, 2.0, 3.0 },
		{  2.0, 3.0, 3.0, 5.0 },
		{ -2, 9, 3, 2 },
		{  8,-8,10,15 }
	};
	int i,j;
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("%4.f ", m1[i][j]);
			gsl_matrix_set(A, i, j, m1[i][j]);
		}
		printf("\n");
	}
	lu(A, 4, 4);
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			printf("%4.f ", gsl_matrix_get(A, i, j));
		}
		printf("\n");
	}
	gsl_matrix_free(A);
	return 0;
}
