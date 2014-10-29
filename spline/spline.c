#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define NDEBUG 1
#include "dbg.h"

struct spline_t
{
	size_t n;
	double *x;
	double *a, *b, *c, *d;
	double p0, pn;

	/* Private */
	double *h;
	double *alpha;
	double *l;
	double *mu;
	double *z;
};

void spline(struct spline_t *s)
{
	ssize_t i, j;
	size_t n = s->n;
	double *x = s->x, *a = s->a, *b = s->b, *c = s->c, *d = s->d;
	double p0 = s->p0, pn = s->pn;
	double *h = s->h, *alpha = s->alpha, *l = s->l, *mu = s->mu, *z = s->z;

	for(i = 0; i<= n-1; i++)
	{
		h[i] = x[i+1] - x[i];
	}

	alpha[0] = 3*(a[1] - a[0])/h[0] - 3*p0;
	alpha[n] = 3*pn - 3*(a[n] - a[n-1])/h[n-1];

	for(i = 1; i <= n-1; i++)
	{
		alpha[i] = 3/h[i]*(a[i+1] - a[i]) - 3/h[i-1]*(a[i] - a[i-1]);
	}

	/* Factorización de Crout */

	l[0] = 2*h[0];
	mu[0] = 0.5;
	z[0] = alpha[0]/l[0];

	for(i = 1; i<= n-1; i++)
	{
		l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
		mu[i] = h[i]/l[i];
		z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
	}

	l[n] = h[n-1]*(2 - mu[n-1]);
	z[n] = (alpha[n] - h[n-1]*z[n-1])/l[n];
	c[n] = z[n];

	/* Fin de la factorización de Crout */

	for(j = n-1; j >= 0; j--)
	{
		c[j] = z[j] - mu[j]*c[j+1];
		b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3;
		d[j] = (c[j+1] - c[j])/(3*h[j]);
	}
}

double fun(double x)
{
	//return sqrt(x);
	return sin(x) * sin(x) * cos(x);
}

void test()
{
	struct spline_t s;
	size_t i, j, n = 10, f = 5;
	double x, y, Sj;
	char ch;

	s.x = malloc(sizeof(double) * (n+1));
	s.a = malloc(sizeof(double) * (n+1));
	s.b = malloc(sizeof(double) * (n+1));
	s.c = malloc(sizeof(double) * (n+1));
	s.d = malloc(sizeof(double) * (n+1));

	s.h = malloc(sizeof(double) * (n+1));
	s.alpha = malloc(sizeof(double) * (n+1));
	s.l = malloc(sizeof(double) * (n+1));
	s.mu = malloc(sizeof(double) * (n+1));
	s.z = malloc(sizeof(double) * (n+1));
	
	s.n = n;
	
	for(i = 0; i<=n; i++)
	{
		s.x[i] = (double) i;
		s.a[i] = fun((double) i);
	}
	s.pn = 0;
	s.p0 = 1;

	spline(&s);

	for(i = 0; i<n*f; i++)
	{
		if(i%f == 0)
			ch = '*';
		else
			ch = ' ';

		j = i/f;
		x = ((double) i) / f;
		y = fun(x);
		Sj = s.a[j] + s.b[j]*(x - s.x[j]) + s.c[j]*pow((x - s.x[j]),2) +
			s.d[j]*pow((x-s.x[j]),3);

		printf("%cf(%f)=%f, S[%d]=%f, a[%d]=%f, b[%d]=%f, c[%d]=%f, d[%d]=%f\n", 
		        ch,  x,  y,    j,  Sj,     j,s.a[j], j,s.b[j], j,s.c[j], j,s.d[j]);
	}
}

int main(int argc, char *argv[])
{
	test();
	return 0;
}
