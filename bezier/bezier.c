#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define NDEBUG 1
#include "dbg.h"

struct bezier_t
{
	/* Bezier control points */
	double *x;
	double *y;
	double *bx;
	double *by;
	size_t n;

	/* Number of values */
	size_t v;

	/* Private */
	double *tx;
	double *ty;
};

void bezier(struct bezier_t *b)
{
	double *x = b->x, *y = b->y;
	double *tx = b->tx, *ty = b->ty;
	double *bx = b->bx, *by = b->by;

	size_t n = b->n, v = b-> v;

	double t = 0;
	size_t i, j, k;
	
	for(i = 0; i < v; i++)
	{
		t = (double)i/(v-1);

		memcpy(tx, x, sizeof(double) * n);
		memcpy(ty, y, sizeof(double) * n);

		for(j = 0; j < n-1; j++)
		{
			for(k = 0; k < n-j-1; k++)
			{
				tx[k] = (1 - t) * tx[k] + tx[k+1] * t;
				ty[k] = (1 - t) * ty[k] + ty[k+1] * t;
			}

		}

		bx[i] = tx[0];
		by[i] = ty[0];
	}
}

void test()
{
	struct bezier_t b;
	b.x = (double []){ 1, 3, 4, 8 };
	b.y = (double []){ 1, 3, 0, 2 };
	b.n = 4;
	b.v = 40;
	b.bx = malloc(sizeof(double) * b.v);
	b.by = malloc(sizeof(double) * b.v);
	b.tx = malloc(sizeof(double) * b.n);
	b.ty = malloc(sizeof(double) * b.n);

	bezier(&b);

	size_t i;
	printf("x = [ ");
	for(i = 0; i < b.v-1; i++)
	{
		printf("%f, ", b.bx[i]);
	}
	printf("%f ]\n", b.bx[b.v-1]);

	printf("y = [ ");
	for(i = 0; i < b.v-1; i++)
	{
		printf("%f, ", b.by[i]);
	}
	printf("%f ]\n", b.by[b.v-1]);
}

int main(int argc, char *argv[])
{
	test();
	return 0;
}
