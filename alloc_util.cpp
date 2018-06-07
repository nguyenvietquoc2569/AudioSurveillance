
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include "alloc_util.h"


char *G_malloc(int n)
{
    char *b;

    b = (char*) malloc((unsigned)n);
    if (b || !n) return(b);

    fprintf (stderr, "Out Of Memory\n");
    exit(1);
}


char *G_calloc(int n,int m)
{
    char *b;

    b =(char*) calloc((unsigned)n,(unsigned)m);
    if (b || !n || !m) return(b);

    fprintf (stderr, "Out Of Memory\n");
    exit(1);
}


char *G_realloc(char *b,int n)
{
    if (b == NULL) b = (char*) malloc ((unsigned)n);
    else b = (char*) realloc(b, (unsigned)n);
    if (b || !n) return(b);

    fprintf (stderr, "Out Of Memory\n");
    exit(1);
}

void G_dealloc(char *b)
{
    free( b );
} 

double *G_alloc_vector(int n)
{
    return (double *) G_calloc (n, sizeof(double));
}


double **G_alloc_matrix(int rows,int cols)
{
    double **m;
    int i;

    m = (double **) G_calloc (rows, sizeof(double *));
    m[0] = (double *) G_calloc (rows*cols, sizeof(double));
    for (i = 1; i < rows; i++)
	m[i] = m[i-1] + cols;
    return m;
}


void G_free_vector(double *v)
{
    if(v!=NULL) free ((char *)v);
}


void G_free_matrix(double **m)
{
    if(m!=NULL) {
      free ((char *)(m[0]));
      free ((char *)m);
    }
}


int *G_alloc_ivector(int n)
{
    return (int *) G_calloc (n, sizeof(int));
}


int **G_alloc_imatrix(int rows,int cols)
{
    int **m;
    int i;

    m = (int **) G_calloc (rows, sizeof(int *));
    m[0] = (int *) G_calloc (rows*cols, sizeof(int));
    for (i = 1; i < rows; i++)
	m[i] = m[i-1] + cols;
    return m;
}


void G_free_ivector(int *v)
{
    if(v!=NULL) {
      free ((char *)v);
    }
}


void G_free_imatrix(int **m)
{
    if(m!=NULL) {
      free ((char *)(m[0]));
      free ((char *)m);
    }
}
