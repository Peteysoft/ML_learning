#include <gsl/gsl_linalg.h>

int svd(double **a, int m, int n, double *s, double **vt) {
  int err;
  gsl_matrix u;
  gsl_matrix vt1;
  gsl_vector s1;
  gsl_vector *work;
  gsl_block ablock;
  gsl_block vblock;
  gsl_block sblock;

  //fill GSL matrices and vectors with pointers to input arguments:
  //(this is more efficient than allocating them specially and copying
  //the data even though it breaks the abstraction)
  u.size1=m;
  u.size2=n;
  u.tda=n;
  u.data=a[0];
  //is this part necessary? lets test...
  ablock.data=a[0];
  ablock.size=sizeof(double)*m*n;
  //u.block=NULL;
  u.block=&ablock;
  u.owner=0;

  s1.size=n;
  printf("%d\n", n);
  s1.stride=1;
  s1.data=s;
  sblock.data=s;
  sblock.size=sizeof(double)*n;
  //s1.block=NULL;
  s1.block=&sblock;
  s1.owner=0;

  gsl_vector_fprintf(stdout, &s1, "%g");

  vt1.size1=n;
  vt1.size2=n;
  vt1.tda=n;
  vt1.data=vt[0];
  vblock.data=vt[0];
  vblock.size=sizeof(double)*n*n;
  //vt1.block=NULL;
  vt1.block=&vblock;
  vt1.owner=0;

  work=gsl_vector_alloc(n);

  printf("Performing SVD:\n");
  err=gsl_linalg_SV_decomp(&u, &vt1, &s1, work);

  gsl_vector_free(work);

  return err;

}

