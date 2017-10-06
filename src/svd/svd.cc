#include <gsl/gsl_linalg.h>

int svd(double **a, int m, int n, double *s, double **vt) {
  int err;
  gsl_matrix u;
  gsl_matrix vt1;
  gsl_vector s1;
  gsl_vector *work;

  //fill GSL matrices and vectors with pointers to input arguments:
  //(this is more efficient than allocating them specially and copying
  //the data even though it breaks the abstraction)
  u.size1=m;
  u.size2=n;
  u.tda=n;
  u.data=a[0];
  //is this part necessary? lets test...
  u.block=(gsl_block *) a[0];
  u.owner=0;

  s1.size=n;
  s1.stride=1;
  s1.data=s;
  s1.block=(gsl_block *) s;
  s1.owner=0;

  vt1.size1=n;
  vt1.size2=n;
  vt1.tda=n;
  vt1.data=vt[0];
  vt1.block=(gsl_block *) vt[0];
  vt1.owner=0;

  work=gsl_vector_alloc(n);

  printf("Performing SVD:\n");
  err=gsl_linalg_SV_decomp(&u, &vt1, &s1, work);

  gsl_vector_free(work);

  return err;

}

