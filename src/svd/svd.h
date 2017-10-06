#ifndef SVD_H
#define SVD_H

  //subroutine for singular value decomposition:
  int				//returns an error code (0 for success)
	svd (double **a,	//input matrix--replaced by U on output
                int m,		//number of rows
                int n,		//number of columns
                double *s,	//singular values
                double **vt);	//V--right singular vectors

#endif

