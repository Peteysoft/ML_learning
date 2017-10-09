#include <stdio.h>
#include <stdlib.h>

//subroutine header for performing singular value decomposition:
#include "svd.h"

//subroutine header for performing cluster analysis:
#include "cluster.h"

//maximum number of clusters:
#define MAX_CLUSTER 10

int main(int argc, char **argv) {
  char *infile;            //input file
  char *outfile;            //output file
  FILE *fs;                //file pointer
  double **a;            //matrix of training data/U
  int m;                //number of rows in matrix
  int n;                //number of columns in matrix
  int nsv;                //number of singular values

  if (argc!=4) {
    printf("syntax: cluster_svd nsv train centers\n");
    printf("  where:\n");
    printf("nsv      = number of singular values (0 = use untransformed data)\n");
    printf("infile   = ASCII input file containing training data\n");
    printf("output   = ASCII output file containing cluster centers\n");
    printf("\n");
    printf("file format:\n");
    printf("- one line header containing number of rows and number of columns\n");
    printf("- row major list of each matrix element\n");
    exit(1);
  }

  if (sscanf(argv[1], "%d", &nsv)!=1) {
    fprintf(stderr, "Error parsing first command line argument\n");
    exit(1);
  }
  infile=argv[2];
  outfile=argv[3];

  fs=fopen(infile, "r");
  if (fs==NULL) {
    fprintf(stderr, "Error opening input file, %s\n", infile);
    exit(1);
  }

  if (fscanf(fs, "%d %d", &m, &n)!=2) {
    fprintf(stderr, "Format error in input file: %s line 0", infile);
    exit(1);
  }
  if (nsv>n) {
    fprintf(stderr, "Command line parameter nsv=%d out of range\n", nsv);
    exit(1);
  }

  a=new double *[m];
  a[0]=new double[m*n];
  for (int i=1; i<m; i++) a[i]=a[0]+i*n;
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      if (fscanf(fs, "%lg", a[i]+j)!=1) {
        fprintf(stderr, "Format error in input file, %s, line %d\n", infile, i);
        exit(1);
      }
    }
  }

  fclose(fs);

  /*********************** end of first part *********************************/

  double *ave;
  double *s;                 //singular values
  double **vt;                //right singular vectors

  //first we calculate and remove the arithmetic means:
  ave=new double[n];
  for (int i=0; i<n; i++) ave[i]=0;
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      ave[j]+=a[i][j];
    }
  }
  for (int i=0; i<n; i++) ave[i]/=m;
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      a[i][j]-=ave[j];
    }
  }

  if (nsv>0) {
    //make space for singular values:
    s=new double[n];

    //make space for right singular vectors:
    vt=new double *[n];
    vt[0]=new double[n*n];
    for (int i=1; i<n; i++) vt[i]=vt[0]+i*n;

    //perform the decomposition:
    int err=svd(a, m, n, s, vt);
    if (err!=0) {
      fprintf(stderr, "Error in svd subroutine\n");
      exit(err);
    }
  }

  /*********************** end of second part ********************************/

  double **mu_p;        //matrix of cluster centers
  int nc;        //number of cluster centers

  //make space for cluster centers:
  mu_p=new double *[MAX_CLUSTER];
  mu_p[0]=new double[MAX_CLUSTER*n];
  for (int i=1; i<MAX_CLUSTER; i++) mu_p[i]=mu_p[0]+i*n;

  //perform the cluster analysis:
  if (nsv>0) {
    nc=cluster(a, m, nsv, MAX_CLUSTER, mu_p);
  } else {
    nc=cluster(a, m, n, MAX_CLUSTER, mu_p);
  }

  if (nc <= 0) {
    fprintf(stderr, "Cluster algorithm failed");
    exit(-1);
  }

  /*********************** end of third part ********************************/

  double **mu;        //cluster centers in un-transformed coords

  //allocate space for the un-transformed cluster centers:
  mu=new double *[nc];
  mu[0]=new double[nc*n];
  for (int i=1; i<nc; i++) mu[i]=mu[0]+i*n;

  //perform the coordinate transformation:
  for (int i=0; i<nc; i++) {
    for (int j=0; j<n; j++) {
      //mu[i][j]=ave[j];
      mu[i][j]=0;
      if (nsv>0) {
        for (int k=0; k<nsv; k++) mu[i][j]+=vt[k][j]*s[k]*mu_p[i][k];
      } else {
        mu[i][j]+=mu_p[i][j];
      }
    }
  }

  //write the results to a file:
  fs=fopen(outfile, "w");
  if (fs==NULL) {
    fprintf(stderr, "Error opening output file, %s\n", outfile);
    exit(1);
  }

  fprintf(fs, "%d %d\n", nc, n);
  for (int i=0; i<nc; i++) {
    //for (int j=0; j<n; j++) printf("%16.8lg", mu[i][j]);
    for (int j=0; j<n; j++) fprintf(fs, "%16.8lg", mu[i][j]);
    fprintf(fs, "\n");
  }

  fclose(fs);

  //clean up:
  delete [] mu[0];
  delete [] mu;

  delete [] mu_p[0];
  delete [] mu_p;

  delete [] ave;
  delete [] a[0];
  delete [] a;
  if (nsv>0) {
    delete [] s;
    delete [] vt[0];
    delete [] vt;
  }

  return 0;
}


