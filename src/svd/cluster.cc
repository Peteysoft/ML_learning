#include <gsl/gsl_rng.h>
#include <sys/timeb.h>

#define MAXITER 100

#define LOGGING 1

//calculate the distance squared:
double cluster_metric2(double *v1, double *v2, int n) {
  double diff;
  double d2=0;
  for (int i=0; i<n; i++) {
    diff=v1[i]-v2[i];
    d2+=diff*diff;
  }
  return d2;
}

//k-means:
int cluster(double **x, int m, int n, int nc, double **mu) {
  timeb now;		//for random number generator seed
  gsl_rng *rng;		//random number generator
  double min[n];	//minimum of each feature
  double max[n];	//maximum of each feature
  int code[m];		//which cluster does each point belong to?
  int np[nc];		//number of pointer belonging to each cluster
  double d2;		//distance squared
  double d2min;		//minimum distance to cluster center
  int cmin;		//index of cluster center with min. dist.
  double diff;		//difference between two coords
  double **mu2;		//revised cluster centers
  int finish;		//flag for convergence

  printf("Performing cluster analysis:\n");
  //find the ranges for the data:
  for (int j=0; j<n; j++) {
    min[j]=x[0][j];
    max[j]=x[0][j];
  }
  for (int i=1; i<m; i++) {
    for (int j=0; j<n; j++) {
      if (x[i][j] < min[j]) {
        min[j]=x[i][j];
      } else if (x[i][j] > max[j]) {
        max[j]=x[i][j];
      }
    }
  }

  //so we can initialize the cluster centers with a bunch of random points:
  //initialize and seed the random number generator:
  rng=gsl_rng_alloc(gsl_rng_mt19937);
  ftime(&now);
  gsl_rng_set(rng, (long) now.time + ((long) now.millitm) << 7);
  for (int i=0; i<nc; i++) {
    for (int j=0; j<n; j++) {
      mu[i][j]=min[j]+gsl_rng_uniform(rng)*(max[j]-min[j]);
#ifdef LOGGING
      printf("%20.8g", mu[i][j]);
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
  }
#ifdef LOGGING
  printf("\n");
#endif

  //start k-means algorithm:
  //allocate revised cluster centers
  mu2=new double *[nc];
  mu2[0]=new double [nc*n];
  for (int i=0; i<nc; i++) {
    np[i]=0;
    mu2[i]=mu2[0]+i*n;
  }
  //initialize to out-of-range value to eliminate repeat code:
  for (int i=0; i<m; i++) code[i]=-1;

  //do this until we reach convergence
  //(here we define it as points don't change cluster centers)
  finish=1;
  for (int iter=0; iter<MAXITER; iter++) {
    //initialize to zero:
    for (int i=0; i<nc; i++) {
      np[i]=0;
      for (int j=0; j<n; j++) mu2[i][j]=0;
    }
    //find nearest cluster center to each point
    //add this point to revised cluster center
    for (int i=0; i<m; i++) {
      cmin=0;
      d2min=0;
      d2min=cluster_metric2(x[i], mu[0], n);
      for (int j=1; j<nc; j++) {
        d2=cluster_metric2(x[i], mu[j], n);
        if (d2 < d2min) {
          d2min=d2;
	  cmin=j;
        }
      }
      if (code[i]!=cmin) finish=0;
      code[i]=cmin;
      for (int j=0; j<n; j++) mu2[cmin][j]+=x[i][j];
      np[cmin]++;
    }
    for (int i=0; i<nc; i++) {
      for (int j=0; j<n; j++) {
        mu[i][j]=mu2[i][j]/np[i];
#ifdef LOGGING
        printf("%20.8g", mu[i][j]);
#endif
      }
#ifdef LOGGING
      printf("\n");
#endif
    }
#ifdef LOGGING
    printf("\n");
#endif
    if (finish) break;
  }

  delete [] mu2[0];
  delete [] mu2;

  return nc;
}
