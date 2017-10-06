#ifndef CLUSTER_H
#define CLUSTER_H

  int				//returns number of cluster centers
    cluster (double ** x,	//training vectors
                int m,		//number of training vectors
                int n,		//dimension of each vector
                int max_nc,	//maximum number of cluster centers
                double **mu);	//returned cluster centers

#endif


