CPP = g++

GSL_INCLUDE=/usr/local/include
GSL_LIB=/usr/local/lib

CPPFLAGS=-g -I$(GSL_INCLUDE)

LDFLAGS=-L$(GSL_LIB) -lgsl -lgslcblas

cluster_svd: cluster_svd.cc svd.o cluster.o svd.h cluster.h
	g++ $(CPPFLAGS) cluster_svd.cc -o cluster_svd svd.o cluster.o $(LDFLAGS)

svd.o: svd.cc
	g++ -c $(CPPFLAGS) svd.cc

cluster.o: cluster.cc
	g++ -c $(CPPFLAGS) cluster.cc

