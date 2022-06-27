objects = getCommonSeq/ComSeq.o getCommonSeq/ClusterToFeature.o 
ARMA_INCLUDE_FLAG = -I  LinearAlgebra/include
CXXFLAGS = $(ARMA_INCLUDE_FLAG)
LIB_FLAGS = -lblas -llapack   
CC = g++  -std=c++11  -g   
all:    prepareall
	$(CC)$(CXXFLAGS)  -o   Cluster     main_Cluster.cpp   $(objects) $(LIB_FLAGS)
	$(CC)$(CXXFLAGS)  -o   FeatureGen  main_FeatureGen.cpp  $(objects) $(LIB_FLAGS)
prepareall:    subsystem
subsystem:
	$(MAKE) -C getCommonSeq/
clean :  cleansub   
	rm     Cluster
	rm     FeatureGen
cleansub :
	rm  $(objects)

