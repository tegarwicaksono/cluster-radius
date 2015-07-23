#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "clusterradius-input.c"
#include "clusterradius-buildarray.c"
#include "clusterradius-readatoms.c"
#include "clusterradius-performanalysis.c"
#include "clusterradius-performvalidity.c"

int main(void) {
	AskForInput();
  for (timestep = inistep; timestep <= finstep; timestep += intstep) {
		printf("timestep = %10d done\n",timestep);
		OpenHeAtomicFile();
	  BuildArrayClusters();
		ExtractHeliumCoordinate();
		BuildDistanceMatrix();
		PerformAnalysis();
		printf("\tcluster analysis done\n");
		PerformAgglomerateCluster();
		printf("\tcluster agglomerate, NCLUSTERS = %d done\n", NCLUSTERS);
		CalculateRadiusCluster();
		printf("\tcluster radius done\n");
		BinRadiusCluster();
		PrintRadiusCluster();
		PrintBinRadiusCluster();
//		PrintAtomLead();
		PrintCFGClusterCentreOfMass();
		PrintClusterPopulation();
		PrintCFGCluster();
		FreeMalloc();
	}
	return 0;
}

