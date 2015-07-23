void StartValidity(void) {
  char nameDist[] = "./Appendix/Valid-Report/ValidRep.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.txt";
  char charDist[sizeof nameDist+200];

  sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD);
  validity = fopen(charDist,"w");

	fprintf(validity,"Timestep\tCorrelation Coefficient\n");
}

void CloseValidity(void) {
  fclose(validity);
}

void PerformValidity(void) {
  double SumAtomDist = 0.0;
	double SumCophDist = 0.0;
	double SumAtomDistSquare = 0.0;
	double SumCophDistSquare = 0.0;
	double SumAtomCophDist	 = 0.0;

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			SumAtomDist += atomDIST[ij][kl];
			SumCophDist += cophDIST[ij][kl];
			SumAtomDistSquare += (atomDIST[ij][kl]*atomDIST[ij][kl]);
			SumCophDistSquare += (cophDIST[ij][kl]*cophDIST[ij][kl]);
			SumAtomCophDist		+= (atomDIST[ij][kl]*cophDIST[ij][kl]);
		}
	}

	CorrCoefficient  = SumAtomCophDist - (SumAtomDist*SumCophDist)/NPAIRDIST;
	CorrCoefficient /= sqrt( (SumAtomDistSquare - (SumAtomDist*SumAtomDist)/NPAIRDIST));
	CorrCoefficient /= sqrt( (SumCophDistSquare - (SumCophDist*SumCophDist)/NPAIRDIST));
	fprintf(validity,"%d\t",timestep);
	fprintf(validity,"%.6lf\t",CorrCoefficient);
	fprintf(validity,"\n");
}

void PerformAgglomerateCluster(void) {
	double act, ref, dif;
	NCLUSTERS = 0;
	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		if (leadSORT[ij] == ij) {
			for (kl = (ij+1); kl < NHELIUM; kl++) {
				if (cophDIST[ij][kl] < CUTOFFRAD) {
					leadSORT[kl] = ij;
				}
			}
		}
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		if (leadSORT[ij] == ij) {
			revAtomName[ij] = NCLUSTERS;	//revAtomName range from 0 to NHELIUM
			clusSeq[NCLUSTERS] = ij;			//clusSeq range from 0 to NCLUSTERS
			for (kl = 0; kl < 3; kl++) {
				clusCOM[NCLUSTERS][kl] += atomCOOR[ij][kl];	//clusCentreofMass range from 0 to NCLUSTERS
			}
			
			NCLUSTERS++;
		} else {
			revAtomName[ij] = revAtomName[leadSORT[ij]];

			for (kl = 0; kl < 3; kl++) {
				act = atomCOOR[ij][kl];
				ref = atomCOOR[leadSORT[ij]][kl];
				dif = act - ref;

				if (dif > 0.5*LENGTH[kl])      	 	 act  = act - LENGTH[kl];
				else if (dif < (-0.5*LENGTH[kl]) ) act  = act + LENGTH[kl];

				clusCOM[revAtomName[leadSORT[ij]]][kl] += act;
			}
		}
	}

					
	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		if (clusSIZE[ij] > 0) {
			for (kl = (ij+1); kl < NHELIUM; kl++) {
				if (cophDIST[ij][kl] < CUTOFFRAD) {			
					clusSIZE[ij] += 1;
					clusSIZE[kl]  = 0;
				}
			}
		}
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		if (clusSIZE[ij] > 0) {
			clusPOPU[clusSIZE[ij]-1] += 1;
			clusPOPU[0] -= clusSIZE[ij];
		}
	}

	for (ij = 0; ij < NCLUSTERS; ij++) {
		for (kl = 0; kl < 3; kl++) {
			clusCOM[ij][kl] = clusCOM[ij][kl] / (double)clusSIZE[clusSeq[ij]];
		}
		
	}
}

void CalculateRadiusCluster(void) {
	int nclust = 0, clustid, sizeid, clus_size, member, member2, idx, leader;
	double triside[3], sum1, sum2, sum3, sum4, denom, volumecluster, radius, length, halfperim;
	double PI  = 4.0*atan(1.0);
	double act[3], ref[3], dif, trimer[3][3];
	double side_a, side_b, side_c;
	
	FILE *opentempfile, *pipe;

	clusMEMBER  = (	  int **)malloc( NCLUSTERS*sizeof( int *) );
	clusCNT     = (	  int  *)malloc( NCLUSTERS*sizeof( int ) );
	clusRAD			= ( double *)malloc( NCLUSTERS*sizeof( double ) );
	for (ij = 0; ij < NCLUSTERS; ij++) {
		clusCNT[ij] = 0;
		clusMEMBER[ij]  = ( int *)malloc(clusSIZE[clusSeq[ij]]*sizeof( int  ) );
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		clustid = revAtomName[ij];
		sizeid  = clusCNT[clustid];
		clusMEMBER[clustid][sizeid] = ij;
		clusCNT[clustid]++;
	}

	for (ij = 0; ij < NCLUSTERS; ij++) {
		leader    = clusSeq[ij];
		clus_size = clusSIZE[leader];
//		printf("clusid = %d, leader = %d, clussize = %d, ", ij, leader, clus_size);
		if (clus_size >= 4) {

			pipe = popen("python","w");
			fprintf(pipe, "import numpy\n");
			fprintf(pipe, "from scipy.spatial import Delaunay\n");
			fprintf(pipe, "from scipy.spatial import ConvexHull\n");
			fprintf(pipe, "from numpy import *\n");
			fprintf(pipe, "from clusterradius_3dconvexhull import tetrahedron_volume\n");
			fprintf(pipe, "from clusterradius_3dconvexhull import convex_hull_volume\n");
			fprintf(pipe, "from clusterradius_3dconvexhull import convex_hull_volume_bis\n");

			fprintf(pipe, "data = [");
			for (kl = 0; kl < clus_size; kl++) {
				member = clusMEMBER[ij][kl];
				for (idx = 0; idx < 3; idx++) {
					act[idx] = atomCOOR[member][idx];
					ref[idx] = atomCOOR[leader][idx];
					dif = act[idx] - ref[idx];
					if 			(dif >   0.5*LENGTH[idx]  ) act[idx]  = act[idx] - LENGTH[idx];
					else if (dif < (-0.5*LENGTH[idx]) ) act[idx]  = act[idx] + LENGTH[idx];
				}
				fprintf(pipe, "(%.6lf,%.4lf,%.6lf),", act[0], act[1], act[2]);
			}
			fprintf(pipe, "]\n");
			fprintf(pipe, "points = numpy.array(data)\n");
			fprintf(pipe, "res = convex_hull_volume(points)\n");
			fprintf(pipe, "f = open(\'__temp_python_file.txt\',\'w\')\n");
			fprintf(pipe, "f.write(repr(res))\n");
			fprintf(pipe, "f.close()\n");
			pclose(pipe);

			opentempfile = fopen("__temp_python_file.txt","r");
			fscanf(opentempfile, "%lf", &volumecluster);
			fclose(opentempfile);

			radius = pow(3.0*volumecluster/(4.0*PI), 1.0/3.0);
		}

		if (clus_size == 3) {
			for (kl = 0; kl < 3; kl++) {
				length = 0.0;
				member  = clusMEMBER[ij][kl];
				for (idx = 0; idx < 3; idx++) {
					act[idx] = atomCOOR[member][idx];
					ref[idx] = atomCOOR[leader][idx];
					dif = act[idx] - ref[idx];
					if 			(dif >   0.5*LENGTH[idx]  ) act[idx]  = act[idx] - LENGTH[idx];
					else if (dif < (-0.5*LENGTH[idx]) ) act[idx]  = act[idx] + LENGTH[idx];
					trimer[kl][idx] = act[idx];
				}
			}

			for (kl = 0; kl < 3; kl++) {
				length = 0.0;
				for (idx = 0; idx < 3; idx++) {
					length += (trimer[kl][idx] - trimer[(kl+1)%3][idx])*(trimer[kl][idx] - trimer[(kl+1)%3][idx]);
				}
				triside[kl] = sqrt(length);
			}
			side_a = triside[0]; side_b = triside[1]; side_c = triside[2];

			sum1 = side_a + side_b + side_c;
			sum2 = side_b + side_c - side_a;
			sum3 = side_c + side_a - side_b;
			sum4 = side_a + side_b - side_c;

			halfperim = 0.5*sum1;
			
			radius = sqrt( (sum2*sum3*sum4) / (4.0*sum1) );

			denom = sqrt(sum1*sum2*sum3*sum4);
		//	radius = side_a*side_b*side_c / denom;
		}

		if (clus_size == 2) {
			length = 0.0;
			member  = clusMEMBER[ij][0];
			member2 = clusMEMBER[ij][1];
			for (idx = 0; idx < 3; idx++) {
				dif = atomCOOR[member][idx] - atomCOOR[member2][idx];
				if 			(dif >   0.5*LENGTH[idx]  ) dif  = dif - LENGTH[idx];
				else if (dif < (-0.5*LENGTH[idx]) ) dif  = dif + LENGTH[idx];				
				length += (dif*dif);
			}
			radius = sqrt(length);
		}

		if (clus_size == 1) {
			radius = 0.0;
		}

		clusRAD[ij] = radius;

	}

	free(clusMEMBER);
	free(clusCNT);
}

void BinRadiusCluster(void) {
	double ratiorad = maxbinRadius / intbinRadius;
	double invint = 1.0/ intbinRadius;
	double actrad;
	int binidx, clussize;

	NBINRAD = (int)ratiorad + 1;

	cntBinRAD     = (	  int  *)malloc( NBINRAD*sizeof( int ) );
	sizeAvRAD     = (	  int  *)malloc( NBINRAD*sizeof( int ) );
	radAvSIZE			= (	double *)malloc( NHELIUM*sizeof( double ) );
	for (ij = 0; ij < NBINRAD; ij++) {
		cntBinRAD[ij] = 0;
		sizeAvRAD[ij] = 0;
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		radAvSIZE[ij] = 0.0;
	}

	for (ij = 0; ij < NCLUSTERS; ij++) {
		actrad = floor( clusRAD[ij] * invint );
		clussize = clusSIZE[clusSeq[ij]];
		binidx = (int)actrad;
		cntBinRAD[binidx] += 1;
		sizeAvRAD[binidx] += clussize;
		radAvSIZE[clussize-1] += clusRAD[ij];

	}
		
}
		

void PrintDistMatrix(void) {
  char nameDist[] = "./Appendix/Dist-Matrix/DistanceMatrix.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			fprintf(distmatrix,"%4d\t%4d\t%.8lf\n", ij, kl, atomDIST[ij][kl]);
		}
	}

	fclose(distmatrix);
}

void PrintCophMatrix(void) {
  char nameDist[] = "./Appendix/Coph-Matrix/CophMatrix.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			fprintf(distmatrix,"%4d\t%4d\t%.8lf\t%.8lf\n", ij, kl, cophDIST[ij][kl], cophSIGN[ij][kl]);
		}
	}

	fclose(distmatrix);
}

void PrintRadiusCluster(void) {
  char nameDist[] = "./Clust-Radius/ClustRadius.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	fprintf(distmatrix,"ClusterID\tClustSize\tClustRad\n");
	for (ij = 0; ij < NCLUSTERS; ij++) {
		fprintf(distmatrix,"%d\t", ij);
		fprintf(distmatrix,"%d\t", clusSIZE[clusSeq[ij]]);
		fprintf(distmatrix,"%.6e\t", clusRAD[ij]);
		fprintf(distmatrix,"\n");
	}
	fclose(distmatrix);
}

void PrintBinRadiusCluster(void) {
  char nameDist[] = "./Bin-Radius/BinRadius.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	double aveSize = 0.0;
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	fprintf(distmatrix,"Radius\tNumber of Clusters\tAveSize\n");
	for (ij = 0; ij < NBINRAD; ij++) {
		aveSize = (cntBinRAD[ij] != 0) ? (double)sizeAvRAD[ij] / (double)cntBinRAD[ij] : 0.0;
		fprintf(distmatrix,"%d\t", ij);
		fprintf(distmatrix,"%d\t", cntBinRAD[ij]);
		fprintf(distmatrix,"%.6lf\t", aveSize);
		fprintf(distmatrix,"\n");
	}
	fclose(distmatrix);
}


void PrintClusterPopulation(void) {
  char nameDist[] = "./Clust-Distribution/ClustDistribution.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	double avrad;
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");
	fprintf(distmatrix,"Size\tTotalPopulation\tAverage Radius\n");

	for (ij = 0; ij < NHELIUM; ij++) {
		avrad = (clusPOPU[ij] != 0) ? (radAvSIZE[ij] / (double)clusPOPU[ij]) : 0.0; 
		fprintf(distmatrix,"%4d\t", ij+1);
		fprintf(distmatrix," %d\t", clusPOPU[ij]);
		fprintf(distmatrix,"%.4e\t", avrad);
		fprintf(distmatrix,"\n");
	}

	fclose(distmatrix);
}

void PrintAtomSize(void) {
  char nameDist[] = "./Appendix/Atom-Size/ATOMSIZE.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t%4d\n", ij,atomSIZE[ij]);
	}

	fclose(distmatrix);
}

void PrintAtomCoordinate(void) {
  char nameDist[] = "./Appendix/Atom-Coordinate/ATOMCoordinate.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t%.8lf\t%.8lf\t%.8lf\n", ij,atomCOOR[ij][0],atomCOOR[ij][1],atomCOOR[ij][2]);
	}

	fclose(distmatrix);
}

void PrintAtomLead(void) {
  char nameDist[] = "./Appendix/Atom-Lead/ATOMLEAD.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	int yes;
	FILE *distmatrix;

	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");
	fprintf(distmatrix,"Index\tX\tY\tZ\tLeadIndx\tSIZE\tLeader?\n");
	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t",ij);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][0]);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][1]);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][2]);
		fprintf(distmatrix,"%4d\t",leadSORT[ij]);
		fprintf(distmatrix,"%4d\t",clusSIZE[leadSORT[ij]]);
		yes = (leadSORT[ij] == ij) ? 1 : 0;
		fprintf(distmatrix,"%4d\t",yes);
		fprintf(distmatrix,"\n");
	}

	fclose(distmatrix);
}

void PrintCFGClustCOMHeader(void) {
  fprintf(cfgnew,"Number of particles = %d\n", NCLUSTERS);
  fprintf(cfgnew,"A = 1 Angstrom (basic length-scale)\n");
  fprintf(cfgnew,"H0(1,1) = %lf A\n",LENGTH[0]);
  fprintf(cfgnew,"H0(1,2) = 0 A\n");
  fprintf(cfgnew,"H0(1,3) = 0 A\n");
  fprintf(cfgnew,"H0(2,1) = 0 A\n");
  fprintf(cfgnew,"H0(2,2) = %lf A\n",LENGTH[1]);
  fprintf(cfgnew,"H0(2,3) = 0 A\n");
  fprintf(cfgnew,"H0(3,1) = 0 A\n");
  fprintf(cfgnew,"H0(3,2) = 0 A\n");
  fprintf(cfgnew,"H0(3,3) = %lf A\n",LENGTH[2]);
  fprintf(cfgnew,".NO_VELOCITY.\n");
  fprintf(cfgnew,"entry_count = 7\n");
  fprintf(cfgnew,"auxiliary[0] = index\n");
  fprintf(cfgnew,"auxiliary[1] = leader\n");
  fprintf(cfgnew,"auxiliary[2] = size\n");
  fprintf(cfgnew,"auxiliary[3] = radius\n");
  fprintf(cfgnew,"4.008\n");
  fprintf(cfgnew,"Ga\n");
}

void PrintCFGClusterCentreOfMass(void) {
  char nameDist[] = "./Appendix/CFG-CentOfMass/CFG.COM.Clust.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	sprintf(charDist,nameDist,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  cfgnew = fopen(charDist,"w");

	PrintCFGClustCOMHeader();
	for (ij = 0; ij < NCLUSTERS; ij++) {
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][0]/LENGTH[0]);
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][1]/LENGTH[1]);
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][2]/LENGTH[2]);

		fprintf(cfgnew,"%d\t",ij);
		fprintf(cfgnew,"%d\t",clusSeq[ij]);
		fprintf(cfgnew,"%d\t",clusSIZE[clusSeq[ij]]);
//		fprintf(cfgnew,"0.0\t");
		fprintf(cfgnew,"%.4e\t",clusRAD[ij]);
		fprintf(cfgnew,"\n");	
	}	

	fclose(cfgnew);
}

void PrintCFGHeliumHeader(void) {
  fprintf(heliumcfg,"Number of particles = %d\n", NHELIUM);
  fprintf(heliumcfg,"A = 1 Angstrom (basic length-scale)\n");
  fprintf(heliumcfg,"H0(1,1) = %lf A\n",LENGTH[0]);
  fprintf(heliumcfg,"H0(1,2) = 0 A\n");
  fprintf(heliumcfg,"H0(1,3) = 0 A\n");
  fprintf(heliumcfg,"H0(2,1) = 0 A\n");
  fprintf(heliumcfg,"H0(2,2) = %lf A\n",LENGTH[1]);
  fprintf(heliumcfg,"H0(2,3) = 0 A\n");
  fprintf(heliumcfg,"H0(3,1) = 0 A\n");
  fprintf(heliumcfg,"H0(3,2) = 0 A\n");
  fprintf(heliumcfg,"H0(3,3) = %lf A\n",LENGTH[2]);
  fprintf(heliumcfg,".NO_VELOCITY.\n");
  fprintf(heliumcfg,"entry_count = 8\n");
	fprintf(heliumcfg, "auxiliary[0] = clustersize\n");
	fprintf(heliumcfg, "auxiliary[1] = clusterid\n");
	fprintf(heliumcfg, "auxiliary[2] = radius\n");
  fprintf(heliumcfg,"4.0026\n");
  fprintf(heliumcfg,"Ga\n");
}

void PrintLineHeAtom(int index) {
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][0] / LENGTH[0]);
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][1] / LENGTH[1]);
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][2] / LENGTH[2]);
	fprintf(heliumcfg,"%d\t",	clusSIZE[leadSORT[index]]);
	fprintf(heliumcfg,"%d\t", revAtomName[leadSORT[index]]);
	fprintf(heliumcfg,"%.6lf\n",	clusRAD[revAtomName[leadSORT[index]]] );
}		

void PrintCFGCluster(void) {
	int i;
  char nameHeCFG[]  = "./CFG-cluster/CFG.Cluster.%s%d.T%dK.He-%.1lfppt.CUTOFFRAD%.0lf.%d.cfg";
  char charHeCFG[sizeof nameHeCFG+200];
	sprintf(charHeCFG,nameHeCFG,NAMECRYSTAL,DEGREE,TEMP,SolConc,CUTOFFRAD,timestep);
  heliumcfg = fopen(charHeCFG,"w");

	PrintCFGHeliumHeader();

	for (i = 0; i < NHELIUM; i++) {
		PrintLineHeAtom(i);
	}

	fclose(heliumcfg);

}

