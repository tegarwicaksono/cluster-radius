void OpenHeAtomicFile(void) {
	int i;

	if (CRYSTALTYPE  == 1) {
		sprintf(charHeInputSingle,nameHeInputSingle,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputSingle,"r");
    for (i = 0; i < 2; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[0]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputSingle,nameHeInputSingle,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputSingle,"r");
    for (i = 0; i < 6; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[1]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputSingle,nameHeInputSingle,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputSingle,"r");
    for (i = 0; i < 10; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[2]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputSingle,nameHeInputSingle,DEGREE,SolConc,TEMP,inistep);
	  heliuminput = fopen(charHeInputSingle,"r");
	  for (i = 0; i < LINEHESINGLE; i++) {
	    do
	      skip = fgetc(heliuminput);
	    while (skip != '\n');
	  }
	}

	if (CRYSTALTYPE == 2) {
		sprintf(charHeInputBoundary,nameHeInputBoundary,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputBoundary,"r");
    for (i = 0; i < 2; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[0]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputBoundary,nameHeInputBoundary,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputBoundary,"r");
    for (i = 0; i < 6; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[1]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputBoundary,nameHeInputBoundary,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputBoundary,"r");
    for (i = 0; i < 10; i++) {
      do
        skip = fgetc(heliuminput);
      while (skip != '\n');
    }
    fscanf(heliuminput,"%s %s", &ig1, &ig2);
    fscanf(heliuminput,"%lf", &LENGTH[2]);
    fscanf(heliuminput,"%s", &ig3);
    fclose(heliuminput);

		sprintf(charHeInputBoundary,nameHeInputBoundary,DEGREE,SolConc,TEMP,timestep);
	  heliuminput = fopen(charHeInputBoundary,"r");
	  for (i = 0; i < LINEHEBOUNDARY; i++) {
	    do
	      skip = fgetc(heliuminput);
	    while (skip != '\n');
	  }
	}
}

void	BuildArrayClusters(void) {
	NPAIRDIST = NHELIUM*(NHELIUM - 1) / 2;
	CUTOFFRADSQ = CUTOFFRAD*CUTOFFRAD;
	NHELIUMMINUSONE = NHELIUM - 1;

	for (ij = 0; ij < 3; ij++)
		HALFLENGTH[ij] = 0.5*LENGTH[ij];

  atomNAME = (		int  	 *)malloc(	NHELIUM*sizeof(	int	)	);
	for (ij = 0; ij < NHELIUM; ij++)
		atomNAME[ij] = 0;

  atomLEAD = (		int		 *)malloc(	NHELIUM*sizeof(	int	)	);
	for (ij = 0; ij < NHELIUM; ij++)
		atomLEAD[ij] = ij;

	leadSORT = (		int		 *)malloc(	NHELIUM*sizeof(	int ) );
	for (ij = 0; ij < NHELIUM; ij++)
		leadSORT[ij] = ij;
	
	atomCOOR = (		double **)malloc(	NHELIUM*sizeof(	double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		atomCOOR[ij] = ( double *)malloc(3*sizeof(double));
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < 3; kl++)
			atomCOOR[ij][kl] = 0.0;

	atomDIST = (		double **)malloc(	NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		atomDIST[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			atomDIST[ij][kl] = 0.0;

	clusCOM = (		double **)malloc(	NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		clusCOM[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			clusCOM[ij][kl] = 0.0;

	cophDIST = (	  double **)malloc( NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++)
		cophDIST[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );

	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			cophDIST[ij][kl] = 0.0;	

	cophSIGN = (	  double **)malloc( NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++)
		cophSIGN[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			cophSIGN[ij][kl] = 0.0;

  clusPOPU 		= (		int	*) malloc( NHELIUM*sizeof( int ) );
	clusSeq			= (		int *) malloc( NHELIUM*sizeof( int ) );
	revAtomName = (		int	*) malloc( NHELIUM*sizeof( int ) );

	for (ij = 0; ij < NHELIUM; ij++) {
		clusPOPU[ij] = 0;
		clusSeq[ij] = 0;
		revAtomName[ij] = 0;
	}
	clusPOPU[0] = NHELIUM;

	atomSIZE = (		int *) malloc( NHELIUM*sizeof( int ) );
	for (ij = 0; ij < NHELIUM; ij++) {
		atomSIZE[ij] = 1;
	}

	clusSIZE = (		int *) malloc( NHELIUM*sizeof( int ) );
	for (ij = 0; ij < NHELIUM; ij++)
		clusSIZE[ij] = 1;
}



void BuildDistanceMatrix(void) {
	double dist = 0.0, del[3];
  for (ij = 0; ij < NHELIUM; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			dist = 0.0;
			for (mn = 0; mn < 3; mn++) {
				del[mn] = fabs(atomCOOR[ij][mn] - atomCOOR[kl][mn]);
				if (del[mn] > HALFLENGTH[mn])
					del[mn] = LENGTH[mn] - del[mn];
				dist += (del[mn]*del[mn]);
			}
			atomDIST[ij][kl] = dist;
			cophDIST[ij][kl] = atomDIST[ij][kl];
//			if (ij == 0) {
//				if (kl < 50) {
//					printf("ij,kl = %d, %d, atomDIST = %.4E\n", ij, kl, atomDIST[ij][kl]);
//				}
//			}
		}
	}			
}

void FreeArrayVertex(void) {
	free(LEFTCAND);
	free(RIGHTCAND);
	free(LOOPCAND);
	free(ZLEFT);
	free(ZRIGHT);
}

void FreeMalloc(void) {

	free(atomNAME);
	free(atomSIZE);
	free(atomLEAD);
	free(atomCOOR);
	free(atomDIST);
	free(revAtomName);
	free(clusSeq);
	free(clusCOM);
	free(clusRAD);
	free(leadSORT);
	free(cophDIST);
	free(cophSIGN);
	free(clusPOPU);
	free(clusSIZE);

	free(labelBULK);
	free(labelEDGE);
	free(labelLOOP);
	free(labelBULKOUT);

	free(cntBinRAD);
	free(sizeAvRAD);
	free(radAvSIZE);
}
