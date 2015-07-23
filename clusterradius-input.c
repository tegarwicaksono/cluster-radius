	char nameHeInputSingle[] = "./CFG-source/CFG.HeBulk.SIZE%d.Occ-%.0lfppt.T%dK.NRAND1.%d.cfg";
	char charHeInputSingle[sizeof nameHeInputSingle+200];

	char nameHeInputBoundary[] 	= "./CFG-source/CFG.ASYM.FLAT.%dDEG.He-%.0lfppt.T%dK.%d.cfg";
	char charHeInputBoundary[sizeof nameHeInputBoundary+200];

  char skip;
	int CRYSTALTYPE;
	int NHELIUM = 500; //In absolute number
	int	DEGREE	= 5;
  int TEMP		= 1000;
  int SIGMA   = 3;
	int NIRON		= 1;
	int NLCAND	= 0;
	int NRCAND  = 0;
	int NMCAND	= 0;
	int	NCAND		= 0;
	int NVERTEX = 1;
  int LINEFE  = 15;
	int LINEHESINGLE		= 17;
	int LINEHEBOUNDARY	= 16;
	int ENTRYHESINGLE 	= 2;
	int ENTRYHEBOUNDARY = 1;

  int inistep  =   100000;
  int finstep  =  5000000;
	int intstep	 =   100000;
	int	NPAIRDIST, NCLUSTERS;
	int	NHELIUMBULK, NHELIUMBULKOUT, NHELIUMEDGE, NHELIUMLOOP, NHELIUMMINUSONE;

  int    	step, atom, timestep, c, lines = 0;
	int idMAXLEFT, idMINLEFT, idMAXRIGHT, idMINRIGHT;

	double	int1lo = 0.31;
	double	int1hi = 0.19;
	double	int2lo = 0.81;
	double	int2hi = 0.69;
	double	cutoffPOTEN = 1.00;
	double 	SolConc;
//	int			NATOMS, NATOMSM;

  double 	x, y, z, ix, iy, iz;
  double 	lx = 0.0, ly = 0.0, lz = 0.0;

	double	CUTOFFDISTANCE = 10.0;	//in Angstrom
	double 	CorrCoefficient;
	double 	CUTOFFRAD = 120.0; //In Angstrom
	double 	CUTOFFRADSQ;

	double	TOPLEFT[3], TOPRIGHT[3], BOTLEFT[3], BOTRIGHT[3];
	double	LEFTEDGE[2] = {0.138, 0.152};		//For the left edge of U-loop
	double	RIGHTEDGE[2] = {0.850, 0.864};	//For the right edge of U-loop
	double	LEFTLOC  = 0.140;
	double	RIGHTLOC = 0.857;
	double	upperBorder  = 1.0;

  double 	LENGTH[3], HALFLENGTH[3], ATOMCOOR[3], ACTCOOR[3], ENTRYCFG;

	int 		target1, target2;

	int 		*atomNAME, *atomLEAD, *atomSIZE, *leadSORT, *revAtomName, **clusMEMBER;
	int 		*labelBULK, *labelEDGE, *labelLOOP, *labelBULKOUT;
	int 		*clusPOPU, *clusSIZE,  *clusSeq, *clusCNT;
	double 	**atomCOOR, **atomDIST, **cophDIST, **cophSIGN, **clusCOM, *clusRAD;
	double	*entryCFG, *poteng;
	double	*LEFTCAND, *RIGHTCAND, **LOOPCAND, *ZLEFT, *ZRIGHT, **VXCAND;
	double	**VXPOINT, **HULLEDGE, **distHelium, **ClosestEdge;
	int			*LabelHeAtoms;

	double  maxbinRadius = 4.00; //in Angstrom
	double	intbinRadius = 0.05; //in Angstrom
	int			NBINRAD;
	int			*cntBinRAD, *sizeAvRAD;
	double	*radAvSIZE;

	char		NAME1[10] = "BULK.SIZE";
	char		NAME2[10] = "GB.DEGREE";
	char		NAMECRYSTAL[10];
  char 		xa, ya, za;
  char 		ig1[7], ig2[1], ig3[1];
	char		w1[6], w2[2], w3[9], w4[1];
  int  		ch, ij, kl, mn;

  char nameHeEvolSizeBulkInside[] = "./He-evolution/BulkInside/Evol.BulkIn.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeBulkInside[sizeof nameHeEvolSizeBulkInside+500];

  char nameHeEvolSizeBulkOutside[] = "./He-evolution/BulkOutside/Evol.BulkOut.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeBulkOutside[sizeof nameHeEvolSizeBulkOutside+500];

  char nameHeEvolSizeEdge[] = "./He-evolution/Edge/Evol.Edge.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeEdge[sizeof nameHeEvolSizeEdge+500];

  char nameHeEvolSizeLoop[] = "./He-evolution/Loop/Evol.Loop.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeLoop[sizeof nameHeEvolSizeLoop+500];


	FILE *input, *validity, *cfgnew;
	FILE *ironinput, *ironCand, *heliuminput, *ironvert, *evol, *heliumcfg;
	FILE *evolsizebulkin, *evolsizebulkout, *evolsizeedge, *evolsizeloop;


void FindNumberHeliumAtoms(void) {

	if (CRYSTALTYPE == 1) {
		if (SolConc == 1.0) NHELIUM = 250;
		if (SolConc == 5.0) NHELIUM = 1256;
		if (SolConc == 10.0) NHELIUM = 2525;
		if (SolConc == 12.0) NHELIUM = 3036;
	} 

	if (CRYSTALTYPE == 2) {
		if (DEGREE == 201) {
			if (SolConc == 1.0) NHELIUM = 647;
			if (SolConc == 5.0) NHELIUM = 3256;
			if (SolConc == 10.0) NHELIUM = 6543;
		}
		if (DEGREE == 44) {
			if (SolConc == 1.0) NHELIUM = 709;
			if (SolConc == 5.0) NHELIUM = 3560;
			if (SolConc == 10.0) NHELIUM = 7157;
		}
		if (DEGREE == 70) {
			if (SolConc == 1.0) NHELIUM = 648;
			if (SolConc == 5.0) NHELIUM = 3230;
			if (SolConc == 10.0) NHELIUM = 6485;
		}
		if (DEGREE == 90) {
			if (SolConc == 1.0) NHELIUM = 645;
			if (SolConc == 5.0) NHELIUM = 3230;
			if (SolConc == 10.0) NHELIUM = 6485;
		}
	}
}

void AskForInput(void) {
	double ini, fin, interval;
	printf("(1) if single crystal, (2) if bicrystals: ");
	scanf("%d", &CRYSTALTYPE);

	if (CRYSTALTYPE == 1) {
		DEGREE = 50;
		strcpy(NAMECRYSTAL, NAME1);
	}

	if (CRYSTALTYPE == 2) {
		printf("Enter the DEGREE (201, 44, 70, 90): ");
		scanf("%d", &DEGREE);
		strcpy(NAMECRYSTAL, NAME2);
	}
		
	printf("Enter He conc (in parts per thousand): ");
	scanf("%lf", &SolConc);

	printf("Enter Temperature (in K): ");
	scanf("%d", &TEMP);

	printf("Enter the cutoff radius for cluster (in Angstrom^2): ");
	scanf("%lf", &CUTOFFRAD);

	FindNumberHeliumAtoms();

	printf("Enter the timestep of the first index (in ns): ");
	scanf("%lf", &ini);
	ini = ini*1.0e6;
	inistep = (int)ini;
	
	printf("Enter the timestep of the last index (in ns): ");
	scanf("%lf", &fin);
	fin = fin*1.0e6;
	finstep = (int)fin;

	printf("Enter the interval timestep (in ns): ");
	scanf("%lf", &interval);
	interval = interval*1.0e6;
	intstep = (int)interval;
/*
	printf("inistep = %d\n", inistep);
	printf("finstep = %d\n", finstep);
	printf("intstep = %d\n", intstep);
*/
}
