void LabelFreeze(void) {
	int lo, hi;

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			lo = (atomLEAD[ij] < atomLEAD[kl]) ? atomLEAD[ij] : atomLEAD[kl];
			hi = (atomLEAD[ij] < atomLEAD[kl]) ? atomLEAD[kl] : atomLEAD[ij];
			if ( (lo == target1) && (hi == target2))
				cophSIGN[ij][kl] = -1.00;
		}
	}
}

void FindSmallestDist(void) {
	double smallest = 1.0E08;
  for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			if (cophSIGN[ij][kl] >= 0.0) {
				if (cophDIST[ij][kl] < smallest) {
					smallest = cophDIST[ij][kl];
					target1  = ij;
					target2  = kl;
				}
			}
		}
	}
}

void CalcClusDist(int id0) {
	int	t0, t2, t1;
	double alph1, alph2, beta, sum, temp = 0.0;

	sum   = atomSIZE[target1] + atomSIZE[target2] + atomSIZE[id0];
  beta	= -1.0*(double)atomSIZE[id0]*cophDIST[target1][target2] / sum;
	temp += beta;

	alph2 = (atomSIZE[id0] + atomSIZE[target2]) / sum;
	if (target2 < id0) 	{t0 = target2; t2 = id0; }
	else								{t0 = id0; t2 = target2; }
	temp += (alph2*cophDIST[t0][t2]);

	alph1 = (atomSIZE[id0] + atomSIZE[target1]) / sum;
	if (target1 < id0) 	{t0 = target1; t1 = id0; }
	else								{t0 = id0; t1 = target1; }
	temp += (alph1*cophDIST[t0][t1]);

	cophDIST[t0][t1] = temp;
}

void UpdateCophDist(void) {
  int lo, hi;
	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			if (cophSIGN[ij][kl] >= 0.0) {
				if 				( atomLEAD[ij] == target1 ) {
					if  		( (atomLEAD[ij] == ij) && (atomLEAD[kl] == kl)) CalcClusDist(kl);
					else 		cophDIST[ij][kl] = (atomLEAD[ij] < atomLEAD[kl]) ? cophDIST[atomLEAD[ij]][atomLEAD[kl]] : cophDIST[atomLEAD[kl]][atomLEAD[ij]];
				} else if ( atomLEAD[kl] == target1 )  { 
					if  		( (atomLEAD[ij] == ij) && (atomLEAD[kl] == kl)) CalcClusDist(ij);
					else cophDIST[ij][kl] = (atomLEAD[ij] < atomLEAD[kl]) ? cophDIST[atomLEAD[ij]][atomLEAD[kl]] : cophDIST[atomLEAD[kl]][atomLEAD[ij]];
				} else if ( (ij == target2) || (atomLEAD[ij] == target2) )		{
					cophDIST[ij][kl] = (target1 < kl) ? cophDIST[target1][kl] : cophDIST[kl][target1];
				} else if ( (kl == target2) || (atomLEAD[kl] == target2) )		{
					cophDIST[ij][kl] = (target1 < ij) ? cophDIST[target1][ij] : cophDIST[ij][target1];
				}
			}
		}
	}
}

void UpdateClusSize(void) {
	int ini2 = atomSIZE[target2];

	atomSIZE[target1] += ini2;
	atomSIZE[target2] -= ini2;
	
	for (ij = 0; ij < NHELIUM; ij++) {
		if (atomLEAD[ij] == target2)
			atomLEAD[ij] = target1;
	}
}

void PerformAnalysis(void) {
	int loop = 0, i;
	do {
		FindSmallestDist();
		LabelFreeze();
		UpdateCophDist();
		UpdateClusSize();
		loop++;
	} while (atomSIZE[0] < NHELIUM);

}
