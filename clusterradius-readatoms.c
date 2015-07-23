void ScanForHeliumCoordinate(int index) {
	int i, nentry;
	fscanf(heliuminput,"%lf", &atomCOOR[index][0]);
  fscanf(heliuminput,"%lf", &atomCOOR[index][1]);
  fscanf(heliuminput,"%lf", &atomCOOR[index][2]);


	for (i = 0; i < 3; i++) {
		atomCOOR[index][i] = atomCOOR[index][i]*LENGTH[i];
	}
	if (CRYSTALTYPE == 1) nentry = ENTRYHESINGLE;
	if (CRYSTALTYPE == 2) nentry = ENTRYHEBOUNDARY;

	for (i = 0; i < nentry; i++) {
	  fscanf(heliuminput,"%lf", &ENTRYCFG);
	}
}

void ExtractHeliumCoordinate(void) {
	int i;

	for (i = 0; i < NHELIUM; i++) {
		ScanForHeliumCoordinate(i);
	}

	fclose(heliuminput);
}

int atomFOLL(int i, int j) {
  if (i == j) return 0;
	else 				return 1;
}


