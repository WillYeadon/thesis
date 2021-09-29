
fullWeld = min(max((alpha2Accumulation + alpha2), 0.0), 1.0);	

float interfaceCells = 0.0;
float weldCells = 0.0;

forAll(alpha2, i)
{
	if (interfaceCapture[i] == 1.0) {interfaceCells += 1;}
	if (interfaceCapture[i] == 1.0 && fullWeld[i] >= 0.95) {weldCells += 1;}
}
if (interfaceCells == 0) {interfaceCells = 1.0;}

//Info << "A " << interfaceCells << "B " << weldCells << endl;
float percentage = (weldCells) / (interfaceCells);
//Info << "percentage " << percentage << endl;

if (percentage >= 0.9)
{
	if (globalPenn == 0)
	{
		int countTwo = 0;
		while (countTwo < 50)
		{
			Info << "P\tE\tN\tN\tE\tT\tR\tA\tT\tE\tD" << endl;
			countTwo += 1;
		}
		globalPenn = 1;
	}
}


forAll(alpha2, i)
{
	if (globalCollapse == 0)
	{
		if (airTop ? (collapseCapture[i] == 1.0 && alpha1[i] == 1.0) : (collapseCapture[i] == 1.0 && alpha2[i] == 1.0))
		{
			int count = 0;
			while (count < 50)
			{
				Info << "C\tO\tL\tL\tA\tP\tS\tE" << endl;
				count += 1;
			}
			tubeCaptured = true;
		    alpha1.write();
		    alpha2.write();
		    alpha3.write();
			globalCollapse = 1;
		}
	}
}

forAll(alpha1, i)
{
	if (globalBurst == 0)
	{
		if (airTop ? (burstCapture[i] == 1.0 && alpha2[i] == 1.0) : (burstCapture[i] == 1.0 && alpha1[i] == 1.0))
		{
			int countTwo = 0;
			while (countTwo < 50)
			{
				Info << "B\tU\tR\tS\tT" << endl;
				countTwo += 1;
			}
			tubeCaptured = true;
		    alpha1.write();
		    alpha2.write();
		    alpha3.write();
			globalBurst = 1;
		}
	}	
}

Info << "interfaceCells " << interfaceCells << " weldCells " << weldCells << " percentage " << percentage << endl;
