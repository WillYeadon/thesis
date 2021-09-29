volScalarField thermalDiff(mixture.k()/mixture.rhoCp());
volScalarField RayleighNumber(alpha2*charLength3*(beta*(T - meltingT)*9.81)/(mixture.nu()*thermalDiff));

/*
volScalarField one(mixture.k());
volScalarField two(mixture.rhoCp());
volScalarField three(mixture.rho());
volScalarField four(mixture.nu());

Info << "k value \t" << one.weightedAverage(mesh.V()).value() << endl;
Info << "rhoCp value \t" << two.weightedAverage(mesh.V()).value() << endl;
Info << "rho value \t" << three.weightedAverage(mesh.V()).value() << endl;
Info << "thermalDiff value \t" << thermalDiff.weightedAverage(mesh.V()).value() << endl;
Info << "nu value \t" << four.weightedAverage(mesh.V()).value() << endl;
Info << "beta value \t" << beta.value() << endl;
Info << "meltingT value \t" << meltingT.value() << endl;
Info << "T value \t" << T.weightedAverage(mesh.V()).value() << endl;
Info << "RayleighNumber value \t" << RayleighNumber.weightedAverage(mesh.V()).value() << endl;
*/

double r_total = 0;
int r_count = 0;

forAll(alpha2, i) 
{
	if (alpha2[i] > 0.5)
	{
		r_total += RayleighNumber[i];
		r_count += 1;
	}
}

if (r_count != 0)
	r_total /= r_count;

double a_total = 0;
int a_count = 0;

volScalarField boundedAlpha2 = min(max(alpha2, scalar(0)), scalar(1));

forAll(boundedAlpha2, i)
{
	a_total += boundedAlpha2[i];
	a_count += 1;
}

if (a_count != 0)
	a_total /= a_count;

if(runTime.write())
{
	output.open("RayleighNumber.txt", std::ios_base::app);
	output << runTime.value() << "\t" << r_total << "\t" << a_total << "\n";
	output.close();
}
