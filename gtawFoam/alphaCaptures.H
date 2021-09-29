    alphaMetalTotal = alpha2 + alpha3;
	alphaMetalTotal = max(min(alphaMetalTotal, scalar(1)), scalar(0));

	Info<< "\n alpha2 Weighted Fraction: " << (alpha2.weightedAverage(mesh.V()).value() / alphaMetalTotal.weightedAverage(mesh.V()).value()) << endl;

	if ((alpha2.weightedAverage(mesh.V()).value() / alphaMetalTotal.weightedAverage(mesh.V()).value()) > alpha2Capture.value())
	{
		alpha2WeightedAv = alpha2;			
	}
	else
	{
		alpha2WeightedAv = 0.0;
	}
