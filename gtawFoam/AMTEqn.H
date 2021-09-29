{ 
	Info << "\n TEqn Disperse \n" << endl;
//	alpha2Accumulation += alpha2;
	alpha1Accumulation += alpha1;
	alpha1Accumulation = min(max(alpha1Accumulation, 0.0), 1.0);	
	alpha1AccumulationNeg = (1.0 - alpha1Accumulation);

	correctedAlpha1Accumulation = alpha1Accumulation - alpha2 - alpha3;
	correctedAlpha1Accumulation = min(max(correctedAlpha1Accumulation, 0.0), 1.0);	
	correctedAlpha1AccumulationNeg = 1.0 - correctedAlpha1Accumulation;

	alpha2Accumulation += max((alpha2.oldTime() - alpha2), scalar(0));
	alpha2Accumulation = min(alpha2Accumulation, 1.0);

	alpha3Accumulation += alpha3;
	alpha3Accumulation = min(alpha3Accumulation, 1.0);	
//	alpha3 = alpha3Accumulation;

	alphaMetalTotal = (alpha2 + alpha3);//*rho0;
	volScalarField alphaTotal = (alpha1 + alpha2 + alpha3);
	
	surfaceScalarField alpha2fAc = min(max(fvc::interpolate(alpha2Accumulation), scalar(0)), scalar(1));

	volScalarField boundedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
	volScalarField boundedAlpha2 = min(max(alpha2, scalar(0)), scalar(1));
	volScalarField boundedAlpha3 = min(max(alpha3, scalar(0)), scalar(1));
	volScalarField boundedAlphaTotal = min(max(alpha1 + alpha2 + alpha3, scalar(0)), scalar(1));

	k_field = k1*boundedAlpha1 + k2*boundedAlpha2 + k3*boundedAlpha3;

	Info << "\n New alphaf terms" << endl;

	surfaceScalarField alpha1f = min(max(fvc::interpolate(boundedAlpha1), scalar(0)), scalar(1));
	surfaceScalarField alpha2f = min(max(fvc::interpolate(boundedAlpha2), scalar(0)), scalar(1));
	surfaceScalarField alpha3f = min(max(fvc::interpolate(boundedAlpha3), scalar(0)), scalar(1));

	Info << "\n New kf terms" << endl;

	kEff = mixture.kf();
//		+ fvc::interpolate( mixture.rho()*mixture.cp()*turbulence->nut() );

    Info << "Applying sauce" << endl;
    qSource = arcSource.applyField(alphaMetalTotal)*qSourceUnits;
	
//    rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2 + alpha3*rho3*cp3;
//    rhoCp = alpha1*rho1*cp1 + mixture.rhoBouss()*cp2 + mixture.rhoExp()*cp3;
    rhoCp = mixture.rhoCp();

	Info << "\n Latent heat \n" << endl;
	
	// LfLimit limits the update of Lf so that it is updated only once
	if (!LfLimit)
//		rhoLf = mixture.rho()*Lf;
		rhoLf = LfSmoothing ? mixture.rho()*Lf : mixture.rhoBoussRho()*Lf;

    // Isothermal so div term dropped unless smoothing.
	// rhoBoussRho will mean that the entire domain is like rho2 density as no alpha multiplication
	
	S_latent = LfSmoothing ? rhoLf*fvc::ddt(alpha2) + alpha2*Lf*fvc::div(rhoPhi) : rhoLf*fvc::ddt(alpha2);// + alpha2*Lf*fvc::div(rhoPhi);
	S_latent *= LfSmoothing ? 1.0 * boundedAlphaTotal : (1.0 - alpha1.oldTime())*(1.0 - alpha1);
	S_latent *= phaseChange ? 1.0 : 0.0;

	sFillerAcc += sFiller;
	sFillerAcc = min(max(sFillerAcc, (0/oneSecond)), (1/oneSecond));  
//	volScalarField alphaMetalTotalFiller = min(max((alphaMetalTotal - sFillerAcc)*boundedAlpha1, scalar(0)), scalar(1)); 
	S_latent *= (1.0 - (sFillerAcc*oneSecond)); 	
	/*
	if (LfSmoothing)
	{
		forAll(S_latent, i)
		{
			if ((alpha1[i] + alpha2[i]) >= 1.0)
				S_latent[i] *= 0.0;
		}
	}
	*/

	forAll(T, i)
	{
		if (sFiller[i] > 0.2)
			T[i] = 1.75*meltingT.value();
	}

	AMqSource = AMarcSource.applyField(alphaMetalTotal)*qSourceUnits;
     
	fvScalarMatrix TEqn
	(
		fvm::ddt(rhoCp, T)
		+ fvm::div(rhoPhiCp, T)
// What is this for?
	    - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
		- fvm::laplacian(kEff, T)
		==
		- S_latent
		+ qSource
		+ AMqSource
//		- sFiller*meltingT*sFillerUnitsCorr
	);

    TEqn.relax();
	Info << "\n TEqn End \n" << endl;
	TEqn.solve();
	Info << "\n TEqn Solved \n" << endl;

/*	
// For H vs T comparison
	Info << "\n New alphaEff terms" << endl;
	surfaceScalarField alphaEffRho
	(
		"alphaEffRho",
		fvc::interpolate(rho)*(mixture.alphaf()
				+ fvc::interpolate( turbulence->nut() ) )
	);

    H = (T - T_0)*mixture.cp();
    volScalarField tempCp(mixture.cp());

    T.correctBoundaryConditions();
    H.correctBoundaryConditions();

	Info << "\n TEqn \n" << endl;
    label nEnergyLoops(readLabel(pimple.dict().lookup("nEnergyLoops")));
    for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
		fvScalarMatrix HEqn
		(
			fvm::ddt(rho, H)
			+ fvm::div(rhoPhi, H)
			- fvc::laplacian(kEff, T)
			- (fvm::laplacian(alphaEffRho, H)
				- fvc::laplacian(alphaEffRho, H))
			==
			- S_latent
		);

		HEqn.solve();   

		T = T_0 + H/tempCp;
	}
*/
}