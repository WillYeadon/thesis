//    MRF.correctBoundaryVelocity(U);

    boundedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
    boundedAlpha2 = min(max(alpha2, scalar(0)), scalar(1));
    boundedAlpha3 = min(max(alpha3, scalar(0)), scalar(1));

//    surfaceScalarField snv = mixture.nHatfMetal();
    surfaceScalarField snv = mixture.nHatf12();
    surfaceScalarField amg = mixture.nHatMetalMG();
    stFactor = 2*rho/(rho2 + rho3);
    alphaMetalTotal = alpha2 + alpha3;
    alphaMetalTotal = max(min(alphaMetalTotal, scalar(1)), scalar(0));

    S_d = -multiplierBuoy*(sqr(alpha3)/(pow3(1.0 - alpha3) + 1e-3))*unitsBuoy;

    S_dg = outerTubeKill ? -multiplierBuoyG*(sqr(1.0 - coreRegion)/(pow3(1.0 - (1.0 - coreRegion)) + 1e-3))*unitsBuoy 
    : -multiplierBuoyG*(sqr(alpha1)/(pow3(1.0 - alpha1) + 1e-3))*unitsBuoy;// 

//    S_Marg = (stOn ? 1.0 : 0.0)*dsigmadT*(fvc::grad(T) - fvc::reconstruct(snv)*(fvc::reconstruct(snv) & fvc::grad(T)))*mag(fvc::grad(alphaMetalTotal))*(2*rho/ (rho1 + rho2));//*pow(alpha2,2);
    S_Marg = (stOn ? 1.0 : 0.0)*dsigmadT*(fvc::grad(T) - fvc::reconstruct(snv)*(fvc::reconstruct(snv) & fvc::grad(T)))*mag(fvc::grad(alpha2))*(2*rho/ (rho1 + rho2));//*alpha2;
    /*
    forAll(S_Marg, i)
    {
        if (alpha2[i] < 1e-3)
            S_Marg[i] *= 0.0;       
    }
    */
//    surfaceScalarField snv = mixture.nHatf12(alpha1, alpha2);
//    surfaceScalarField snv = mixture.nHatf12();

    stOn ? Info << "Marg On" << endl : Info << "Marg Off" << endl;
    Info << "UEqn.H" << endl;
    velSource.sourceOnCheck() ? uSource = velSource.applyField()*velSourceUnits : uSource = velSourceZeroField;
    Info << "UEqn" << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
// Continuity error is applied twice?
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
//      - fvm::Sp(contErr, U)
//      + MRF.DDt(rho, U)
//      + turbulence->divDevRhoReff(rho, U)
//      - (stOn ? 1.0 : 0.0)*dsigmadT*(fvc::grad(T) - fvc::reconstruct(snv)*(fvc::reconstruct(snv) & fvc::grad(T)))*mag(fvc::grad(alphaMetalTotal))*(2*rho/ (rho1 + rho2))*alpha2
//      - (stOn ? 1.0 : 0.0)*dsigmadT*(fvc::grad(T) - fvc::reconstruct(snv)*(fvc::reconstruct(snv) & fvc::grad(T)))*amg*(mixture.rhoBouss()/(mixture.rhoBouss() + 1e-6*rho2))
      - S_Marg
      - fvm::Sp(S_d, U)
      - fvm::Sp(S_dg, U)
      - lorentz*lorentzUnits
     ==
        fvOptions(rho, U)
        + uSource
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  mixture.surfaceTensionForce()//*fvc::interpolate(topMetal)//*fvc::interpolate(stFactor)
//                  - (vsfRho ? ghf*fvc::snGrad(rho) : ghf*fvc::snGrad(rho*(scalar(1.0) - beta*boundedAlpha2*(T - Tref))))
                - (vsfRho ? (rhofT ? ghf*fvc::snGrad((boundedAlpha1*rho1 + mixture.rhoBouss() - boundedAlpha2*rhoRef + mixture.rhoExp())*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2))) 
                        : ghf*fvc::snGrad(rho*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2)))) 
                        : ghf*fvc::snGrad(rho*(scalar(1.0) - beta*boundedAlpha2*(T - Tref))*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2))))    
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    Info << "UEqn End" << endl;
        
        fvOptions.correct(U);
    }
