#include "MULES.H"

//{
//    rho = mixture.rho();
//    #include "alphaInterfaceCapture.H"

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    PtrList<surfaceScalarField> alphaPhiCorrs(3);

    phic = min(mixture.cAlpha()*phic, max(phic));

//    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
//    {
        #include "alpha1PhiGen.H"
        #include "alpha2PhiGen.H"
        #include "alpha3PhiGen.H"

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi = dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0);
    rhoPhiCp *= 0.0;

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("sumAlpha", dimless, 0)
    );

    sFiller = alphaSource.applyField()*alphaSourceUnits;
    sFiller *= (1.0 - baseMetal);

        surfaceScalarField& alpha1Phi = alphaPhiCorrs[0];
        alpha1Phi += upwind<scalar>(mesh, phi).flux(alpha1);

        if (MULESCorr)
        {
            solve
            (
                fvm::ddt(alpha1) + fvc::div(alpha1Phi)
                ==
//                - Su*alpha1
                - fvOptions(alpha1)
                - sFiller
            );
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                alpha1Phi
            );
        }       

//        rhoPhi += alpha1Phi*(rho1);
        rhoPhi += alpha1Phi*rho1;
        rhoPhiCp += alpha1Phi*rho1*cp1;
        sumAlpha += alpha1;

        surfaceScalarField& alpha2Phi = alphaPhiCorrs[1];
        alpha2Phi += upwind<scalar>(mesh, phi).flux(alpha2);

        #include "alphaSuSp.H"

        if (MULESCorr)
        {
            solve
            (
                fvm::ddt(alpha2) + fvc::div(alpha2Phi)
                ==
                Su
                + sFiller
//                (Su + fvm::Sp(Sp + divU, alpha2))
//                + fvOptions(alpha2)
            );
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha2,
                phi,
                alpha2Phi,
                divU,//Sp,
                Su,
                oneField(),
                zeroField()
            );
        }
/*
        volScalarField rho2Vsf = mixture.rhoBouss();
        volScalarField rho2CpVsf = mixture.rhoBouss()*mixture.cp2TDepVSF();

        forAll(rhoPhi, i)
        {
            if(rho2Vsf[i] > 1e-6)
                rhoPhi[i] += alpha2Phi[i]*rho2Vsf[i];
            if(rho2CpVsf[i] > 1e-6)
                rhoPhiCp[i] += alpha2Phi[i]*rho2CpVsf[i];
        }
*/
//        rhoPhi += alpha2Phi*fvc::interpolate(mixture.rhoBouss());
//        rhoPhiCp += alpha2Phi*fvc::interpolate(mixture.rhoBouss())*fvc::interpolate(mixture.cp2TDepVSF());
        rhoPhi += alpha2Phi*rho2;
        rhoPhiCp += alpha2Phi*rho2*cp2;
        sumAlpha += alpha2;
        
        surfaceScalarField& alpha3Phi = alphaPhiCorrs[2];
        alpha3Phi += upwind<scalar>(mesh, phi).flux(alpha3);

        if (MULESCorr)
        {
            solve
            (
                fvm::ddt(alpha3) + fvc::div(alpha3Phi)
                ==
                dmdt/oneSecond
//                (dmdt/oneSecond + fvm::Sp(dmdtSp/oneSecond + divU, alpha3))
            );
        }
        else
        {
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha3,
            phi,
            alpha3Phi,
            divU,//dmdtSp,
            dmdt,
            oneField(),
            zeroField()
        );
        }
/*
        volScalarField rho3Vsf = mixture.rhoExp();
        volScalarField rho3CpVsf = mixture.rhoExp()*mixture.cp3TDepVSF();

        forAll(rhoPhi, i)
        {
            if(rho3Vsf[i] > 1e-6)
                rhoPhi[i] += alpha3Phi[i]*rho3Vsf[i];
            if(rho3CpVsf[i] > 1e-6)
                rhoPhiCp[i] += alpha3Phi[i]*rho3CpVsf[i];
        }
*/
        rhoPhi += alpha3Phi*(rho3);
        rhoPhiCp += alpha3Phi*(rho3)*cp3;
//        rhoPhi += alpha3Phi*fvc::interpolate(mixture.rhoExp());
        sumAlpha += alpha3;
/*
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha2,
            alpha2Phi
        );            

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha3,
            alpha3Phi
        );
*/

/*
*/
        volScalarField sumCorr(1.0 - sumAlpha);
        alpha1 += alpha1*sumCorr;
        alpha2 += alpha2*sumCorr;
        alpha3 += alpha3*sumCorr;

        alpha1 = max(min(alpha1, scalar(1)), scalar(0));
        alpha2 = max(min(alpha2, scalar(1)), scalar(0));
        alpha3 = max(min(alpha3, scalar(1)), scalar(0));

//        #include "alphaClean.H"

        forAll(alpha1, i)
        {
            if (alpha1[i] < 1e-3)
                alpha1[i] = 0.0;
            if (alpha2[i] < 1e-3)
                alpha2[i] = 0.0;
            if (alpha3[i] < 1e-3)
                alpha3[i] = 0.0;
        }

//        alpha2 = 1.0 - alpha1 - alpha3;

        Info << "2" << endl;
//        alphaMetalTotal = alpha2 + alpha3;
//        alphaMetalTotal = max(min(alphaMetalTotal, scalar(1)), scalar(0));
/*
        rhoPhi = 
            alpha1Phi*(rho1)
          + alpha2Phi*(rho2)
          + alpha3Phi*(rho3);
    rhoPhiCp = 
            alpha1Phi*(rho1)*cp1
          + alpha2Phi*(rho2)*cp2
          + alpha3Phi*(rho3)*cp3;
*/

//}