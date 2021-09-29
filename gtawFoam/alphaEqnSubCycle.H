(void)alphaApplyPrevCorr;
(void)icAlpha;
(void)scAlpha;

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    tmp<volScalarField> trSubDeltaT;
    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }
/*  
    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
*/
    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    if (AM)
    {
        #include "createAMsources.H"
        #include "AMalphaEqn.H"
    }
    else
    {    
        #include "alphaEqn.H"
    }
}

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
    )()
);

Info << "alphaEqn Complete" << endl;

        