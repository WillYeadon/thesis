{

emSource = emSourceDec.applyField()*emSourceUnits;

if (threePhaseEM)
{
  forAll(alphaEM, i)
  {
    ((alpha1[i] > 0.25) && (alpha1[i] < 0.75) && (alphaMetalTotal[i] > 0.25) && (alphaMetalTotal[i] < 0.75))
    ? alphaEM = 1.0 : alphaEM = 0.0; 
  }

  emSource *= alphaEM;
}

UxB = U ^ emSource;
US = fvc::div(U);
UxBS = fvc::div(UxB);
psiUB = linearInterpolate(UxB) & mesh.Sf(); 

Info<< " Solving Elpot " << endl;
while (bpiso.correct())
{
  fvScalarMatrix ElPotEqn
  {
    fvm::laplacian(ElPot) 
    ==
    fvc::div(UxB)   
  };

  ElPotEqn.setReference(ElPotRefCell, ElPotRefValue);
  ElPotEqn.solve();
};

Jn = -(fvc::snGrad(ElPot) * mesh.magSf()) + psiUB;
Jnv = Jn * mesh.Cf();
J0 = fvc::surfaceIntegrate(Jnv) - (fvc::surfaceIntegrate(Jn) * mesh.C());
J0.correctBoundaryConditions();

lorentz = (sigma_l * (J0 ^ emSource));
/*
if (threePhaseEM)
{
  forAll(alphaEM, i)
  {
    ((alpha1[i] > 0.25) && (alpha1[i] < 0.75) && (alphaMetalTotal[i] > 0.25) && (alphaMetalTotal[i] < 0.75))
    ? alphaEM = 1.0 : alphaEM = 0.0; 
  }

  lorentz = alphaEM * (sigma_l * (J0 ^ emSource));
}
else
{
  lorentz = alphaEM * (sigma_l * (J0 ^ emSource));
}
*/


//HartmanNumber = B*ell*Foam::sqrt((sigma_l)/(rho*nuEM));
//PrandtlNumber = (rhoCp/k_EM)*nuEM;
//GrashofNumber = ell*ell*ell*mag(g.value())*beta*(T - Tref)/Foam::sqr(nuEM);
//RayleighNumber = GrashofNumber*PrandtlNumber;

/*
if(runTime.outputTime())
{
  ElPot.write();
  B.write();
  J0.write();
  HartmanNumber.write(); 
  PrandtlNumber.write(); 
  GrashofNumber.write(); 
  RayleighNumber.write(); 
  UxB.write();
  US.write();
  UxBS.write();
  lorentz.write();
}
*/
}