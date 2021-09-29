Info << "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
	IOobject
	(
		"physicalProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

bool threePhaseEM(physicalProperties.lookupOrDefault<bool>("threePhaseEM", false));
// Warning removal as variable used for comparisons
(void)threePhaseEM;

pisoControl piso(mesh);
pisoControl bpiso(mesh, "BPISO");

//const dictionary& BControls = mesh.solverDict(B.name());

//label nBCorr(readLabel(BControls.lookup("nBCorr")));

// k is thermal conductivity
dimensionedScalar k_EM
(
    "k_EM",
    dimensionSet(1,1,-3,-1,0,0,0),
    physicalProperties
);

// Magnetic Permeability
/*
dimensionedScalar muMag
(
	"muMag",
	dimensionSet( 1, 1, -2, 0, 0, -2, 0),
	physicalProperties // reads from dictionary physicalProperties
);
*/
dimensionedScalar d2dm2
(
	"d2dm2",
	dimensionSet( 0, -2, 0, 0, 0, 0, 0),
	1.0
);

// Permeativity
dimensionedScalar epsilon
(
	"epsilon",
	dimensionSet( -1, -3, 4, 0, 0, 2, 0),
	physicalProperties
);

// Electrical conductivity Ohms^-1
dimensionedScalar sigma_l
(
	"sigma_l",
	dimensionSet( -1, -3, 3, 0, 0, 2, 0),
	physicalProperties // reads from dictionary physicalProperties
);

//Magnetic field
dimensionedVector B0
(
	"B0",
	dimensionSet(1, 0, -2, 0, 0, -1, 0),
	physicalProperties
);

Info<< "Reading field sigma\n" << endl;
volScalarField muMagR
(
	IOobject
	(
	"muMagR",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("muMagR", dimensionSet( 0, 0, 0, 0, 0, 0, 0), physicalProperties)
);

Info << "Reading field ElPot\n" << endl; // this is read from the mesh
volScalarField ElPot
(
	IOobject
	(
		"ElPot",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
//	dimensionedScalar("ElPot",dimensionSet(1,2,-3,0,0,-1,0), 0.0)
);
// N.B can get the units from 0/ElPot etc...

label ElPotRefCell = 0;
scalar ElPotRefValue = 0.0;
setRefCell(ElPot, mesh.solutionDict().subDict("ElPot"), ElPotRefCell, ElPotRefValue);
mesh.setFluxRequired(ElPot.name());


Info << "Reading field A\n" << endl;
volVectorField A
(
	IOobject
	(
		"A",
		runTime.timeName(),
		mesh,
//		IOobject::MUST_READ,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh
//	dimensionedVector("A",dimensionSet(1,1,-1,0,0,-1,0), Zero)
);

Info << "Calculating Magnetic Field B\n" << endl;
volVectorField B
(
	IOobject
	(
		"B",
		runTime.timeName(),
		mesh,
//		IOobject::NO_READ, // A-V interpretation
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
//	fvc::curl(A) // A-V interpretation
	mesh
);

//surfaceScalarField psiub = fvc::interpolate(U ^ B0) & mesh.Sf();
Info<< "Calculating Field UxB\n" << endl;
volVectorField UxB
(
	IOobject
	(
		"UxB",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	U ^ B0
);
Info<< "Calculating Field UxBS\n" << endl;
volScalarField UxBS
(
	IOobject
	(
		"UxBS",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	fvc::div(UxB)
);
Info<< "Calculating Field US\n" << endl;
volScalarField US
(
	IOobject
	(
		"US",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	fvc::div(U)
);
Info<< "Calculating Field Psi UxB\n" << endl;
surfaceScalarField psiUB
(
    IOobject
    (
    "psiUB",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
	fvc::interpolate(UxB) & mesh.Sf()
//    linearInterpolate(UxB) & mesh.Sf()
);

Info << "Calc Jn\n" << endl;
surfaceScalarField Jn
(
	IOobject
	(
		"Jn",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	-(fvc::snGrad(ElPot) * mesh.magSf()) + psiUB
);
Info << "Calc Jnv\n" << endl;
surfaceVectorField Jnv
(
	IOobject
	(
		"Jnv",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	Jn * mesh.Cf()
);

Info<< "Calculating Field J0\n" << endl;
volVectorField J0
(
	IOobject
	(
		"J0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	fvc::surfaceIntegrate(Jnv) - (fvc::surfaceIntegrate(Jn) * mesh.C())
);
/*
*/
/*
Info<< "Calculating Field JxB\n" << endl;
volVectorField JxB
(
	IOobject
	(
		"JxB",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	J0 ^ B
);
*/
dimensionedScalar nuEM
(
	"nuEM",
	dimensionSet( 0, 2, -1, 0, 0, 0, 0),
	physicalProperties // reads from dictionary physicalProperties
);

dimensionedScalar ell
(
	"ell",
	dimensionSet( 0, 1, 0, 0, 0, 0, 0),
	1.0
);

Info<< "Calculating Probes\n" << endl;
/*
Info<< "Calculating HartmanNumber\n" << endl;
volVectorField HartmanNumber
(
    IOobject
    (
        "HartmanNumber",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    B*ell*Foam::sqrt((sigma_l)/(rho*nuEM))
);

Info<< "Calculating PrandtlNumber\n" << endl;
volScalarField PrandtlNumber
(
    IOobject
    (
        "PrandtlNumber",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rhoCp/k_EM)*nuEM
);

Info<< "Calculating GrashofNumber\n" << endl;
volScalarField GrashofNumber
(
    IOobject
    (
        "GrashofNumber",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ell*ell*ell*mag(g.value())*beta*(T - Tref)/Foam::sqr(nuEM)
);

Info<< "Calculating RayleighNumber\n" << endl;
volScalarField RayleighNumber
(
    IOobject
    (
        "RayleighNumber",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    GrashofNumber*PrandtlNumber
);
*/
Info << "Calc lorentz\n" << endl;
volVectorField lorentz
(
	IOobject
	(
		"lorentz",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, // it is calculated
		IOobject::AUTO_WRITE
	),
	sigma_l*((fvc::grad(ElPot) ^ B0) + ((U ^ B0) ^ B0))
);

volScalarField alphaEM
(
    IOobject
    (
        "alphaEM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha1 + alpha2 + alpha3)
);


Info<< "Making Source\n" << endl;
emSourceClass emSourceDec(physicalProperties, U, phi);

dimensionedScalar emSourceUnits
(
    "emSourceUnits",
    dimensionSet(1,-1,-2,0,0,-1,0),
    1.0
);
dimensionedScalar lorentzUnits
(
    "lorentzUnits",
    dimensionSet(0,0,0,0,0,0,0),
    1.0
);
/*
*/
volVectorField emSource
(
    IOobject
    (
        "emSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    emSourceDec.applyField()*emSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);
