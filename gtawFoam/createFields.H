#include "createRDeltaT.H"
#include "createUnitCorr.H"

Info << "Reading Transport Properties\n" << endl;
IOdictionary transportPropertiesDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Flags
//scalar phaseChange(readScalar(transportPropertiesDict.lookup("phaseChange"))); 
scalar alphaInterfaceCut(transportPropertiesDict.lookupOrDefault<scalar>("alphaInterfaceCut", 1e-6));
bool phaseChange(transportPropertiesDict.lookupOrDefault<bool>("phaseChange", true));
bool tubeCapture(transportPropertiesDict.lookupOrDefault<bool>("tubeCapture", false));
bool airTop(transportPropertiesDict.lookupOrDefault<bool>("airTop", true));
bool tubeCaptured = false;
bool vsfRho(transportPropertiesDict.lookupOrDefault<bool>("vsfRho", false));
bool rhofT(transportPropertiesDict.lookupOrDefault<bool>("rhofT", false));
bool LfSmoothing(transportPropertiesDict.lookupOrDefault<bool>("LfSmoothing", false));
bool LfLimit(transportPropertiesDict.lookupOrDefault<bool>("LfLimit", false));
bool alphaGeo(transportPropertiesDict.lookupOrDefault<bool>("alphaGeo", true));
bool amGeo(transportPropertiesDict.lookupOrDefault<bool>("amGeo", false));
bool outerTubeKill(transportPropertiesDict.lookupOrDefault<bool>("outerTubeKill", false));
bool probesOn(transportPropertiesDict.lookupOrDefault<bool>("probesOn", false));
scalar charLength3(transportPropertiesDict.lookupOrDefault<scalar>("charLength3", 1e-6));


// Compilation Flags
bool AM = false;
(void)AM;
if (AM){for (int i = 0; i < 10; ++i){Info << "AM ON\tAM ON\tAM ON\tAM ON\tAM ON\t" << endl;}}

// Warning removal as variable used for comparisons
(void)alphaInterfaceCut;
(void)phaseChange;
(void)tubeCapture;
(void)airTop;
(void)tubeCaptured;
(void)vsfRho;
(void)rhofT;
(void)LfSmoothing;
(void)LfLimit;
(void)amGeo;
(void)outerTubeKill;
(void)alphaGeo;
(void)probesOn;
(void)charLength3;

dimensionedScalar dmdtRelaxAlpha
(
    "dmdtRelaxAlpha",
    dimensionSet(0,0,0,0,0,0,0),
    transportPropertiesDict.lookup("dmdtRelaxAlpha")
);

dimensionedScalar dmdtRelaxPressure
(
    "dmdtRelaxPressure",
//    dimensionSet(0,0,0,0,0,0,0),
    dimensionSet(-1,5,-3,1,0,0,0),
    transportPropertiesDict.lookup("dmdtRelaxPressure")
);


dimensionedScalar meltingT
(
    "meltingT",
    dimensionSet(0,0,0,1,0,0,0),
    transportPropertiesDict.lookup("meltingT")
);

dimensionedScalar sFillerUnitsCorr
(
    "sFillerUnitsCorr",
    dimensionSet(1,-1,-2,-1,0,0,0),
    10.0
);

dimensionedScalar cp
(
    "cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    1.0
);

dimensionedScalar Lf
(
    "Lf",
    dimensionSet(0,2,-2,0,0,0,0),
    transportPropertiesDict.lookup("Lf")
);

dimensionedScalar rhoRef
(
    "rhoRef",
    dimensionSet(1,-3,0,0,0,0,0),
    transportPropertiesDict.lookupOrDefault("rhoRef", 1000.0)
);

dimensionedScalar beta
(
    "beta",
    dimensionSet(0,0,0,-1,0,0,0),
    transportPropertiesDict.lookup("beta")
);

dimensionedScalar beta_2
(
    "beta_2",
    dimensionSet(0,0,0,-1,0,0,0),
    transportPropertiesDict.lookupOrDefault("beta_2", 0.0)
);

dimensionedScalar sourcedeltaT
(
    "sourcedeltaT",
    dimensionSet(0,0,1,0,0,0,0),
    1.0
);

dimensionedScalar sourceTime
(
    "sourceTime",
    dimensionSet(0,0,1,0,0,0,0),
    1.0
);

dimensionedScalar AMtriggerTime
(
    "sourceTime",
    dimensionSet(0,0,1,0,0,0,0),
    transportPropertiesDict.lookupOrDefault("AMtriggerTime", 0.0)
);
dimensionedScalar Tref
(
    "Tref",
    dimensionSet(0,0,0,1,0,0,0),
    transportPropertiesDict.lookup("Tref")
);

dimensionedScalar Tref_2
(
    "Tref_2",
    dimensionSet(0,0,0,1,0,0,0),
    transportPropertiesDict.lookupOrDefault("Tref_2", 0.0)
);


// Used for H comparison
dimensionedScalar T_0
(
    "T_0",
    dimensionSet(0,0,0,1,0,0,0),
    transportPropertiesDict.lookup("T_0")
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Making mixture\n" << endl;
immiscibleIncompressibleThreePhaseMixture mixture(U, phi);//, Tref, beta);
volScalarField& T(mixture.T());
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

const dimensionedScalar& cp1 = mixture.cp1();
const dimensionedScalar& cp2 = mixture.cp2();
const dimensionedScalar& cp3 = mixture.cp3();

const dimensionedScalar& k1 = mixture.k1();
const dimensionedScalar& k2 = mixture.k2();
const dimensionedScalar& k3 = mixture.k3();

mixture.callUpdateCpMembers();
Info<< "Mixture made\n" << endl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Mass flux
/*
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0)
);
*/
/*
Info<< "Reading field H\n" << endl;
volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
//        IOobject::AUTO_WRITE
    ),
    (T - T_0)*(alpha1*cp1*rho1 + alpha2*cp2*rho2 + alpha3*cp3*rho3)/rho,
    T.boundaryField().types()
//    mesh
);
H.oldTime();
*/
// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
//#include "ghMod.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
#include "createMRF.H"
#include "createFvOptions.H"

dimensionedScalar multiplierBuoy
(
    "multiplierBuoy",
    dimensionSet(0,0,0,0,0,0,0),
    transportPropertiesDict.lookup("multiplierBuoy")//1.6e3//1e6
);

dimensionedScalar multiplierBuoyG
(
    "multiplierBuoyG",
    dimensionSet(0,0,0,0,0,0,0),
    transportPropertiesDict.lookup("multiplierBuoyG")//1.6e3//1e6
);

dimensionedScalar dsigmadT
(
    "dsigmadT",
    dimensionSet(1,0,-2,-1,0,0,0),
    transportPropertiesDict.lookup("dsigmadT")
);


dimensionedScalar epsilonBuoy
(
    "epsilonBuoy",
    dimensionSet(0,0,0,0,0,0,0),
    1e-3
);

Info<< "Solid Mom kill\n" <<endl;
volScalarField S_d
(
    IOobject
    (
        "S_d",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -multiplierBuoy*(sqr(alpha3)/(pow3(1.0 - alpha3) + 1e-3))*unitsBuoy
);

// Dev, just for gas phase
volScalarField S_dg
(
    IOobject
    (
        "S_dg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -multiplierBuoyG*(sqr(1.0 - alpha1)/(pow3(alpha1) + 1e-3))*unitsBuoy
);

volVectorField S_Marg
(
    IOobject
    (
        "S_Marg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),   
    mesh,
    dimensionedVector("S_Marg",dimensionSet(1,-2,-2,0,0,0,0), Zero)
);


Info<< "Reading / calculating rho*cp\n" <<endl;
volScalarField rhoCp
(
    IOobject
    (
    "rhoCp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    (alpha1*rho1 + alpha2*rho2 + alpha3*rho3)*cp//,
);
rhoCp.oldTime();

dimensionedScalar rhoCpSmall
(
    "rhoCpSmall",
    rhoCp.dimensions(),
    1e-3
);


volScalarField rhoLf
(
    IOobject
    (
    "rhoLf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
//    ((alpha1*rho1/cp1) + (alpha2*rho2/cp2) + (alpha3*rho3/cp3))*Lf
    (alpha1*rho1 + alpha2*rho2 + alpha3*rho3)*Lf//,
//    (alpha1*rho1 + alpha2*rho2Bouss + alpha3*rho3)*Lf//,
//    alpha2.boundaryField().types()
);
rhoLf.oldTime();

dimensionedScalar rhoCp_units
(
    "rhoCp_units",
    rhoCp.dimensions(),
    1.0
);

volScalarField rhoCp_units_object
(
    IOobject
    (
    "rhoCp_units_object",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    mesh,
    rhoCp_units
);

Info<< "Reading / calculating rhoPhi*cp\n" <<endl;
surfaceScalarField rhoPhiCp
(
    IOobject
    (
    "rhoPhiCp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

volScalarField changeAlpha1
(
    IOobject
    (
        "changeAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("changeAlpha1",dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField interfaceAlpha2
(
    IOobject
    (
        "interfaceAlpha2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("interfaceAlpha2",dimensionSet(0,0,0,0,0,0,0), 0.0)
);

// Phase Change
volScalarField dmdt
(
    IOobject
    (
        "dmdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*(T - meltingT)*cp/Lf
);

volScalarField dmdtSp
(
    IOobject
    (
        "dmdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*(T - meltingT)*cp/Lf
);

volScalarField negDmdt
(
    IOobject
    (
        "negDmdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("negDmdt",dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField boundedAlpha1
(
    IOobject
    (
        "boundedAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    min(max(alpha1, scalar(0)), scalar(1))
);

volScalarField boundedAlpha2
(
    IOobject
    (
        "boundedAlpha2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    min(max(alpha2, scalar(0)), scalar(1))
);

volScalarField boundedAlpha3
(
    IOobject
    (
        "boundedAlpha3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    min(max(alpha3, scalar(0)), scalar(1))
);

volScalarField bodyAlpha2
(
    IOobject
    (
        "bodyAlpha2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha2
);

Info<< "Creating field k_field\n" << endl;
volScalarField k_field
(
    IOobject
    (
        "k_field",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    k1*boundedAlpha1 + k2*boundedAlpha2 + k3*boundedAlpha3
);

volScalarField Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sp",dimVelocity/dimLength, 0)
);

volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Su",dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField boundedAlpha
(
    IOobject
    (
        "boundedAlpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1
);

volScalarField noInterfaceAlpha2
(
    IOobject
    (
        "noInterfaceAlpha2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha2
);



volScalarField boundedAlphaChange
(
    IOobject
    (
        "boundedAlphaChange",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("boundedAlphaChange",dimensionSet(0,0,0,0,0,0,0), 0.0)
);

dimensionedScalar k_weld
(
    "k_weld",
    dimensionSet(1,1,-3,-1,0,0,0),
    transportPropertiesDict.lookup("k_weld")
);
volScalarField alpha1Accumulation
(
    IOobject
    (
        "alpha1Accumulation",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1
);

volScalarField alpha1AccumulationNeg
(
    IOobject
    (
        "alpha1AccumulationNeg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1.0 - alpha1Accumulation
);

volScalarField correctedAlpha1Accumulation
(
    IOobject
    (
        "correctedAlpha1Accumulation",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1Accumulation - alpha2 - alpha3
);

volScalarField correctedAlpha1AccumulationNeg
(
    IOobject
    (
        "correctedAlpha1AccumulationNeg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1.0 - correctedAlpha1Accumulation
);

volScalarField alpha2Accumulation
(
    IOobject
    (
        "alpha2Accumulation",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    alphaLiquid
    alpha2
);

volScalarField alpha3Accumulation
(
    IOobject
    (
        "alpha3Accumulation",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    alphaLiquid
    alpha3
);

volScalarField alphaMetalTotal
(
    IOobject
    (
        "alphaMetalTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha2 + alpha3)
);

volScalarField S_ld
(
    IOobject
    (
        "S_ld",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    -(rho*Lf/(oneSecond))*(alphaLiquid/oneKelvin)
    -(rho*Lf/(oneSecond))*(alpha2/oneKelvin)
);

volScalarField S_lod
(
    IOobject
    (
        "S_lod",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    -S_ld*(T.oldTime()) + (rho*Lf/oneSecond)*(alpha1.oldTime() - alpha1)
    -S_ld*(T.oldTime()) - rho*Lf*fvc::ddt(boundedAlphaChange)
);

volScalarField factor
(
    IOobject
    (
        "factor",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (scalar(3)*rho)/(rho1 + rho2 + rho3)
);



volScalarField channel
(
    IOobject
    (
        "channel",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1
);

volScalarField rhoLfCp
(
    IOobject
    (
    "rhoLfCp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    ((boundedAlpha1*rho1/cp1) + (boundedAlpha2*rho2/cp2) + (boundedAlpha3*rho3/cp3))*Lf
);
rhoLfCp.oldTime();


volScalarField S_latent
(
    IOobject
    (
        "S_latent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    rhoLf*fvc::ddt(alpha2)
//    rhoLfCp*fvc::ddt(alpha2)
);

surfaceScalarField kEff
(
    IOobject
    (
        "kEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate( mixture.k() )
//    + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
);

//#include "createIncompressibleRadiationModel.H"


#include "createEMFields.H"


volScalarField alpha2WeightedAv
(
    IOobject
    (
        "alpha2WeightedAv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    alpha2*0.0
);

dimensionedScalar alpha2Capture
(
    "alpha2Capture",
    dimensionSet(0,0,0,0,0,0,0),
    transportPropertiesDict.lookupOrDefault("alpha2Capture", 1.0)//1.6e3//1e6
);

dimensionedScalar nu1_shadow
(
    "nu1_shadow",
    dimensionSet( 0, 2, -1, 0, 0, 0, 0),
    transportPropertiesDict.lookupOrDefault("nu1_shadow", 1e-06)
);

dimensionedScalar nu2_shadow
(
    "nu2_shadow",
    dimensionSet( 0, 2, -1, 0, 0, 0, 0),
    transportPropertiesDict.lookupOrDefault("nu2_shadow", 1e-06)
);

dimensionedScalar nu3_shadow
(
    "nu3_shadow",
    dimensionSet( 0, 2, -1, 0, 0, 0, 0),
    transportPropertiesDict.lookupOrDefault("nu3_shadow", 1e-06)
);

volScalarField nuField
(
    IOobject
    (
        "PrandtlNumber",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*nu1_shadow + alpha2*nu2_shadow + alpha3 * nu3_shadow
);

dimensionedScalar catchcatch
(
    "catchcatch",
    dimensionSet( 1, 1, -3, -1, 0, 0, 0),
    1e-06
);


Info<< "Calculating PrandtlNumber\n" << endl;
/*
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
    ((alpha1*rho1*cp1 + alpha2*rho2*cp2 + alpha3*rho3*cp3)/(alpha1*k1 + alpha2*k2 + alpha3*k3 + catchcatch))/nuField
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
    ell*ell*ell*mag(g.value())*beta*(T - Tref)/Foam::sqr(nuField)
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
dimensionedScalar rhoNoUpdateScalar
(
    "rhoNoUpdateScalar",
    dimensionSet(1,-3,0,0,0,0,0),
    1.0
);
volScalarField rhoNoUpdate
(
    IOobject
    (
        "rhoNoUpdate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoNoUpdateScalar*(alpha1 + alpha2 + alpha3)
);

volScalarField cpLiquidField
(
    IOobject
    (
        "cpLiquidField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    cp2*min(max((alpha1 + (1.0 - alpha1)), scalar(1)), scalar(1))
//    (cp_T_control) ? cp2_shadow.weightedAverage(mesh.V()).value()*min(max((alpha1 + (1.0 - alpha1)), scalar(1)), scalar(1))
//    : cp2*min(max((alpha1 + (1.0 - alpha1)), scalar(1)), scalar(1))
);

dimensionedScalar rhoOne
(
    "rho0",
    dimensionSet(1,-3,0,0,0,0,0),
    1.0
);



Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// 5283 Below works but just one phase, needs to be after rest so T updates
/*
#include "compressibleTransportModel.H"
Info << "\t\t SPAC \t\t" << endl;
autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
#include "createRadiationModel.H"
*/
// 5283 End 


// 24601 FOR RADIATION
/*
#include "createRadiationGIF.H"
*/
//\24601

volScalarField alphaDyn
(
    IOobject
    (
        "alphaDyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha2
);

volScalarField stFactor
(
    IOobject
    (
        "stFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    2*rho/(rho2 + rho3)
);

volScalarField sFillerAcc
(
    IOobject
    (
        "sFillerAcc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.0*(alpha1 + alpha2 + alpha3)/oneSecond
);

volScalarField baseMetal
(
    IOobject
    (
        "baseMetal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (alpha2 + alpha3)
);

volScalarField AMbaseMetal
(
    IOobject
    (
        "baseMetal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (alpha2 + alpha3)
);

volScalarField coreRegion
(
    IOobject
    (
        "coreRegion",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1 + alpha2 + alpha3
);


volScalarField alpha1OldTime
(
    IOobject
    (
        "alpha1OldTime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1.oldTime()
);

volScalarField alpha3OldTime
(
    IOobject
    (
        "alpha3OldTime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha3.oldTime()
);


volScalarField collapseCapture
(
    IOobject
    (
        "collapseCapture",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    coreRegion*(alpha1OldTime - 1e5*(alpha1OldTime*alpha3OldTime))
);

volScalarField burstCapture
(
    IOobject
    (
        "burstCapture",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    coreRegion*(alpha3OldTime - 1e5*(alpha1OldTime*alpha3OldTime))
);

volScalarField interfaceCapture
(
    IOobject
    (
        "interfaceCapture",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    coreRegion*(alpha1 + alpha2 + alpha3)
);

volScalarField fullWeld
(
    IOobject
    (
        "fullWeld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha2
);


if (airTop)
{
    burstCapture = coreRegion*(alpha1OldTime - 1e5*(alpha1OldTime*alpha3OldTime));
    collapseCapture = coreRegion*(alpha3OldTime - 1e5*(alpha1OldTime*alpha3OldTime));
}



/*
volScalarField topMetal
(
    IOobject
    (
        "topMetal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (alpha1 + alpha2 + alpha3)
);
*/
bool stOn = false;
if (mixture.sigma12Access().value() != 0 && dsigmadT.value() != 0)
    stOn = true;

(void)stOn;

int globalPenn{0};
int globalCollapse{0};
int globalBurst{0};
(void)globalPenn;
(void)globalCollapse;
(void)globalBurst;