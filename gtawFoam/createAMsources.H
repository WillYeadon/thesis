IOdictionary AMsourcePropertiesDict
(
    IOobject
    (
        "AMsourceProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

AMheatSourceClass AMarcSource(AMsourcePropertiesDict, alphaMetalTotal, U, phi);

volScalarField AMqSource
(
    IOobject
    (
        "AMqSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    AMarcSource.applyField(alphaMetalTotal)*qSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);

IOdictionary AMfillerPropertiesDict
(
    IOobject
    (
        "AMfillerProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

AMalphaSourceClass AMalphaSource(AMfillerPropertiesDict, U, phi);

volScalarField AMsFiller
(
    IOobject
    (
        "AMsFiller",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    AMalphaSource.applyField()*alphaSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);
