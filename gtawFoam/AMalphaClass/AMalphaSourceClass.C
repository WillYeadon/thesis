/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Alpha CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

    Class to create up to three heat AMalphas with a elliptic paraboloid shape
    Modified from http://dx.doi.org/10.1016/j.ijheatmasstransfer.2016.04.064

\*---------------------------------------------------------------------------*/

#include "AMalphaClass/AMalphaSourceClass.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

Foam::AMalphaSourceClass::AMalphaSourceClass(
        const dictionary& AMalphaProperties,
        const volVectorField& U,
        const surfaceScalarField& phi
):

    AMalphaPropertiesDict(AMalphaProperties),
    U_(U),
    phi_(phi),
    
    // N.B this assummes mesh starts from (0 0 0) if it starts from (10 10 10) or something you'll probs get a fpe
    AMalphaVel("AMalphaVel", dimLength, AMalphaPropertiesDict.lookupOrDefault<vector>("AMalphaVel", Zero)),
    AMalphaPos0("AMalphaPos0", dimLength, AMalphaPropertiesDict.lookupOrDefault<vector>("AMalphaPos0", Zero)),
    AMalphaOmega("AMalphaOmega", dimLength, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaOmega", 1.0)),
    AMalphaMult("AMalphaMult", dimless, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaMult", 1.0)),
    
    AMalphaStartTime("AMalphaStartTime", dimTime, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaStartTime", 0.0)),
    AMalphaPauseTime("AMalphaPauseTime", dimTime, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaPauseTime", 360.0)),
    AMalphaEndTime("AMalphaEndTime", dimTime, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaEndTime", 720.0)),
    
    AMalphaFieldCut("AMalphaFieldCut", dimless, AMalphaPropertiesDict.lookupOrDefault<scalar>("AMalphaFieldCut", 0.01)),
    AMalphaPause(AMalphaPropertiesDict.lookupOrDefault<bool>("AMalphaPause", false)),
    AMalphaOn(AMalphaPropertiesDict.lookupOrDefault<bool>("AMalphaOn", false)),
    AMalphaMod(AMalphaPropertiesDict.lookupOrDefault<bool>("AMalphaMod", false)),

    AMalphaFieldSize
    (
        IOobject
        (
            "AMalphaFieldSize",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMalphaFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    AMalphaFieldOutput
    (
        IOobject
        (
            "AMalphaFieldOutput",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMalphaFieldOutput", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
  
}

Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::geoField(const dimensionedVector AMalphaPos, const dimensionedScalar omega)
{
    const volVectorField& cellCentre = U_.mesh().C();

    scalar C_d = log(AMalphaFieldCut.value());

    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    scalar currentTime = U_.mesh().time().value();

    bool AMalphaPauseSwitch;
    (currentTime < AMalphaPauseTime.value()) ? AMalphaPauseSwitch = true : AMalphaPauseSwitch = false;

    volScalarField expX
    (
//        exp(C_d*(sqr((cellCentre.component(vector::X)-(AMalphaPos.component(vector::X)))/oneMeter)/sqr(omega/oneMeter)))
        exp(C_d*(sqr((cellCentre.component(vector::X)-(AMalphaPos.component(vector::X))
        + (AMalphaVel.component(vector::X)*(currentTime  - AMalphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    );

    volScalarField expY
    (
        exp(C_d*(sqr((cellCentre.component(vector::Y)-(AMalphaPos.component(vector::Y))
        + (AMalphaVel.component(vector::Y)*(currentTime  - AMalphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 

    volScalarField expZ
    (
        exp(C_d*(sqr((cellCentre.component(vector::Z)-(AMalphaPos.component(vector::Z))
        + (AMalphaVel.component(vector::Z)*(currentTime  - AMalphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 


    if (AMalphaPause && AMalphaPauseSwitch)      
    {
      expX = exp(C_d*(sqr((cellCentre.component(vector::X)-(AMalphaPos.component(vector::X)))/oneMeter)/sqr(omega/oneMeter)));
      expY = exp(C_d*(sqr((cellCentre.component(vector::Y)-(AMalphaPos.component(vector::Y)))/oneMeter)/sqr(omega/oneMeter)));
      expZ = exp(C_d*(sqr((cellCentre.component(vector::Z)-(AMalphaPos.component(vector::Z)))/oneMeter)/sqr(omega/oneMeter)));
    }

    volScalarField expTot
    (
      expX*expY*expZ
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "geoField",
            expTot
    )
  );  
}


Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::genField(const dimensionedVector AMalphaPos, const dimensionedScalar omega)
{
    volScalarField genField
    (
        geoField(AMalphaPos, omega)
    );

    forAll(genField, i)
    {
        (genField[i] >= AMalphaFieldCut.value()) ? genField[i] = scalar(1) : genField[i] = scalar(0);
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "genField",
            genField
        )
  );  

}

Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::normField(const volScalarField& AMalphaField)
{
    scalar fieldMax = max(AMalphaField).value();  
    
    volScalarField AMalphaFieldNorm
    (
        AMalphaField/fieldMax
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMalphaFieldNorm",
            AMalphaFieldNorm
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::binField(const volScalarField& AMalphaField)
{
    volScalarField AMalphaFieldBin(AMalphaField);

    forAll(AMalphaFieldBin, i)
    {
        if (AMalphaFieldBin[i] >= 1.0)
            AMalphaFieldBin[i] = 1.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMalphaFieldBin",
            AMalphaFieldBin
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::floorField(const volScalarField& AMalphaField)
{
    volScalarField AMalphaFieldFloor(AMalphaField);

    forAll(AMalphaFieldFloor, i)
    {
        if (AMalphaFieldFloor[i] <= 1.0)
            AMalphaFieldFloor[i] = 0.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMalphaFieldFloor",
            AMalphaFieldFloor
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::AMalphaSourceClass::applyField()
{
    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    Info << " Generating Alpha Field " << endl;
    AMalphaFieldSize = 0.0;

    if ((U_.mesh().time() > AMalphaStartTime) && (U_.mesh().time() < AMalphaEndTime)) 
        AMalphaFieldSize = binField(genField(AMalphaPos0, AMalphaOmega));

    if (U_.mesh().time() > AMalphaEndTime)
      AMalphaFieldSize = 0.0;

    scalar heatAlphaFieldMax = max(AMalphaFieldSize).value();
    scalar heatAlphaFieldMin = min(AMalphaFieldSize).value();
    Info << "heatAlphaField max: " << heatAlphaFieldMax 
    << " min: " << heatAlphaFieldMin << endl;

    if (!AMalphaOn)
        AMalphaMult = 0.0;

    volScalarField AMalpha0 = normField(geoField(AMalphaPos0, AMalphaOmega));    

    scalar AMalpha0Max = max(AMalpha0).value();
    scalar AMalpha0Min = min(AMalpha0).value();
    Info <<  "AMalpha0 max: " << AMalpha0Max << " min: " << AMalpha0Min << endl;

    AMalphaFieldOutput = AMalphaFieldSize*AMalphaMult;

    if (AMalphaMod)
        AMalphaFieldOutput = AMalphaFieldSize*AMalpha0*AMalphaMult;

    return AMalphaFieldOutput;
}