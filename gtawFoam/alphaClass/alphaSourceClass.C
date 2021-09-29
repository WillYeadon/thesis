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

	Class to create up to three heat alphas with a elliptic paraboloid shape
	Modified from http://dx.doi.org/10.1016/j.ijheatmasstransfer.2016.04.064

\*---------------------------------------------------------------------------*/

#include "alphaClass/alphaSourceClass.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

Foam::alphaSourceClass::alphaSourceClass(
        const dictionary& alphaProperties,
        const volVectorField& U,
        const surfaceScalarField& phi
):

    alphaPropertiesDict(alphaProperties),
    U_(U),
    phi_(phi),
    
    // N.B this assummes mesh starts from (0 0 0) if it starts from (10 10 10) or something you'll probs get a fpe
    alphaVel("alphaVel", dimLength, alphaPropertiesDict.lookupOrDefault<vector>("alphaVel", Zero)),
    alphaPos0("alphaPos0", dimLength, alphaPropertiesDict.lookupOrDefault<vector>("alphaPos0", Zero)),
    alphaOmega("alphaOmega", dimLength, alphaPropertiesDict.lookupOrDefault<scalar>("alphaOmega", 1.0)),
    alphaMult("alphaMult", dimless, alphaPropertiesDict.lookupOrDefault<scalar>("alphaMult", 1.0)),
    
    alphaStartTime("alphaStartTime", dimTime, alphaPropertiesDict.lookupOrDefault<scalar>("alphaStartTime", 0.0)),
    alphaPauseTime("alphaPauseTime", dimTime, alphaPropertiesDict.lookupOrDefault<scalar>("alphaPauseTime", 360.0)),
    alphaEndTime("alphaEndTime", dimTime, alphaPropertiesDict.lookupOrDefault<scalar>("alphaEndTime", 720.0)),
    
    alphaFieldCut("alphaFieldCut", dimless, alphaPropertiesDict.lookupOrDefault<scalar>("alphaFieldCut", 0.01)),
    alphaPause(alphaPropertiesDict.lookupOrDefault<bool>("alphaPause", false)),
    alphaOn(alphaPropertiesDict.lookupOrDefault<bool>("alphaOn", false)),
    alphaMod(alphaPropertiesDict.lookupOrDefault<bool>("alphaMod", false)),

    alphaFieldSize
    (
        IOobject
        (
            "alphaFieldSize",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("alphaFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaFieldOutput
    (
        IOobject
        (
            "alphaFieldOutput",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("alphaFieldOutput", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
  
}

Foam::tmp<Foam::volScalarField>
Foam::alphaSourceClass::geoField(const dimensionedVector alphaPos, const dimensionedScalar omega)
{
    const volVectorField& cellCentre = U_.mesh().C();

    scalar C_d = log(alphaFieldCut.value());

    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    scalar currentTime = U_.mesh().time().value();

    bool alphaPauseSwitch;
    (currentTime < alphaPauseTime.value()) ? alphaPauseSwitch = true : alphaPauseSwitch = false;

    volScalarField expX
    (
        exp(C_d*(sqr((cellCentre.component(vector::X)-(alphaPos.component(vector::X))
        + (alphaVel.component(vector::X)*(currentTime  - alphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    );

    volScalarField expY
    (
        exp(C_d*(sqr((cellCentre.component(vector::Y)-(alphaPos.component(vector::Y))
        + (alphaVel.component(vector::Y)*(currentTime  - alphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 

    volScalarField expZ
    (
        exp(C_d*(sqr((cellCentre.component(vector::Z)-(alphaPos.component(vector::Z))
        + (alphaVel.component(vector::Z)*(currentTime  - alphaPauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 


    if (alphaPause && alphaPauseSwitch)      
    {
      expX = exp(C_d*(sqr((cellCentre.component(vector::X)-(alphaPos.component(vector::X)))/oneMeter)/sqr(omega/oneMeter)));
      expY = exp(C_d*(sqr((cellCentre.component(vector::Y)-(alphaPos.component(vector::Y)))/oneMeter)/sqr(omega/oneMeter)));
      expZ = exp(C_d*(sqr((cellCentre.component(vector::Z)-(alphaPos.component(vector::Z)))/oneMeter)/sqr(omega/oneMeter)));
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
Foam::alphaSourceClass::genField(const dimensionedVector alphaPos, const dimensionedScalar omega)
{
	volScalarField genField
	(
		geoField(alphaPos, omega)
	);

    forAll(genField, i)
    {
        (genField[i] >= alphaFieldCut.value()) ? genField[i] = scalar(1) : genField[i] = scalar(0);
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
Foam::alphaSourceClass::normField(const volScalarField& alphaField)
{
	scalar fieldMax = max(alphaField).value();	
	
	volScalarField alphaFieldNorm
	(
		alphaField/fieldMax
	);

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "alphaFieldNorm",
			alphaFieldNorm
    	)
  	);
}

Foam::tmp<Foam::volScalarField>
Foam::alphaSourceClass::binField(const volScalarField& alphaField)
{
    volScalarField alphaFieldBin(alphaField);

    forAll(alphaFieldBin, i)
    {
        if (alphaFieldBin[i] >= 1.0)
            alphaFieldBin[i] = 1.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "alphaFieldBin",
            alphaFieldBin
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::alphaSourceClass::floorField(const volScalarField& alphaField)
{
    volScalarField alphaFieldFloor(alphaField);

    forAll(alphaFieldFloor, i)
    {
        if (alphaFieldFloor[i] <= 1.0)
            alphaFieldFloor[i] = 0.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "alphaFieldFloor",
            alphaFieldFloor
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::alphaSourceClass::applyField()
{
    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    Info << " Generating Heat Alpha Field " << endl;
    alphaFieldSize = 0.0;

    if ((U_.mesh().time() > alphaStartTime) && (U_.mesh().time() < alphaEndTime)) 
    	alphaFieldSize = binField(genField(alphaPos0, alphaOmega));

    if (U_.mesh().time() > alphaEndTime)
      alphaFieldSize = 0.0;

    scalar heatAlphaFieldMax = max(alphaFieldSize).value();
    scalar heatAlphaFieldMin = min(alphaFieldSize).value();
    Info << "heatAlphaField max: " << heatAlphaFieldMax 
    << " min: " << heatAlphaFieldMin << endl;

    if (!alphaOn)
        alphaMult = 0.0;

    volScalarField alpha0 = normField(geoField(alphaPos0, alphaOmega));    

    scalar alpha0Max = max(alpha0).value();
    scalar alpha0Min = min(alpha0).value();
    Info <<  "alpha0 max: " << alpha0Max << " min: " << alpha0Min << endl;

    alphaFieldOutput = alphaFieldSize*alphaMult;

    if (alphaMod)
        alphaFieldOutput = alphaFieldSize*alpha0*alphaMult;

    return alphaFieldOutput;
}