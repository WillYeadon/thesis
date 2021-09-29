/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
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

	Class to create up to three heat sources with a elliptic paraboloid shape
	Modified from http://dx.doi.org/10.1016/j.ijheatmasstransfer.2016.04.064

\*---------------------------------------------------------------------------*/

#include "emClass/emSourceClass.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

Foam::emSourceClass::emSourceClass(
        const dictionary& emProperties,
        const volVectorField& U,
        const surfaceScalarField& phi
):

    emPropertiesDict(emProperties),
    U_(U),
    phi_(phi),
    
    sourceVel("sourceVel", dimLength, emPropertiesDict),
    sourcePos0("sourcePos0", dimLength, emPropertiesDict.lookupOrDefault<vector>("sourcePos0", Zero)),
    sourceOmega("sourceOmega", dimLength, emPropertiesDict),
    sourcePenn("sourcePenn", dimLength, emPropertiesDict),
    sourceStartTime("sourceStartTime", dimTime, emPropertiesDict),
    sourcePauseTime("sourcePauseTime", dimTime, emPropertiesDict), 
    sourceEndTime("sourceEndTime", dimTime, emPropertiesDict), 
    sourceB0("sourceB0", dimLength, emPropertiesDict),
    
    sourceC_C("sourceC_C", dimless, emPropertiesDict),
    sourceFieldCut("sourceFieldCut", dimless, emPropertiesDict),

    sourcePause(emPropertiesDict.lookup("sourcePause")),
    sourceOn(emPropertiesDict.lookup("sourceOn")), 

    emFieldSize
    (
        IOobject
        (
            "emFieldSize",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("emFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    emFieldOutput
    (
        IOobject
        (
            "emFieldOutput",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedVector("emFieldOutput", dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{
  
}

bool Foam::emSourceClass::sourceOnCheck()
{
    return (sourceOn ? true : false);
}

Foam::tmp<Foam::volScalarField>
Foam::emSourceClass::geoField(const dimensionedVector sourcePos)
{
    const volVectorField& cellCentre = U_.mesh().C();

    scalar C_d = log(sourceFieldCut.value());

    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    scalar currentTime = U_.mesh().time().value();

    bool sourcePauseSwitch;
    (currentTime < sourcePauseTime.value()) ? sourcePauseSwitch = true : sourcePauseSwitch = false;

    volScalarField expX
    (
        exp(C_d*(sqr((cellCentre.component(vector::X)-(sourcePos.component(vector::X))
        + (sourceVel.component(vector::X)*(currentTime  - sourcePauseTime.value())))/oneMeter)/sqr(sourceOmega/oneMeter)))
    );

    volScalarField expY
    (
        exp(C_d*(sqr((cellCentre.component(vector::Y)-(sourcePos.component(vector::Y))
        + (sourceVel.component(vector::Y)*(currentTime  - sourcePauseTime.value())))/oneMeter)/sqr(sourceOmega/oneMeter)))
    ); 

    volScalarField expZ
    (
        exp(C_d*(sqr((cellCentre.component(vector::Z)-(sourcePos.component(vector::Z)))/oneMeter)/sqr(sourcePenn/oneMeter)))
    );


    if (sourcePause && sourcePauseSwitch)      
    {
      expX = exp(C_d*(sqr((cellCentre.component(vector::X)-(sourcePos.component(vector::X)))/oneMeter)/sqr(sourceOmega/oneMeter)));
      expY = exp(C_d*(sqr((cellCentre.component(vector::Y)-(sourcePos.component(vector::Y)))/oneMeter)/sqr(sourceOmega/oneMeter)));
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
Foam::emSourceClass::genField(const dimensionedVector sourcePos)
{
	volScalarField genField
	(
		geoField(sourcePos)
	);

    forAll(genField, i)
    {
        (genField[i] >= sourceFieldCut.value()) ? genField[i] = scalar(1) : genField[i] = scalar(0);
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
Foam::emSourceClass::normField(const volScalarField& sourceField)
{
	scalar fieldMax = max(sourceField).value();	
	
	volScalarField sourceFieldNorm
	(
		sourceField/fieldMax
	);

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "sourceFieldNorm",
			sourceFieldNorm
    	)
  	);
}

Foam::tmp<Foam::volVectorField>
Foam::emSourceClass::applyField()
{
    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    Info << " Generating EM Source Field " << endl;
    emFieldSize = 0.0;

    if ((sourceStartTime < U_.mesh().time()) && (U_.mesh().time() < sourceEndTime)) 
    {
    	emFieldSize = genField(sourcePos0);
    }

    if ((U_.mesh().time() > sourceEndTime) || (U_.mesh().time() < sourceStartTime))
    {
        emFieldSize = 0.0;
    }

    scalar emSourceFieldMax = max(emFieldSize).value();
    scalar emSourceFieldMin = min(emFieldSize).value();
    Info << "emSourceField max: " << emSourceFieldMax 
    << " min: " << emSourceFieldMin << endl;

    volScalarField source0 = normField(geoField(sourcePos0));
    
    scalar source0Max = max(source0).value();
	scalar source0Min = min(source0).value();

    Info <<  "emSource0 max: " << source0Max 
    << " min: " << source0Min << endl;

    emFieldOutput = (emFieldSize*source0)*sourceB0;

    if (!sourceOn)  
        emFieldOutput *= 0.0;            

    return emFieldOutput;
}