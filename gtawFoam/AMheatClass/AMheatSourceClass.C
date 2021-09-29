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

	Class to create up to three heat AMsources with a elliptic paraboloid shape
	Modified from http://dx.doi.org/10.1016/j.ijheatmasstransfer.2016.04.064

\*---------------------------------------------------------------------------*/

#include "AMheatClass/AMheatSourceClass.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

Foam::AMheatSourceClass::AMheatSourceClass(
        const dictionary& AMsourceProperties,
        const volScalarField& alphaMetalTotal,
        const volVectorField& U,
        const surfaceScalarField& phi
):

    AMsourcePropertiesDict(AMsourceProperties),
    alphaMetalTotal_(alphaMetalTotal),
    U_(U),
    phi_(phi),
    
    AMsourceVel("AMsourceVel", dimLength, AMsourcePropertiesDict),
    AMsourcePos0("AMsourcePos0", dimLength, AMsourcePropertiesDict.lookupOrDefault<vector>("AMsourcePos0", Zero)),
    AMsourcePos1("AMsourcePos1", dimLength, AMsourcePropertiesDict.lookupOrDefault<vector>("AMsourcePos1", Zero)),
    AMsourcePos2("AMsourcePos2", dimLength, AMsourcePropertiesDict.lookupOrDefault<vector>("AMsourcePos2", Zero)),
    AMsourceOmega("AMsourceOmega", dimLength, AMsourcePropertiesDict),
    AMsourcePenn("AMsourcePenn", dimLength, AMsourcePropertiesDict),
    AMsourceModOmega("AMsourceModOmega", dimLength, AMsourcePropertiesDict.lookupOrDefault<scalar>("AMsourceModOmega", 1e-3)),
    AMsourceModPenn("AMsourceModPenn", dimLength, AMsourcePropertiesDict.lookupOrDefault<scalar>("AMsourceModPenn", 1e-3)),
    AMsourceModFraction("AMsourceModFraction", dimless, AMsourcePropertiesDict.lookupOrDefault<scalar>("AMsourceModFraction", 0.1)),
    AMsourceStartTime("AMsourceStartTime", dimTime, AMsourcePropertiesDict),
    AMsourcePauseTime("AMsourcePauseTime", dimTime, AMsourcePropertiesDict), 
    AMsourceEndTime("AMsourceEndTime", dimTime, AMsourcePropertiesDict), 
    AMsourcePower("AMsourcePower", dimLength, AMsourcePropertiesDict),
    
    AMsourceC_C("AMsourceC_C", dimless, AMsourcePropertiesDict),
    AMsourceFieldCut("AMsourceFieldCut", dimless, AMsourcePropertiesDict),

    AMsourcePause(AMsourcePropertiesDict.lookup("AMsourcePause")),
    AMsourceMod(AMsourcePropertiesDict.lookupOrDefault<bool>("AMsourceMod", false)),  
    AMsourceOn(AMsourcePropertiesDict.lookupOrDefault<bool>("AMsourceOn", true)),
    AMsourceTwo(AMsourcePropertiesDict.lookupOrDefault<bool>("AMsourceTwo", false)),
    AMsourceThree(AMsourcePropertiesDict.lookupOrDefault<bool>("AMsourceThree", false)),  

    AMsourceFieldSize
    (
        IOobject
        (
            "AMsourceFieldSize",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMsourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    AMsourceFieldSizeTwo
    (
        IOobject
        (
            "AMsourceFieldSizeTwo",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMsourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    AMsourceFieldSizeThree
    (
        IOobject
        (
            "AMsourceFieldSizeThree",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMsourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    AMsourceFieldOutput
    (
        IOobject
        (
            "AMsourceFieldOutput",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("AMsourceFieldOutput", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
  
}

Foam::tmp<Foam::volScalarField>
Foam::AMheatSourceClass::geoField(const dimensionedVector AMsourcePos, const dimensionedScalar omega, const dimensionedScalar penn)
{
    const volVectorField& cellCentre = U_.mesh().C();

    scalar C_d = log(AMsourceFieldCut.value());

    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    scalar currentTime = U_.mesh().time().value();

    bool AMsourcePauseSwitch;
    (currentTime < AMsourcePauseTime.value()) ? AMsourcePauseSwitch = true : AMsourcePauseSwitch = false;

    volScalarField expX
    (
        exp(C_d*(sqr((cellCentre.component(vector::X)-(AMsourcePos.component(vector::X))
        + (AMsourceVel.component(vector::X)*(currentTime  - AMsourcePauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    );

    volScalarField expY
    (
        exp(C_d*(sqr((cellCentre.component(vector::Y)-(AMsourcePos.component(vector::Y))
        + (AMsourceVel.component(vector::Y)*(currentTime  - AMsourcePauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 

    volScalarField expZ
    (
      exp(C_d*sqrt(sqr((cellCentre.component(vector::Z) - (AMsourcePos.component(vector::Z)))/oneMeter)) / (penn/oneMeter))
    );


    if (AMsourcePause && AMsourcePauseSwitch)      
    {
      expX = exp(C_d*(sqr((cellCentre.component(vector::X)-(AMsourcePos.component(vector::X)))/oneMeter)/sqr(omega/oneMeter)));
      expY = exp(C_d*(sqr((cellCentre.component(vector::Y)-(AMsourcePos.component(vector::Y)))/oneMeter)/sqr(omega/oneMeter)));
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
Foam::AMheatSourceClass::genField(const volScalarField& alphaMetalTotal, const dimensionedVector AMsourcePos, const dimensionedScalar omega, const dimensionedScalar penn)
{
	volScalarField genField
	(
		geoField(AMsourcePos, omega, penn)*alphaMetalTotal
	);

    forAll(genField, i)
    {
        ((genField[i] >= AMsourceFieldCut.value()) && (alphaMetalTotal[i] >= AMsourceFieldCut.value())) ? genField[i] = scalar(1) : genField[i] = scalar(0);
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
Foam::AMheatSourceClass::normField(const volScalarField& AMsourceField)
{
	scalar fieldMax = max(AMsourceField).value();	
	
	volScalarField AMsourceFieldNorm
	(
		AMsourceField/fieldMax
	);

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMsourceFieldNorm",
			AMsourceFieldNorm
    	)
  	);
}

Foam::tmp<Foam::volScalarField>
Foam::AMheatSourceClass::binField(const volScalarField& AMsourceField)
{
    volScalarField AMsourceFieldBin(AMsourceField);

    forAll(AMsourceFieldBin, i)
    {
        if (AMsourceFieldBin[i] >= 1.0)
            AMsourceFieldBin[i] = 1.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMsourceFieldBin",
            AMsourceFieldBin
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::AMheatSourceClass::floorField(const volScalarField& AMsourceField)
{
    volScalarField AMsourceFieldFloor(AMsourceField);

    forAll(AMsourceFieldFloor, i)
    {
        if (AMsourceFieldFloor[i] <= 1.0)
            AMsourceFieldFloor[i] = 0.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "AMsourceFieldFloor",
            AMsourceFieldFloor
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::AMheatSourceClass::applyField(const volScalarField& alphaMetalTotal)
{
    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    Info << " Q_arc " << endl;
    scalar Q_arc = AMsourcePower.value();
    Info << " Generating Heat Source Field " << endl;
    AMsourceFieldSize = 0.0;
    AMsourceFieldSizeTwo = 0.0;
    AMsourceFieldSizeThree = 0.0;

    volScalarField AMsourceModField0 = genField(alphaMetalTotal, AMsourcePos0, AMsourceModOmega, AMsourceModPenn);
    volScalarField AMsourceModField1 = genField(alphaMetalTotal, AMsourcePos1, AMsourceModOmega, AMsourceModPenn);
    volScalarField AMsourceModField2 = genField(alphaMetalTotal, AMsourcePos2, AMsourceModOmega, AMsourceModPenn);

    /*
    Source Mod allows you to overlay two fields to make a sombrero shape 
    */

    if (!AMsourceMod)
    {
        AMsourceModField0 *= 0.0;
        AMsourceModField1 *= 0.0;
        AMsourceModField2 *= 0.0;
    }

    if (AMsourceOn && U_.mesh().time() < AMsourceEndTime) 
    {
    	AMsourceFieldSize = binField(genField(alphaMetalTotal, AMsourcePos0, AMsourceOmega, AMsourcePenn) + AMsourceModField0);
    }
    if (AMsourceTwo && (U_.mesh().time() < AMsourceEndTime)) 
    {
    	AMsourceFieldSizeTwo = binField(genField(alphaMetalTotal, AMsourcePos1, AMsourceOmega, AMsourcePenn) + AMsourceModField1);
    }
    if (AMsourceThree && (U_.mesh().time() < AMsourceEndTime)) 
    {
    	AMsourceFieldSizeThree = binField(genField(alphaMetalTotal, AMsourcePos2, AMsourceOmega, AMsourcePenn) + AMsourceModField2);
    }

    if (U_.mesh().time() > AMsourceEndTime)
    {
      AMsourceFieldSize = 0.0;
      AMsourceFieldSizeTwo = 0.0;
      AMsourceFieldSizeThree = 0.0;
    }

    scalar heatSourceFieldMax = max(AMsourceFieldSize).value();
    scalar heatSourceFieldMin = min(AMsourceFieldSize).value();
    Info << "heatSourceField max: " << heatSourceFieldMax 
    << " min: " << heatSourceFieldMin << endl;

    Info << " delHField " << endl;
    volScalarField delHField = ((AMsourceC_C * Q_arc) / (Foam::sqr(AMsourceOmega/oneMeter)*Foam::sqrt(AMsourcePenn/oneMeter)))*alphaMetalTotal;

    volScalarField AMsource0 = normField(geoField(AMsourcePos0, AMsourceOmega, AMsourcePenn)); 
	volScalarField AMsource1 = normField(geoField(AMsourcePos1, AMsourceOmega, AMsourcePenn)); 
    volScalarField AMsource2 = normField(geoField(AMsourcePos2, AMsourceOmega, AMsourcePenn)); 
    
    if (AMsourceMod)
    {
        AMsource0 = normField(geoField(AMsourcePos0, AMsourceOmega, AMsourcePenn) + AMsourceModFraction*geoField(AMsourcePos0, AMsourceModOmega, AMsourceModPenn));
        AMsource1 = normField(geoField(AMsourcePos1, AMsourceOmega, AMsourcePenn) + AMsourceModFraction*geoField(AMsourcePos0, AMsourceModOmega, AMsourceModPenn));
        AMsource2 = normField(geoField(AMsourcePos2, AMsourceOmega, AMsourcePenn) + AMsourceModFraction*geoField(AMsourcePos0, AMsourceModOmega, AMsourceModPenn));
    }

    scalar AMsource0Max = max(AMsource0).value();
    scalar AMsource0Min = min(AMsource0).value();
    Info <<  "AMsource0 max: " << AMsource0Max << " min: " << AMsource0Min << endl;

    AMsourceFieldOutput = (AMsourceFieldSize*AMsource0 + AMsourceFieldSizeTwo*AMsource1 + AMsourceFieldSizeThree*AMsource2)*delHField;

    return AMsourceFieldOutput;
}