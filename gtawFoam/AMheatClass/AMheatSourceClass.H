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

\*---------------------------------------------------------------------------*/

#ifndef AMheatSourceClass_H
#define AMheatSourceClass_H

#include "uniformDimensionedFields.H"
#include "dimensionedTypes.H"
#include "Time.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "fvMesh.H"

namespace Foam
{

class AMheatSourceClass
{

        const dictionary& AMsourcePropertiesDict;
        const volScalarField& alphaMetalTotal_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dimensionedVector AMsourceVel;
        dimensionedVector AMsourcePos0;
        dimensionedVector AMsourcePos1;
        dimensionedVector AMsourcePos2;
        dimensionedScalar AMsourceOmega;
        dimensionedScalar AMsourcePenn;
        dimensionedScalar AMsourceModOmega;
        dimensionedScalar AMsourceModPenn;
        dimensionedScalar AMsourceModFraction;
        dimensionedScalar AMsourceStartTime;
        dimensionedScalar AMsourcePauseTime;
        dimensionedScalar AMsourceEndTime;
        dimensionedScalar AMsourcePower;

        dimensionedScalar AMsourceC_C;
        dimensionedScalar AMsourceFieldCut;

        bool AMsourcePause;
        bool AMsourceMod;
        bool AMsourceOn;        
        bool AMsourceTwo;
        bool AMsourceThree;

        volScalarField AMsourceFieldSize;
//        volScalarField AMsourceModSize;
        volScalarField AMsourceFieldSizeTwo;
        volScalarField AMsourceFieldSizeThree;
        volScalarField AMsourceFieldOutput;

public:

	AMheatSourceClass(
	        const dictionary& AMsourceProperties,
            const volScalarField& alphaMetalTotal,
	        const volVectorField& U,
	        const surfaceScalarField& phi
		);

	virtual ~AMheatSourceClass(){}

    virtual tmp<volScalarField> geoField(const dimensionedVector AMsourcePos, const dimensionedScalar omega, const dimensionedScalar penn);
    virtual tmp<volScalarField> genField(const volScalarField& alphaMetalTotal, const dimensionedVector AMsourcePos, const dimensionedScalar omega, const dimensionedScalar penn);
    virtual tmp<volScalarField> normField(const volScalarField& AMsourceField);
    virtual tmp<volScalarField> binField(const volScalarField& AMsourceField);
    virtual tmp<volScalarField> floorField(const volScalarField& AMsourceField);
    virtual tmp<volScalarField> applyField(const volScalarField& alphaMetalTotal);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //