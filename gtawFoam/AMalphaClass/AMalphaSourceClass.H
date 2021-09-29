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

\*---------------------------------------------------------------------------*/

#ifndef AMalphaSourceClass_H
#define AMalphaSourceClass_H

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

class AMalphaSourceClass
{

        const dictionary& AMalphaPropertiesDict;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dimensionedVector AMalphaVel;
        dimensionedVector AMalphaPos0;
        dimensionedScalar AMalphaOmega;
        dimensionedScalar AMalphaMult;

        dimensionedScalar AMalphaStartTime;
        dimensionedScalar AMalphaPauseTime;
        dimensionedScalar AMalphaEndTime;

        dimensionedScalar AMalphaC_C;
        dimensionedScalar AMalphaFieldCut;

        bool AMalphaPause;
        bool AMalphaOn;
        bool AMalphaMod;        

        volScalarField AMalphaFieldSize;
//        volScalarField AMalphaModSize;
        volScalarField AMalphaFieldOutput;

public:

	AMalphaSourceClass(
	        const dictionary& AMalphaProperties,
	        const volVectorField& U,
	        const surfaceScalarField& phi
		);

	virtual ~AMalphaSourceClass(){}

    virtual tmp<volScalarField> geoField(const dimensionedVector AMalphaPos, const dimensionedScalar omega);
    virtual tmp<volScalarField> genField(const dimensionedVector AMalphaPos, const dimensionedScalar omega);
    virtual tmp<volScalarField> normField(const volScalarField& AMalphaField);
    virtual tmp<volScalarField> binField(const volScalarField& AMalphaField);
    virtual tmp<volScalarField> floorField(const volScalarField& AMalphaField);
    virtual tmp<volScalarField> applyField();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //