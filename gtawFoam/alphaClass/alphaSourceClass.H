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

#ifndef alphaSourceClass_H
#define alphaSourceClass_H

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

class alphaSourceClass
{

        const dictionary& alphaPropertiesDict;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dimensionedVector alphaVel;
        dimensionedVector alphaPos0;
        dimensionedScalar alphaOmega;
        dimensionedScalar alphaMult;

        dimensionedScalar alphaStartTime;
        dimensionedScalar alphaPauseTime;
        dimensionedScalar alphaEndTime;

        dimensionedScalar alphaC_C;
        dimensionedScalar alphaFieldCut;

        bool alphaPause;
        bool alphaOn;
        bool alphaMod;        

        volScalarField alphaFieldSize;
//        volScalarField alphaModSize;
        volScalarField alphaFieldOutput;

public:

	alphaSourceClass(
	        const dictionary& alphaProperties,
	        const volVectorField& U,
	        const surfaceScalarField& phi
		);

	virtual ~alphaSourceClass(){}

    virtual tmp<volScalarField> geoField(const dimensionedVector alphaPos, const dimensionedScalar omega);
    virtual tmp<volScalarField> genField(const dimensionedVector alphaPos, const dimensionedScalar omega);
    virtual tmp<volScalarField> normField(const volScalarField& alphaField);
    virtual tmp<volScalarField> binField(const volScalarField& alphaField);
    virtual tmp<volScalarField> floorField(const volScalarField& alphaField);
    virtual tmp<volScalarField> applyField();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //