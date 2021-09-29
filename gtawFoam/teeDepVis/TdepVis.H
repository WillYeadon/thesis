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

Class
    Foam::viscosityModels::TdepVis

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    TdepVis.C

\*---------------------------------------------------------------------------*/

#ifndef TdepVis_H
#define TdepVis_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class TdepVis Declaration
\*---------------------------------------------------------------------------*/

class TdepVis
:
    public viscosityModel
{
    // Private data

        dictionary TdepVisCoeffs_;

        dimensionedScalar rho_;
        dimensionedScalar nuM_;
        dimensionedScalar Tmelt_;
        dimensionedScalar Tunits_;
        dimensionedScalar smallT_;
        dimensionedScalar tuner_;
        dimensionedScalar R_;
        dimensionedScalar nu0_;
//        dimensionedScalar nuMin_;
//        dimensionedScalar nuMax_;

        volScalarField nu0V_;
        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
//        tmp<Foam::dimensionedScalar> calcNuM() const;
        tmp<volScalarField> calcNu() const;

public:

    //- Runtime type information
    TypeName("TdepVis");


    // Constructors

        //- Construct from components
        TdepVis
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~TdepVis()
    {}


    // Member Functions

    // returns reference to temporary
/*
        const dimensionedScalar& nu0() const
        {
            return nuM_/exp(2.65*pow(Tmelt_, 1.27)/R_*Tmelt_);
        }
*/
        //- Return the laminar viscosity

//        tmp<volScalarField> outputNuV() const;

        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
