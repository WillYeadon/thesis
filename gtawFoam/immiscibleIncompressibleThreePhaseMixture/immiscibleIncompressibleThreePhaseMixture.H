/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::immiscibleIncompressibleThreePhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    immiscibleIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleThreePhaseMixture_H
#define immiscibleIncompressibleThreePhaseMixture_H

#include "gtawIncompressibleThreePhaseMixture/gtawIncompressibleThreePhaseMixture.H"
#include "gtawInterfaceProperties/gtawInterfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class immiscibleIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleThreePhaseMixture
:
    public gtawIncompressibleThreePhaseMixture,
    public gtawInterfaceProperties
{

public:

    // Constructors

        //- Construct from components
        immiscibleIncompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi//,
//            const volScalarField& T,
//            const dimensionedScalar& Tref,
//            const dimensionedScalar& beta//9,
//9            const word& alpha1Name = "alpha1",
//9            const word& alpha2Name = "alpha2",
//9            const word& alpha3Name = "alpha3"
        /*
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const volScalarField& alpha3
        */
        );


    //- Destructor
    virtual ~immiscibleIncompressibleThreePhaseMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            gtawIncompressibleThreePhaseMixture::correct();
            gtawInterfaceProperties::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
