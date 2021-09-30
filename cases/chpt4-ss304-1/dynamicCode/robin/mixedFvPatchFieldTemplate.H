/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with dynamic code generation of a
    mixed fvPatchField.

    - without state

SourceFiles
    mixedFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef mixedFvPatchTemplateScalarField_H
#define mixedFvPatchTemplateScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     A templated MixedValueFvPatchField
\*---------------------------------------------------------------------------*/

class robinMixedValueFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("robin");


    // Constructors

        //- Construct from patch and internal field
        robinMixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        robinMixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        robinMixedValueFvPatchScalarField
        (
            const robinMixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        robinMixedValueFvPatchScalarField
        (
            const robinMixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField > clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new robinMixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        robinMixedValueFvPatchScalarField
        (
            const robinMixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField >
            (
                new robinMixedValueFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~robinMixedValueFvPatchScalarField();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

