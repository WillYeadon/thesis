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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 276893b7ceb62a059a9a0eb22e4c2c39823c2afe
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void swirlster_276893b7ceb62a059a9a0eb22e4c2c39823c2afe(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    swirlsterFixedValueFvPatchScalarField
);


const char* const swirlsterFixedValueFvPatchScalarField::SHA1sum =
    "276893b7ceb62a059a9a0eb22e4c2c39823c2afe";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swirlsterFixedValueFvPatchScalarField::
swirlsterFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe"
            " from patch/DimensionedField\n";
    }
}


swirlsterFixedValueFvPatchScalarField::
swirlsterFixedValueFvPatchScalarField
(
    const swirlsterFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe"
            " from patch/DimensionedField/mapper\n";
    }
}


swirlsterFixedValueFvPatchScalarField::
swirlsterFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe"
            " from patch/dictionary\n";
    }
}


swirlsterFixedValueFvPatchScalarField::
swirlsterFixedValueFvPatchScalarField
(
    const swirlsterFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe"
            " as copy\n";
    }
}


swirlsterFixedValueFvPatchScalarField::
swirlsterFixedValueFvPatchScalarField
(
    const swirlsterFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

swirlsterFixedValueFvPatchScalarField::
~swirlsterFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void swirlsterFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs swirlster sha1: 276893b7ceb62a059a9a0eb22e4c2c39823c2afe\n";
    }

//{{{ begin code
    #line 30 "/home/will/OpenFOAM/will-6/run/krade-weld-0/0/alpha.Phase1.boundaryField.hot"
const fvPatch& boundaryPatch = patch();
         const vectorField& Cf = boundaryPatch.Cf();
         scalarField& field = *this;

         field = patchInternalField(); 
            forAll(Cf, facei)
            {
                field[facei] = 0.5;
            }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

