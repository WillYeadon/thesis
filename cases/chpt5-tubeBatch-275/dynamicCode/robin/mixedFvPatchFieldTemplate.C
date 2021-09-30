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

#include "mixedFvPatchFieldTemplate.H"
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
    // SHA1 = fa6bcb9352550c3007ffd513d48777a6a88042f2
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void robin_fa6bcb9352550c3007ffd513d48777a6a88042f2(bool load)
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
    robinMixedValueFvPatchScalarField
);


const char* const robinMixedValueFvPatchScalarField::SHA1sum =
    "fa6bcb9352550c3007ffd513d48777a6a88042f2";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

robinMixedValueFvPatchScalarField::
robinMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2"
            " from patch/DimensionedField\n";
    }
}


robinMixedValueFvPatchScalarField::
robinMixedValueFvPatchScalarField
(
    const robinMixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2"
            " from patch/DimensionedField/mapper\n";
    }
}


robinMixedValueFvPatchScalarField::
robinMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2"
            " from patch/dictionary\n";
    }
}


robinMixedValueFvPatchScalarField::
robinMixedValueFvPatchScalarField
(
    const robinMixedValueFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2"
            " as copy\n";
    }
}


robinMixedValueFvPatchScalarField::
robinMixedValueFvPatchScalarField
(
    const robinMixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

robinMixedValueFvPatchScalarField::
~robinMixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void robinMixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs robin sha1: fa6bcb9352550c3007ffd513d48777a6a88042f2\n";
    }

//{{{ begin code
    
//}}} end code

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

