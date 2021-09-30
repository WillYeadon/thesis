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
    // SHA1 = d48b2d67c3c4946f301c470d8c7d7d97d76c1978
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void robin_d48b2d67c3c4946f301c470d8c7d7d97d76c1978(bool load)
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
    robinFixedValueFvPatchScalarField
);


const char* const robinFixedValueFvPatchScalarField::SHA1sum =
    "d48b2d67c3c4946f301c470d8c7d7d97d76c1978";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

robinFixedValueFvPatchScalarField::
robinFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978"
            " from patch/DimensionedField\n";
    }
}


robinFixedValueFvPatchScalarField::
robinFixedValueFvPatchScalarField
(
    const robinFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978"
            " from patch/DimensionedField/mapper\n";
    }
}


robinFixedValueFvPatchScalarField::
robinFixedValueFvPatchScalarField
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
        Info<<"construct robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978"
            " from patch/dictionary\n";
    }
}


robinFixedValueFvPatchScalarField::
robinFixedValueFvPatchScalarField
(
    const robinFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978"
            " as copy\n";
    }
}


robinFixedValueFvPatchScalarField::
robinFixedValueFvPatchScalarField
(
    const robinFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

robinFixedValueFvPatchScalarField::
~robinFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void robinFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs robin sha1: d48b2d67c3c4946f301c470d8c7d7d97d76c1978\n";
    }

//{{{ begin code
    #line 30 "/home/will/OpenFOAM/will-6/run/krade-weld-0/0/T.boundaryField.top"
scalarField field(this->patch().Cf() * 2.0);
    	       field = 320;
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

