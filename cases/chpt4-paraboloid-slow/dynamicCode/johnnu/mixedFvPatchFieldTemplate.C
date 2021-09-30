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
#line 36 "/home/will/OpenFOAM/will-6/run/krade-weld-0/0/T.boundaryField.top"
//         #include "fluidThermo.H"
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
    // SHA1 = 4ae3693dfad00ff072a5737c97689db1a74e7879
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void johnnu_4ae3693dfad00ff072a5737c97689db1a74e7879(bool load)
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
    johnnuMixedValueFvPatchScalarField
);


const char* const johnnuMixedValueFvPatchScalarField::SHA1sum =
    "4ae3693dfad00ff072a5737c97689db1a74e7879";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

johnnuMixedValueFvPatchScalarField::
johnnuMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879"
            " from patch/DimensionedField\n";
    }
}


johnnuMixedValueFvPatchScalarField::
johnnuMixedValueFvPatchScalarField
(
    const johnnuMixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879"
            " from patch/DimensionedField/mapper\n";
    }
}


johnnuMixedValueFvPatchScalarField::
johnnuMixedValueFvPatchScalarField
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
        Info<<"construct johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879"
            " from patch/dictionary\n";
    }
}


johnnuMixedValueFvPatchScalarField::
johnnuMixedValueFvPatchScalarField
(
    const johnnuMixedValueFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879"
            " as copy\n";
    }
}


johnnuMixedValueFvPatchScalarField::
johnnuMixedValueFvPatchScalarField
(
    const johnnuMixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

johnnuMixedValueFvPatchScalarField::
~johnnuMixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void johnnuMixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs johnnu sha1: 4ae3693dfad00ff072a5737c97689db1a74e7879\n";
    }

//{{{ begin code
    #line 47 "/home/will/OpenFOAM/will-6/run/krade-weld-0/0/T.boundaryField.top"
const scalarField& delta = patch().deltaCoeffs();
          const fvMesh& mesh = patch().boundaryMesh().mesh();
//          const fluidThermo& thermo = mesh.lookupObject<fluidThermo>(basicThermo::dictName);
          const label patchi = patch().index();
          this->refValue() = 298;
          this->refGrad() = 0;
          this->valueFraction() = 1.0/(1.0+0.6/(1.0/delta));
//}}} end code

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

