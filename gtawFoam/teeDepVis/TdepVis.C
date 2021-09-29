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

#include "TdepVis.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(TdepVis, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        TdepVis,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
/*
Foam::tmp<Foam::dimensionedScalar>
Foam::viscosityModels::TdepVis::calcNuM() const
{
	return nuM_/exp(2.65*pow(Tmelt_, 1.27)/R_*Tmelt_);
}
*/
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::TdepVis::calcNu() const
{

	const volScalarField& T_ = U_.mesh().lookupObject<volScalarField>("T");

    const volScalarField outputNu
    (
//    	outputNuV()
		tuner_*(nu0V_/rho_)*exp(2.65*pow((max(T_, smallT_)/Tunits_), 0.27)/R_)
    );

	Info<< "calcNu: " << outputNu().weightedAverage(U_.mesh().V()) << endl;

	return tuner_*(nu0V_/rho_)*exp(2.65*pow((max(T_, smallT_)/Tunits_), 0.27)/R_);

/*
    return max
    (
        nuMin_,
        exp(2.65*pow(T, 1.27)/R_*T)*min
        (
            nuMax_,
            k_*pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("small", dimless, small)
                ),
                n_.value() - scalar(1)
            )
        )
    );
*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::TdepVis::TdepVis
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
//    threePhaseMixture(U.mesh(), *this),

    viscosityModel(name, viscosityProperties, U, phi),
    TdepVisCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
/*
    rho_("rho", dimensionSet(1,-3,0,0,0,0,0), TdepVisCoeffs_),
    nuM_("nuM", dimensionSet(1,1,-3,-1,0,0,0), TdepVisCoeffs_),
    Tmelt_("Tmelt", dimensionSet(0,0,0,1,0,0,0), TdepVisCoeffs_),
//    R_("R", dimensionSet(1,2,-2,-1,-1,0,0), 8.31446),
    R_("R", dimensionSet(1,2,-2,-1,0,0,0), 8.31446),
*/
    rho_("rho", dimless, TdepVisCoeffs_),
    nuM_("nuM", dimless, TdepVisCoeffs_),
    Tmelt_("Tmelt", dimless, TdepVisCoeffs_),
	Tunits_("Tunits", dimensionSet(0,0,0,1,0,0,0), 1),
	smallT_("smallT", dimensionSet(0,0,0,1,0,0,0), 1e-6),
	tuner_("tuner", dimensionSet(0,2,-1,0,0,0,0), 1),
//    R_("R", dimensionSet(1,2,-2,-1,-1,0,0), 8.31446),
    R_("R", dimless, 8.31446),
	nu0V_
    (
        IOobject
        (
            "nu0V",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
//        nuM_,
//        (nuM_/exp(2.65*pow(Tmelt_, 1.27)/R_*Tmelt_)),
        (nuM_/exp(2.65*pow(Tmelt_, 0.27)/R_)),
        calculatedFvPatchScalarField::typeName
    ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::TdepVis::outputNuV() const
{
	const volScalarField& T = U_.mesh().lookupObject<volScalarField>("T");

	return tuner_*(nu0V_/rho_)*exp(2.65*pow((T/Tunits_), 0.27)/R_);
}
*/
bool Foam::viscosityModels::TdepVis::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    TdepVisCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    TdepVisCoeffs_.lookup("rho") >> rho_;
    TdepVisCoeffs_.lookup("Tmelt") >> Tmelt_;
    TdepVisCoeffs_.lookup("nuM") >> nuM_;
// /    TdepVisCoeffs_.lookup("nuMax") >> nuMax_;

/*
    TdepVisCoeffs_.lookup("k") >> k_;
    TdepVisCoeffs_.lookup("n") >> n_;
    TdepVisCoeffs_.lookup("nuMin") >> nuMin_;
    TdepVisCoeffs_.lookup("nuMax") >> nuMax_;
*/
    return true;
}


// ************************************************************************* //
