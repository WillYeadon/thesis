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

#include "gtawInterfaceProperties.H"
#include "groundAlphaContact/groundAlphaContactAngleFvPatchScalarField.H"
//#include "weldPhaseProperties/alphaContactAngle/alphaContactAngle/alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::gtawInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::gtawInterfaceProperties::correctContactAngle12
(
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& abf = alpha1_.boundaryField();
    const volVectorField::Boundary& U = U_.boundaryField();

    const fvMesh& mesh = U_.mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            Info << "Bing Bang Bong" << endl;
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

//            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
//            acap.evaluate();
        }
    }
}

void Foam::gtawInterfaceProperties::correctContactAngle13
(
//    const volScalarField& alpha1,
//    const volScalarField& alpha3,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& abf = alpha2_.boundaryField();
//        mixture_.alpha1().boundaryField();

    const volVectorField::Boundary& U = U_.boundaryField();
//        mixture_.U().boundaryField();

    const fvMesh& mesh = U_.mesh();
//    const fvMesh& mesh = mixture_.U().mesh();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

//            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
//            acap.evaluate();
        }
    }
}

void Foam::gtawInterfaceProperties::correctContactAngle23
(
//    const volScalarField& alpha2,
//    const volScalarField& alpha3,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& abf = alpha3_.boundaryField();
//        mixture_.alpha1().boundaryField();

    const volVectorField::Boundary& U = U_.boundaryField();
//        mixture_.U().boundaryField();

    const fvMesh& mesh = U_.mesh();
//    const fvMesh& mesh = mixture_.U().mesh();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

//            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
//            acap.evaluate();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gtawInterfaceProperties::gtawInterfaceProperties
(
	const volScalarField& alpha1,
	const volScalarField& alpha2,
	const volScalarField& alpha3,
	const volVectorField& U,
    const IOdictionary& dict
//    const gtawIncompressibleThreePhaseMixture& mixture
)
:
    alpha1_(alpha1),
    alpha2_(alpha2),
    alpha3_(alpha3),
    U_(U),
    transportPropertiesDict_(dict),

    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(U_.mesh().V()), 1.0/3.0)
    ),

    sigma12_("sigma12", dimensionSet(1, 0, -2, 0, 0), transportPropertiesDict_.lookup("sigma12")),
    sigma13_("sigma13", dimensionSet(1, 0, -2, 0, 0), transportPropertiesDict_.lookup("sigma13")),
    sigma23_("sigma23", dimensionSet(1, 0, -2, 0, 0), transportPropertiesDict_.lookup("sigma23")),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    )
{
//    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceVectorField> 
Foam::gtawInterfaceProperties::nHatfv12() const
{
    surfaceVectorField gradAlpha12f
    (
        fvc::interpolate(alpha2_)*fvc::interpolate(fvc::grad(alpha1_))
      - fvc::interpolate(alpha1_)*fvc::interpolate(fvc::grad(alpha2_))
    );

    // Face unit interface normal
    return gradAlpha12f/(mag(gradAlpha12f) + deltaN_);
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::gtawInterfaceProperties::nHatfv13() const
{

    surfaceVectorField gradAlpha13f
    (
        fvc::interpolate(alpha3_)*fvc::interpolate(fvc::grad(alpha1_))
      - fvc::interpolate(alpha1_)*fvc::interpolate(fvc::grad(alpha3_))
    );

    // Face unit interface normal
    return gradAlpha13f/(mag(gradAlpha13f) + deltaN_);
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::gtawInterfaceProperties::nHatfv23() const
{

    surfaceVectorField gradAlpha23f
    (
        fvc::interpolate(alpha3_)*fvc::interpolate(fvc::grad(alpha2_))
      - fvc::interpolate(alpha2_)*fvc::interpolate(fvc::grad(alpha3_))
    );

    // Face unit interface normal
    return gradAlpha23f/(mag(gradAlpha23f) + deltaN_);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::gtawInterfaceProperties::nHatf12() const
{
    const fvMesh& mesh = alpha1_.mesh();
//    const fvMesh& mesh = mixture_.U().mesh();

    surfaceVectorField gradAlpha12f
    (
        fvc::interpolate(alpha2_)*fvc::interpolate(fvc::grad(alpha1_))
      - fvc::interpolate(alpha1_)*fvc::interpolate(fvc::grad(alpha2_))
    );

    return (nHatfv12() & mesh.Sf());
//    return (gradAlpha12f/(mag(gradAlpha12f) + deltaN_)) & mesh.Sf();
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::gtawInterfaceProperties::nHatf13(//) const
//            const volScalarField& alpha1,
//            const volScalarField& alpha3
        ) const
{
    const fvMesh& mesh = U_.mesh();
//    const fvMesh& mesh = mixture_.U().mesh();

    return (nHatfv13() & mesh.Sf());
//    return (nHatfv13(alpha1_, alpha3_) & mesh.Sf());
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::gtawInterfaceProperties::nHatf23(//) const
//            const volScalarField& alpha2,
//            const volScalarField& alpha3
        ) const
{
    const fvMesh& mesh = U_.mesh();
//    const fvMesh& mesh = mixture_.U().mesh();

    return (nHatfv23() & mesh.Sf());
//    return (nHatfv23(alpha2_, alpha3_) & mesh.Sf());
}

Foam::tmp<Foam::volScalarField> 
Foam::gtawInterfaceProperties::K12(//) const
//            const volScalarField& alpha1,
//            const volScalarField& alpha2
        ) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv12();//alpha1, alpha2);
//    correctContactAngle12(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());
    correctContactAngle12(tnHatfv.ref().boundaryFieldRef());
//    const fvMesh& mesh = mixture_.U().mesh();
    const fvMesh& mesh = U_.mesh();

    return -fvc::div(tnHatfv & mesh.Sf());
//    return -fvc::div(nHatf12());
}

Foam::tmp<Foam::volScalarField> 
Foam::gtawInterfaceProperties::K13(//) const
//            const volScalarField& alpha1,
//            const volScalarField& alpha3
        ) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv13();//alpha1, alpha3);
//    correctContactAngle13(alpha1, alpha3, tnHatfv.ref().boundaryFieldRef());
    correctContactAngle12(tnHatfv.ref().boundaryFieldRef());
//    const fvMesh& mesh = mixture_.U().mesh();
    const fvMesh& mesh = U_.mesh();
    
    return -fvc::div(tnHatfv & mesh.Sf());
//    return -fvc::div(nHatf13());
}

Foam::tmp<Foam::volScalarField> 
Foam::gtawInterfaceProperties::K23(//) const
//            const volScalarField& alpha2,
//            const volScalarField& alpha3
        ) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv23();//alpha2, alpha3);
//    correctContactAngle23(alpha2, alpha3, tnHatfv.ref().boundaryFieldRef());
    correctContactAngle12(tnHatfv.ref().boundaryFieldRef());
//    const fvMesh& mesh = mixture_.U().mesh();
    const fvMesh& mesh = U_.mesh();

    return -fvc::div(tnHatfv & mesh.Sf());
//    return -fvc::div(nHatf23());
}

Foam::tmp<Foam::surfaceScalarField>
Foam::gtawInterfaceProperties::surfaceTensionForce() const
{
    tmp<surfaceScalarField> stf//tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
    	        U_.time().timeName(),
	    		U_.mesh()
            ),
            U_.mesh(),
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    stf = (
            sigma12_*fvc::interpolate(K12())*(
//            sigma12_*fvc::interpolate(K12(alpha1_, alpha2_))*(
                                                fvc::interpolate(alpha2_)*fvc::snGrad(alpha1_)
                                              - fvc::interpolate(alpha1_)*fvc::snGrad(alpha2_))
/*
          + sigma13_*fvc::interpolate(K13())*(
//          + sigma13_*fvc::interpolate(K13(alpha1_, alpha3_))*(
                                                fvc::interpolate(alpha3_)*fvc::snGrad(alpha1_)
                                              - fvc::interpolate(alpha1_)*fvc::snGrad(alpha3_))
		  + sigma23_*fvc::interpolate(K23())*(
//          + sigma23_*fvc::interpolate(K23(alpha2_, alpha3_))*(
                                                fvc::interpolate(alpha3_)*fvc::snGrad(alpha2_)
                                              - fvc::interpolate(alpha2_)*fvc::snGrad(alpha3_))
*/
          );

    return stf;

}

Foam::tmp<Foam::volScalarField>
Foam::gtawInterfaceProperties::nearInterface() const
{
    return max(max
    (
	    pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_),
	    pos0(alpha2_ - 0.01)*pos0(0.99 - alpha2_)),
	    pos0(alpha3_ - 0.01)*pos0(0.99 - alpha3_)
    );
}

// ************************************************************************* //
