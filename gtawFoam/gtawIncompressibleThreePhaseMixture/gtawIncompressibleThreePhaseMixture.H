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
    Foam::gtawIncompressibleThreePhaseMixture

Description

SourceFiles
    gtawIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef gtawIncompressibleThreePhaseMixture_H
#define gtawIncompressibleThreePhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "threePhaseMixture/threePhaseMixture.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "basicThermo.H"
//#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gtawIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class gtawIncompressibleThreePhaseMixture
:
    public IOdictionary,
    public transportModel,
    public threePhaseMixture
{
private:

// 24601 FOR RADIATION
/*
        autoPtr<basicThermo> thermo1_;
        autoPtr<basicThermo> thermo2_;
        autoPtr<basicThermo> thermo3_;
*/
/*
        autoPtr<rhoThermo> thermo1_;
        autoPtr<rhoThermo> thermo2_;
        autoPtr<rhoThermo> thermo3_;
*/
//\24601
protected:
        scalar a = 1.0;
        scalar b = 1.0;
        scalar c = 1.0;
        scalar d = 1.0;
        scalar e = 1.0;


        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
        autoPtr<viscosityModel> nuModel3_;

        word materialNameLiquid_;
        word materialNameSolid_;        

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        bool rho2Opt_;
        dimensionedScalar rho2a_;
        dimensionedScalar rho2b_;
        dimensionedScalar rho2c_;
        dimensionedScalar rho2d_;

        dimensionedScalar rho3_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dimensionedScalar Tref_;
        dimensionedScalar beta_;

        dimensionedScalar CTLE_;

        bool Tdep1_;
        bool Tdep2_;
        bool Tdep3_;
        bool Shomate_;

        dimensionedScalar cp1_;

        dimensionedScalar cp2_;
        dimensionedScalar cp2a_;
        dimensionedScalar cp2b_;
        dimensionedScalar cp2c_;
		dimensionedScalar cp2d_;
        dimensionedScalar cp2e_;
        dimensionedScalar cp2f_;        
        

        dimensionedScalar cp3_;
        dimensionedScalar cp3a_;
        dimensionedScalar cp3b_;
        dimensionedScalar cp3c_;
        dimensionedScalar cp3d_;
        dimensionedScalar cp3e_;
        dimensionedScalar cp3f_;

        dimensionedScalar smallT_;
        dimensionedScalar smallNu_;

        dimensionedScalar k1_;

        dimensionedScalar k2_;
        bool k2Opt_;
        dimensionedScalar k2a_;
        dimensionedScalar k2b_;
        dimensionedScalar k2c_;                

        dimensionedScalar k3_;
        bool k3Opt_;
        dimensionedScalar k3a_;
        dimensionedScalar k3b_;
        dimensionedScalar k3c_;                

        volScalarField nu_;
        volScalarField k_;
        volScalarField cp_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();

        void calcK();

        void updateCpMembers();

public:
   
    // Constructors

        //- Construct from components
        gtawIncompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi//,
        );


    //- Destructor
    virtual ~gtawIncompressibleThreePhaseMixture()
    {}


    // Member Functions

		void calcNuUpdate()
		{
			calcNu();
		}	

// 24601 FOR RADIATION
/*
        //- Return const-access to phase1 viscosityModel
//        const rhoThermo& thermo1() const
        const basicThermo& thermo1() const
        {
            return thermo1_();
        }

        //- Return const-access to phase1 viscosityModel
//        const rhoThermo& thermo2() const
        const basicThermo& thermo2() const
        {
            return thermo2_();
        }

        //- Return const-access to phase1 viscosityModel
//        const rhoThermo& thermo3() const
        const basicThermo& thermo3() const
        {
            return thermo3_();
        }
*/
//\24601

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase3 viscosityModel
        const viscosityModel& nuModel3() const
        {
            return nuModel3_();
        }

        const word& materialNameLiquid() const
        {
            return materialNameLiquid_;
        }

        const word& materialNameSolid() const
        {
            return materialNameSolid_;
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase3 density
        const dimensionedScalar& rho3() const
        {
            return rho3_;
        };

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& cp2() const
        {
            return cp2_;
        }

        //- Return const-access to phase3 density
        const dimensionedScalar& cp3() const
        {
            return cp3_;
        }

        const dimensionedScalar& k1() const
        {
            return k1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& k2() const
        {
            return k2_;
        }

        //- Return const-access to phase3 density
        const dimensionedScalar& k3() const
        {
            return k3_;
        }

        const volScalarField& cp2V() const
        {
            return cp2a_ + cp2b_*T_;
        }      
        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;
/*
        tmp<volScalarField> nuAcess() const;
        {
            return nu_;
        }
*/
        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
		virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        tmp<volScalarField> k() const;
        tmp<volScalarField> k2T() const;
        tmp<volScalarField> k3T() const;        

        tmp<surfaceScalarField> kf() const;

        tmp<volScalarField> cp() const;

        tmp<surfaceScalarField> cpf() const;

        tmp<volScalarField> alpha() const;

        tmp<volScalarField> rhoExp() const;
        tmp<volScalarField> rhoBouss() const;
        tmp<volScalarField> rhoBoussRho() const;
        
		tmp<volScalarField> cp2TDepVSF() const;
        tmp<volScalarField> cp3TDepVSF() const;

        tmp<surfaceScalarField> alphaf() const;
        
        tmp<volScalarField> rho() const;
        tmp<volScalarField> rhoCp() const;

        void callUpdateCpMembers()
        {
            updateCpMembers();
        }

        //- Correct the laminar viscosity
        virtual void correct() // corrects ALL fields
        {
            void calcNu();
            void calcK();
            void updateCpMembers();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
// 24601 FOR RADIATION
/*
        void correctThermo();
*/
//\24601
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //