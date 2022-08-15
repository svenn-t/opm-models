// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
This file is part of the Open Porous Media project (OPM).

OPM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

OPM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OPM.  If not, see <http://www.gnu.org/licenses/>.

Consult the COPYING file in the top-level source directory of this
module for the precise wording of the license and the list of
copyright holders.
*/
/*!
* \file
*
* \copydoc Opm::ReaktoroIntensiveQuantities
*/
#ifndef EWOMS_REAKTORO_INTENSIVE_QUANTITIES_HH
#define EWOMS_REAKTORO_INTENSIVE_QUANTITIES_HH

#include "reaktoroproperties.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/material/constraintsolvers/ReaktoroChemicalEquilibrium.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>

namespace Opm {
/*!
* \ingroup ReaktoroModel \ingroup IntensiveQuantities
*
* \brief Contains the quantities which are are constant within a finite volume in the chemical transport model using
* Reaktoro
*/
template <class TypeTag>
class ReaktoroIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { x0Idx = Indices::x0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { dimWorld = GridView::dimensionworld };

    using ReaktoroChemicalEquilibriumSolver = Opm::ReaktoroChemicalEquilibrium<Scalar, FluidSystem, Evaluation>;
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, /*storeEnthalpy=*/enableEnergy>;
    using DiffusionIntensiveQuantities = Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;

public:
    ReaktoroIntensiveQuantities()
    {}

    ReaktoroIntensiveQuantities(const ReaktoroIntensiveQuantities& other) = default;

    ReaktoroIntensiveQuantities& operator=(const ReaktoroIntensiveQuantities& other) = default;

    /*!
    * \brief IntensiveQuantities::update
    */
    void update(const ElementContext& elemCtx,
                unsigned dofIdx,
                unsigned timeIdx)
    {
        // parent type update
        ParentType::update(elemCtx, dofIdx, timeIdx);
        ParentType::checkDefined();

        //
        // Reaktoro chemical equilibrium solver
        //
        // Get primary variables
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // Set temperature
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        // Set pressure and mole fractions in fluid state
        Evaluation p = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, p);

        Evaluation lastX = 1.0;
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            const Evaluation& x = min(max(priVars.makeEvaluation(x0Idx + compIdx, timeIdx), 0.0), 1.0);
            fluidState_.setMoleFraction(0, compIdx, x);
            lastX -= x;
        }
        fluidState_.setMoleFraction(0, numComponents - 1, lastX);

        // output mole fractions before Reaktoro
        std::cout << "BEFORE Reaktoro : " << std::endl;
        std::cout << "x = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.moleFraction(0, compIdx) << ", ";
        }
        std::cout << std::endl;
        std::cout << "y = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.moleFraction(1, compIdx) << ", ";
        }
        std::cout << std::endl;

        // molarity
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updateAll(fluidState_);
        fluidState_.setDensity(0, FluidSystem::density(fluidState_, paramCache, 0));
        std::cout << "c_w = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.molarity(0, compIdx) << ", ";
        }
        std::cout << std::endl;
        std::cout << "c_g = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.molarity(1, compIdx) << ", ";
        }
        std::cout << std::endl;

        // compute moles in aquous phase (molarity * pore_volume)
        const auto& problem = elemCtx.problem();
        Scalar poreVolume = elemCtx.dofVolume(dofIdx, timeIdx) * problem.porosity(elemCtx, dofIdx, timeIdx);
        std::cout << "PV = " << poreVolume << std::endl;
        std::cout << "poro = " << problem.porosity(elemCtx, dofIdx, timeIdx) << std::endl;
        std::cout << "cell_vol = " << elemCtx.dofVolume(dofIdx, timeIdx) << std::endl;

        ComponentVector mole;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            mole[compIdx] = max(fluidState_.molarity(0, compIdx) * poreVolume, 0.0);
        }
        
        // print moles (in aqueous phase)
        std::cout << "mole = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << mole[compIdx] << ", ";
        }
        std::cout << std::endl;

        // print P and T
        std::cout << "P = " << fluidState_.pressure(0) << std::endl;
        std::cout << "T = " << fluidState_.temperature(0) << std::endl;

        // compute chemical equilibrium and update saturation, density, and mole fractions
        ReaktoroChemicalEquilibriumSolver::solve(fluidState_, paramCache, mole);

        // output mole fractions after Reaktoro
        std::cout << "AFTER Reaktoro : " << std::endl;
        std::cout << "x = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.moleFraction(0, compIdx) << ", ";
        }
        std::cout << std::endl;
        std::cout << "y = ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            std::cout << fluidState_.moleFraction(1, compIdx) << ", ";
        }
        std::cout << std::endl;

        // output density
        std::cout << "rho_w = " << fluidState_.density(0) << std::endl;
        std::cout << "rho_g = " << fluidState_.density(1) << std::endl;

        // output saturations
        std::cout << "Sw = " << fluidState_.saturation(0) << std::endl;
        std::cout << "Sg = " << fluidState_.saturation(1) << std::endl;

        // 
        // Update rest of the quantities
        //
        // rel. perm
        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        // viscosity
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Opm::Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

    }

    /*!
    * \copydoc ImmiscibleIntensiveQuantities::fluidState
    */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
    * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
    */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
    * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
    */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
    * \copydoc ImmiscibleIntensiveQuantities::mobility
    */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
    * \copydoc ImmiscibleIntensiveQuantities::porosity
    */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Evaluation porosity_;
    DimMatrix intrinsicPerm_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];

};  // class ReaktoroIntensiveQuantities

}  // namespace Opm

#endif