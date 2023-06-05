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
* \brief Contains the classes required to extend the black-oil model with (simple) microbe interaction.
*/
#ifndef EWOMS_BLACK_OIL_MICROBES_MODULE_HH
#define EWOMS_BLACK_OIL_MICROBES_MODULE_HH

#include "blackoilproperties.hh"

#include <opm/models/blackoil/blackoilmicrobesparams.hh>
#include <opm/models/io/vtkblackoilmicrobesmodule.hh>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/BactPara.hpp>
#endif

#include <dune/common/fvector.hh>

#include <cmath>
#include <stdexcept>
#include <string>

namespace Opm {
/*!
* \ingroup BlackOil
* \brief Contains the high level supplements required to extend the black oil
*        model with (simple) microbe interaction.
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = MathToolbox<Evaluation>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { bacteriaConcentrationIdx = Indices::bacteriaConcentrationIdx };
    enum { contiBacteriaEqIdx = Indices::contiBacteriaEqIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    static constexpr unsigned enableMicrobes = enableMicrobesV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:

#if HAVE_ECL_INPUT
    /*!
    * \brief Initialize internal data needed for microbes module
    */
    static void initFromState(const EclipseState& eclState)
    {
        // Check if module is enabled but MICROBES is not present
        if (enableMicrobes && !eclState.runspec().microbes()) {
            throw std::runtime_error("Microbes module enabled at compile time, but the deck does not contain "
                                     "MICROBES!");
        }
        // Check for opposite of the above: module disabled but MICROBES is in deck
        else if (!enableMicrobes && eclState.runspec().microbes()) {
            throw std::runtime_error("Microbes module disabled at compile time, but deck contains MICROBES!");
        }
        
        // If MICROBES is not present, this module should not be active regardless
        if (!eclState.runspec().microbes())
            return;

        // Set parameters for bacteria source term, taken from the BACTPARA keyword
        const auto& BactPara = eclState.getBactPara();
            setBactPara(BactPara.getMaxGrowthRate(),
                        BactPara.getHalfVelocityGas(),
                        BactPara.getYieldCoefficient(),
                        BactPara.getStoichiometricCoefficient(),
                        BactPara.getDecayCoefficient());

    }
#endif

    /*!
    * \brief Specify the bacteria properties a single region.
    *
    * The index of specified here must be in range [0, numSatRegions)
    */
    static void setBactPara(const Scalar& maxGrowthRate,
                            const Scalar& halfVelocityGas,
                            const Scalar& yieldCoeff,
                            const Scalar& stoicCoeff,
                            const Scalar& decayCoeff)
    {
        params_.maxGrowthRate_ = maxGrowthRate;
        params_.halfVelocityGas_ = halfVelocityGas;
        params_.yieldCoeff_ = yieldCoeff;
        params_.stoicCoeff_ = stoicCoeff;
        params_.decayCoeff_ = decayCoeff;
    }

    /*!
    * \brief Register all run-time parameters for the black-oil microbes module.
    */
    static void registerParameters()
    {
        if (!enableMicrobes)
            return;

        VtkBlackOilMicrobesModule<TypeTag>::registerParameters();
    }

    /*!
    * \brief Register all microbes specific VTK and ECL output modules.
    */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableMicrobes)
            return;

        model.addOutputModule(new VtkBlackOilMicrobesModule<TypeTag>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableMicrobes)
            return false;
        
        // True if bacteria equation applies
        return eqIdx == contiBacteriaEqIdx;
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // Should be OK if bacteria concentration if in kg/m3 or similar.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableMicrobes)
            return;
        
        // Calculate volume water in SM3
        const auto& fs = intQuants.fluidState();
        LhsEval surfaceVolumeWater =
            Toolbox::template decay<LhsEval>(fs.saturation(oilPhaseIdx))
            * Toolbox::template decay<LhsEval>(fs.invB(oilPhaseIdx))
            * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // Avoid singular matrix if no water is present.
        surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);

        // Microbes suspended in water phase
        const LhsEval massMicrobes = 
            surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.bacteriaConcentration());
        storage[contiBacteriaEqIdx] += massMicrobes;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)
    {
        if (!enableMicrobes)
            return;
        
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(oilPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);


        if (upIdx == inIdx) {
            flux[contiBacteriaEqIdx] = extQuants.volumeFlux(oilPhaseIdx) * up.bacteriaConcentration();
        }
        else {
            flux[contiBacteriaEqIdx] = extQuants.volumeFlux(oilPhaseIdx) * decay<Scalar>(up.bacteriaConcentration());
        }
    }

    static void addSource(RateVector& source,
                            const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)
    {
        if (!enableMicrobes)
            return;
        
        // Get bacteria parameters
        Scalar mu = maxGrowthRate();
        Scalar alpha = halfVelocityGas();
        Scalar Y = yieldCoeff();
        Scalar gamma = stoicCoeff();
        Scalar kd = decayCoeff();

        // Convert RS to mole fraction for use in source term
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        const auto& Rs = fs.Rs();
        unsigned pvtRegionIndex = fs.pvtRegionIndex();

        const auto& xG = massFracToMoleFrac(pvtRegionIndex, RsToMassFraction(pvtRegionIndex, Rs));

        // Get saturation, porosity, bacteria concentration and inverse Bg for convenience
        const Evaluation& cBact = intQuants.bacteriaConcentration();
        const Evaluation& Sw = fs.saturation(oilPhaseIdx);
        const Evaluation& poro = intQuants.porosity();
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);
        // const Evaluation& invBo = fs.invB(oilPhaseIdx);

        // Calculate bacteria growth rate
        Evaluation kg = mu * (xG / (xG + alpha));

        // Compute source terms
        // Microbial growth and decay rate
        source[contiBacteriaEqIdx] += (kg - kd) * cBact;

        // Microbial consumption of dissolved gas is proportional to bacterial growth rate
        source[conti0EqIdx + gasPhaseIdx] += gamma * poro * Sw * cBact * kg / (Y * rho_gRef);
    }

    static const Scalar maxGrowthRate()
    {
        return params_.maxGrowthRate_;
    }

    static const Scalar halfVelocityGas()
    {
        return params_.halfVelocityGas_;
    }

    static const Scalar yieldCoeff()
    {
        return params_.yieldCoeff_;
    }

    static const Scalar stoicCoeff()
    {
        return params_.stoicCoeff_;
    }

    static const Scalar decayCoeff()
    {
        return params_.decayCoeff_;
    }

private:
    static BlackOilMicrobesParams<Scalar> params_;

    static Evaluation RsToMassFraction(unsigned regionIdx, const Evaluation& Rs) {
        Scalar rho_oRef = FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, regionIdx);
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);

        const Evaluation rho_oG = Rs * rho_gRef;

        return rho_oG/(rho_oRef + rho_oG);
    }

    static Evaluation massFracToMoleFrac(unsigned regionIdx, const Evaluation& XoG)
    {
        Scalar M_Gas = FluidSystem::molarMass(FluidSystem::gasCompIdx, regionIdx);
        Scalar M_Brine = FluidSystem::molarMass(FluidSystem::oilCompIdx, regionIdx);

        return XoG*M_Brine / (M_Gas*(1 - XoG) + XoG*M_Brine);
    }
};  // class BlackOilMicrobesModule

template <class TypeTag, bool enableMicrobesV>
BlackOilMicrobesParams< typename BlackOilMicrobesModule<TypeTag, enableMicrobesV>::Scalar >
BlackOilMicrobesModule<TypeTag, enableMicrobesV>::params_;

/*!
* \ingroup BlackOil
* \class Opm::BlackOilMicrobesIntensiveQuantities
*
* \brief Provides the volumetric quantities required for the equations needed by the
*        (simple) microbes extension of the black-oil model.
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using MicrobesModule = BlackOilMicrobesModule<TypeTag>;

    enum { bacteriaConcentrationIdx = Indices::bacteriaConcentrationIdx };

public:
    /*!
    * \brief Update intensive quantities associated with module for (simple) microbial interaction
    */

    void MicrobesPropertiesUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // update bacteria concentration from primary variables
        bacteriaConcentration_ = priVars.makeEvaluation(bacteriaConcentrationIdx, timeIdx, linearizationType);
    }

    const Evaluation& bacteriaConcentration() const
    { return bacteriaConcentration_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation bacteriaConcentration_;

};  // class BlackOilMicrobesIntensiveQuantities

template <class TypeTag>
class BlackOilMicrobesIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

public:
    void MicrobesPropertiesUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {}
    
    const Evaluation& bacteriaConcentration() const
    { throw std::logic_error("bacteriaConcentration() called but MICROBES is disabled!"); }

};  // class BlackOilMicrobesIntensiveQuantities<TypeTag, false>

/*!
* \ingroup BlackOil
* \class Opm::BlackOilMicrobesExtensiveQuantities
*
* \brief Extensive quantities for (simple) microbes interaction in black-oil model
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilMicrobesExtensiveQuantities<TypeTag, false>{};

}  // namespace Opm

#endif