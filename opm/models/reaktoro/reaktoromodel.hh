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
* \copydoc Opm::ReaktoroModel
*/
#ifndef EWOMS_REAKTORO_MODEL_HH
#define EWOMS_REAKTORO_MODEL_HH

#include "reaktoroproperties.hh"
#include "reaktorolocalresidual.hh"
#include "reaktoroprimaryvariables.hh"
#include "reaktororatevector.hh"
#include "reaktoroboundaryratevector.hh"
#include "reaktorointensivequantities.hh"
#include "reaktoroextensivequantities.hh"
#include "reaktoroindices.hh"

#include <opm/models/common/multiphasebasemodel.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <sstream>
#include <string>

namespace Opm {

template <class TypeTag>
class ReaktoroModel;

}  // namespace Opm

namespace Opm::Properties {

namespace TTag {
/*!
* \brief Define the type tag for the Reaktoro model.
*/
struct ReaktoroModel { using InheritsFrom = std::tuple<VtkDiffusion,
                                                       VtkEnergy,
                                                       VtkComposition,
                                                       MultiPhaseBaseModel>; 
};  // struct ReaktoroModel
}  // namespace TTag

//! Use the local Jacobian operator for the Reaktoro model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroLocalResidual<TypeTag>; };

//! the Model property
template<class TypeTag>
struct Model<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroModel<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroPrimaryVariables<TypeTag>; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroBoundaryRateVector<TypeTag>; };

//! the IntensiveQuantities property
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroIntensiveQuantities<TypeTag>; };

//! the ExtensiveQuantities property
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroExtensiveQuantities<TypeTag>; };

//! The indices required by the Reaktoro model
template<class TypeTag>
struct Indices<TypeTag, TTag::ReaktoroModel> { using type = Opm::ReaktoroIndices<TypeTag, /*PVIdx=*/0>; };

//! Disable the energy equation by default
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::ReaktoroModel> { static constexpr bool value = false; };

// disable molecular diffusion by default
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::ReaktoroModel> { static constexpr bool value = false; };

//! The basis value for the weight of the pressure primary variable
template<class TypeTag>
struct ReaktoroPressureBaseWeight<TypeTag, TTag::ReaktoroModel>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

//! The basis value for the weight of the mole fraction primary variables
template<class TypeTag>
struct ReaktoroMoleFractionsBaseWeight<TypeTag, TTag::ReaktoroModel>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

}  // namespace Opm::Properties

namespace Opm {

/*!
* \ingroup ReaktoroModel
* 
* \brief A multiphase chemical transport model using Reaktoro
*/
template <class TypeTag>
class ReaktoroModel
    : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;

public:
    ReaktoroModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
    * \brief Register all run-time parameters for the chemical transport model using Reaktoro.
    */
    static void registerParameters()
    {
        // generic parameters
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        // compositional parameters
        Opm::VtkCompositionModule<TypeTag>::registerParameters();

        // diffusion parameters
        if (enableDiffusion)
            Opm::VtkDiffusionModule<TypeTag>::registerParameters();

        // energy parameters
        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();
            
    }

    /*!
    * \copydoc FvBaseDiscretization::name
    */
    static std::string name()
    { return "reaktoro"; }

    /*!
    * \copydoc FvBaseDiscretization::primaryVarName
    */
    std::string primaryVarName(unsigned pvIdx) const
    {
        // check energy module
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;
        
        // pressure and mole fraction names
        std::ostringstream oss;
        if (pvIdx == Indices::pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (Indices::x0Idx <= pvIdx && pvIdx < Indices::x0Idx + numComponents - 1)
            oss << "moleFrac^" << FluidSystem::componentName(pvIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
    * \copydoc FvBaseDiscretization::eqName
    */
    std::string eqName(unsigned eqIdx) const
    {
        // check energy module
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;
        
        // continuity equations
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents) {
            unsigned compIdx = eqIdx - Indices::conti0EqIdx;
            oss << "continuity^" << FluidSystem::componentName(compIdx);
        }
        else
            assert(false);

        return oss.str();

    }

    /*!
    * \copydoc FvBaseDiscretization::primaryVarWeight
    */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        // check energy module
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        if (tmp > 0)
            return tmp;
        
        // pressure
        if (Indices::pressure0Idx == pvIdx)
            return getPropValue<TypeTag, Properties::ReaktoroPressureBaseWeight>();
        
        // mole fraction
        else
            return getPropValue<TypeTag, Properties::ReaktoroMoleFractionsBaseWeight>();
    
    }

    /*!
    * \copydoc FvBaseDiscretization::eqWeight
    */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        // check energy module
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            return tmp;

        // for continuity equations, convert from mole to mass [kg]
        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(compIdx <= numComponents);
        return FluidSystem::molarMass(compIdx);
    }

    /*!
    * \internal
    */
    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        this->addOutputModule(new Opm::VtkCompositionModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Opm::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }

};  // class ReaktoroModel

}  // namespace Opm

#endif