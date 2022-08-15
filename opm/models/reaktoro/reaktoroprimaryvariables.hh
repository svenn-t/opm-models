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
* \copydoc Opm::ReaktoroPrimaryVariables
*/
#ifndef EWOMS_REAKTORO_PRIMARY_VARIABLES_HH
#define EWOMS_REAKTORO_PRIMARY_VARIABLES_HH

#include "reaktoroindices.hh"
#include "reaktoroproperties.hh"

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/common/energymodule.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <iostream>

namespace Opm {

/*!
* \ingroup ReaktoroModel
*
* \brief Represents the primary variables used in the multiphase chemical transport model using Reaktoro.
*
* This class is basically a Dune::FieldVector which can retrieve its
* contents from an aribitatry fluid state.
*/

template <class TypeTag>
class ReaktoroPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { x0Idx = Indices::x0Idx };

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    using EnergyModule = Opm::EnergyModule<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>()>;

    using Toolbox = Opm::MathToolbox<Evaluation>;

public:
    ReaktoroPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    /*!
    * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
    */
    ReaktoroPrimaryVariables(Scalar value) : ParentType(value)
    {
        Opm::Valgrind::CheckDefined(value);
        Opm::Valgrind::SetDefined(*this);
    }

    /*!
    * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
    * ImmisciblePrimaryVariables& )
    */
    ReaktoroPrimaryVariables(const ReaktoroPrimaryVariables& value) = default;
    ReaktoroPrimaryVariables& operator=(const ReaktoroPrimaryVariables& value) = default;

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams&,
                                /*isInEquilibrium*/bool = false)
    {
        // TODO: ensure chemical equilibrium before assigning primary variables
        assignNaive(fluidState);
    }

    /*!
    * \copydoc ImmisciblePrimaryVariables::assignNaive
    */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // reset everything
        (*this) = 0.0;

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // set mole fractions of first phase
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            (*this)[x0Idx + compIdx] = Toolbox::value(fluidState.moleFraction(0, compIdx));
            Opm::Valgrind::CheckDefined((*this)[x0Idx + compIdx]);
        }

        // set the pressure of the first phase
        (*this)[pressure0Idx] = Toolbox::value(fluidState.pressure(/*phaseIdx=*/0));
        Opm::Valgrind::CheckDefined((*this)[pressure0Idx]);
    }

}; // class ReaktoroPrimaryVariables

}  // namespace Opm

#endif
