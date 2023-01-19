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
* \copydoc Opm::SulfateReducingBacteriaProblem
*/
#ifndef EWOMS_SULFATE_REDUCING_BACTERIA_PROBLEM_HH
#define EWOMS_SULFATE_REDUCING_BACTERIA_PROBLEM_HH

#include <opm/models/io/cubegridvanguard.hh>

#include <opm/material/fluidsystems/SRBFluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Opm {
template <class TypeTag>
class SulfateReducingBacteriaProblem;
}  // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct  SulfateReducingBacteriaBaseProblem {};
}  // namespace TTag

// Runtime options specific for this problem
template<class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct SimulationName { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct InitPressure { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct Injrate { using type = UndefinedProperty;};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { using type = Dune::YaspGrid<1>; };

// Grid vanguard
template<class TypeTag>
struct Vanguard<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { using type = Opm::CubeGridVanguard<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{ using type = Opm::SulfateReducingBacteriaProblem<TypeTag>; };

// Set fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Opm::SRBFluidSystem<Scalar>;
};

// Set fluid-matrix interactions
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static_assert(FluidSystem::numPhases == 2,
                  "A fluid system with two phases is required "
                  "for this problem!");

    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>;

public:
    using type = Opm::LinearMaterial<Traits>;
};

// Disable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { static constexpr bool value = false; };

// Default grid
template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { static constexpr int value = 250; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { static constexpr int value = 1; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { static constexpr int value = 1; };

// Default end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;  // seconds
};

// Default initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1;  // seconds
};

// Default runtime options specific for this problem
template<class TypeTag>
struct Temperature<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 293.15;  // K  (= 25 Celsius)
};

// Default initial pressure
template<class TypeTag>
struct InitPressure<TypeTag, TTag::SulfateReducingBacteriaBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;  // Pa  (= 100 bar)
};

// Default injection rate 
template <class TypeTag>
struct Injrate<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;  // kg / (m3 * s)
};

// Default simulation name
template<class TypeTag>
struct SimulationName<TypeTag, TTag::SulfateReducingBacteriaBaseProblem> { static constexpr auto value = "srbacteria"; };

}  // namespace Opm::Properties

namespace Opm {
/*!
* \ingroup TestProblems
* \brief Problem where hydrogen (H2) gas is injected in a aqueous solution with sulfate-reducing bacteria (with plenty
* of sulfate to feed on in addition to H2)
*
* The domain is 1D for simplicity. H2 is injected at left boundary though boundary conditions, and we have a free-flow
* boundary condition on the right.
*/
template <class TypeTag>
class SulfateReducingBacteriaProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    // Some type definitions
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    // Simulator and model definitions
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // Grid definitions
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    // Index shorthands for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { H2Idx = FluidSystem::H2Idx };
    enum { SO4Idx = FluidSystem::SO4Idx };
    enum { BactIdx = FluidSystem::BactIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    // Fluid-matrix interactions
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    // Vector and matrix definitions
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;

public:
    /*!
    * \copydoc Doxygen::defaultProblemConstructor
    */
    SulfateReducingBacteriaProblem(Simulator& simulator)
        : ParentType(simulator)
    { }
    
    /*!
    * \copydoc FvBaseProblem::finishInit
    */
    void finishInit()
    {
        // Base class init
        ParentType::finishInit();

        // Fluid system init
        FluidSystem::init();

        // Runtime parameters for this problem
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        initPressure_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitPressure);

        // Fluid-matrix parameters
        materialParams_.finalize();

        // Permeability
        K_ = this->toDimMatrix_(1.e-12);  // m2

        // Porosity
        porosity_ = 0.35;
    }
    
    /*!
    * \copydoc FvBaseMultiPhaseProblem::registerParameters
    */
    static void registerParameters()
    {
        // Register parent class parameters
        ParentType::registerParameters();

        // Register runtime options for this problem
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitPressure,
                             "The initial pressure [Pa] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Injrate, 
                            "Injection rate [kg/m3*s] of H2 into the reservoir");
    }

    /*!
    * \copydoc FvBaseProblem::name
    */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName);
        return oss.str();
    }
    
    /*!
    * \copydoc FvBaseProblem::name
    */
    void endTimeStep()
    {
#ifndef NDEBUG
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
#endif // NDEBUG
    }

    /*!
    * \copydoc FvBaseMultiPhaseProblem::temperature
    */
    template <class Context>
    Scalar temperature(const Context& /*context*/,
                       unsigned /*spaceIdx*/,
                       unsigned /*timeIdx*/) const
    { return temperature_; }

    /*!
    * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
    */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& /*context*/,
                                           unsigned /*spaceIdx*/,
                                           unsigned /*timeIdx*/) const
    { return K_; }

    /*!
    * \copydoc FvBaseMultiPhaseProblem::porosity
    */
    template <class Context>
    Scalar porosity(const Context& /*context*/,
                    unsigned /*spaceIdx*/,
                    unsigned /*timeIdx*/) const
    { return porosity_; }

    /*!
    * \copydoc FvBaseMultiPhaseProblem::materialLawParams
    */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& /*context*/,
                      unsigned /*spaceIdx*/,
                      unsigned /*timeIdx*/) const
    { return materialParams_; }

    /*!
    * \copydoc FvBaseProblem::boundary
    */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        // Get boundary coordinates
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        
        // Injection of H2 at left boundary
        if (pos[0] < 1e-6) {
            Scalar injrate = EWOMS_GET_PARAM(TypeTag, Scalar, Injrate);
            RateVector massRate(0.0);
            massRate[conti0EqIdx + H2Idx] = injrate;
            values.setMassRate(massRate);
        }

        // Dirichlet at right boundary
        else if (pos[0] > this->boundingBoxMax()[0] - 1e-6) {
            Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }

        // No-flow on the other boundaries
        else
            values.setNoFlow();
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        // Initialize fluid state
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // Assign fluid state
        values.assignNaive(fs);
    }

    /*!
    * \copydoc FvBaseProblem::source
    *
    * Chemical reactions and microbial activity
    */
    template <class Context>
    void source(RateVector& rate,
                const Context& context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        // Get up-to-date intensive quantities
        const IntensiveQuantities& intQuants = context.intensiveQuantities(spaceIdx, timeIdx);

        // Microbial growth rate is a double Monod model involving H2 and sulfate
        Evaluation xH2 = intQuants.fluidState().moleFraction(liquidPhaseIdx, H2Idx);
        Evaluation xSO4 = intQuants.fluidState().moleFraction(liquidPhaseIdx, SO4Idx);
        Scalar kg_max = 2.246e-5;
        Scalar kH2 = 5.018e-8;
        Scalar kSO4 = 1.8e-6;

        Evaluation kg = kg_max * (xH2 / (kH2 + xH2)) * (xSO4/ (xSO4 + kSO4));

        // Microbial decay (constant)
        Evaluation kd = 5e-7;

        // Insert rates for brine, H2, sulfate, and bacteria
        Evaluation rhoMolar = intQuants.fluidState().molarDensity(liquidPhaseIdx);
        Evaluation poro = intQuants.porosity();
        Evaluation Sw = intQuants.fluidState().saturation(liquidPhaseIdx);
        Evaluation xBact = intQuants.fluidState().moleFraction(liquidPhaseIdx, BactIdx);
        Scalar yield = 8.639e12;
        Evaluation rg = poro * Sw * rhoMolar * xBact * kg / yield;
        rate[conti0EqIdx + BrineIdx] = 4 * rg;
        rate[conti0EqIdx + H2Idx] = -5 * rg;
        rate[conti0EqIdx + SO4Idx] = -1 * rg;
        rate[conti0EqIdx + BactIdx] = (Sw * kg - kd) / rhoMolar;
    }

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        // Set temperature and initial pressure
        fs.setPressure(liquidPhaseIdx, initPressure_ * 1e5);  // [Pa]
        fs.setPressure(gasPhaseIdx, initPressure_ * 1e5);  // [Pa]

        fs.setTemperature(temperature_);

        // Set saturations
        fs.setSaturation(liquidPhaseIdx, 1.0);
        fs.setSaturation(gasPhaseIdx, 0.0);

        // Set liquid mole fractions
        Scalar xH2 = 0.0;
        Scalar xSO4 = 0.05;
        Scalar xBact = 0.01;
        Scalar xBrine = 1 - xH2 - xSO4 - xBact;
        fs.setMoleFraction(liquidPhaseIdx, H2Idx, xH2);
        fs.setMoleFraction(liquidPhaseIdx, SO4Idx, xSO4);
        fs.setMoleFraction(liquidPhaseIdx, BactIdx, xBact);
        fs.setMoleFraction(liquidPhaseIdx, BrineIdx, xBrine);

        // Calculate the rest from liquid phase info
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/false);
    }

    // Define variables
    Scalar temperature_;
    Scalar initPressure_;
    Scalar porosity_;
    DimMatrix K_;
    MaterialLawParams materialParams_;
};  // class SulfateReducingBacteriaProblem

}  // namespace Opm

#endif