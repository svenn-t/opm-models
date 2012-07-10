// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

/*!
 * \file
 *
 * \brief Defines the indices required for the compositional flow model based 
 *        on flash calculations.
 */
#ifndef DUMUX_FLASH_INDICES_HH
#define DUMUX_FLASH_INDICES_HH

namespace Dumux {
/*!
 * \brief Defines the indices required for the compositional flow model based 
 *        on flash calculations.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset = 0>
class FlashIndices
{
public:
    // Primary variable indices
    static const int cTot0Idx = PVOffset; //!< Index of the total concentration of the first component in the pore space.

    // equation indices
    static const int conti0EqIdx = PVOffset; //!< Index of the mass conservation equation for the first component.
};

}

#endif