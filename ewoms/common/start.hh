/*
  Copyright (C) 2010-2013 by Andreas Lauser
  Copyright (C) 2011 by Philipp Nuske
  Copyright (C) 2012 by Bernd Flemisch

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
*/
/*!
 * \file
 * \brief Provides convenience routines to bring up the simulation at runtime.
 */
#ifndef EWOMS_START_HH
#define EWOMS_START_HH

#include <opm/core/utility/PropertySystem.hpp>
#include "parametersystem.hh"
#include <opm/material/Valgrind.hpp>

#include <ewoms/version.hh>
#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/simulator.hh>
#include <ewoms/common/timer.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/common/parametertreeparser.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <locale>

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <termios.h>
#include <signal.h>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Opm {
// forward declaration of property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(ThreadManager);
NEW_PROP_TAG(PrintProperties);
NEW_PROP_TAG(PrintParameters);
NEW_PROP_TAG(ParameterFile);
} // namespace Properties
} // namespace Opm
//! \cond SKIP_THIS

namespace Ewoms {
/*!
 * \brief Register all runtime parameters, parse the command line
 *        arguments and the parameter file.
 *
 * \param argc The number of command line arguments
 * \param argv Array with the command line argument strings
 */
template <class TypeTag>
int setupParameters_(int argc, char **argv)
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParameterMetaData;

    // first, get the MPI rank of the current process
    int myRank = 0;
#if HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

    ////////////////////////////////////////////////////////////
    // Register all parameters
    ////////////////////////////////////////////////////////////
    EWOMS_REGISTER_PARAM(TypeTag, std::string, ParameterFile,
                         "An .ini file which contains a set of run-time "
                         "parameters");
    EWOMS_REGISTER_PARAM(TypeTag, int, PrintProperties,
                         "Print the values of the compile time properties at "
                         "the start of the simulation");
    EWOMS_REGISTER_PARAM(TypeTag, int, PrintParameters,
                         "Print the values of the run-time parameters at the "
                         "start of the simulation");

    Simulator::registerParameters();
    ThreadManager::registerParameters();

    ////////////////////////////////////////////////////////////
    // set the parameter values
    ////////////////////////////////////////////////////////////

    // check whether the user wanted to see the help message
    for (int i = 1; i < argc; ++i) {
        if (std::string("--help") == argv[i] || std::string("-h") == argv[i]) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0], "", /*handleHelp=*/true, std::cout);
            return /*status=*/2;
        }
    }

    // fill the parameter tree with the options from the command line
    std::string s = Parameters::parseCommandLineOptions<TypeTag>(argc, argv);
    if (!s.empty()) {
        return /*status=*/1;
    }

    std::string paramFileName = EWOMS_GET_PARAM_(TypeTag, std::string, ParameterFile);
    if (paramFileName != "") {
        ////////////////////////////////////////////////////////////
        // add the parameters specified using an .ini file
        ////////////////////////////////////////////////////////////

        // check whether the parameter file is readable.
        std::ifstream tmp;
        tmp.open(paramFileName.c_str());
        if (!tmp.is_open()) {
            std::ostringstream oss;
            if (myRank == 0) {
                oss << "Parameter file \"" << paramFileName
                    << "\" does not exist or is not readable.";
                Parameters::printUsage<TypeTag>(argv[0], oss.str());
            }
            return /*status=*/1;
        }

        // read the parameter file.
        Dune::ParameterTreeParser::readINITree(paramFileName,
                                               ParameterMetaData::tree(),
                                               /*overwrite=*/false);
    }

    EWOMS_END_PARAM_REGISTRATION(TypeTag);

    return /*status=*/0;
}

/*!
 * \brief Resets the current TTY to a usable state if the program was interrupted by
 *        SIGABRT or SIGINT.
 */
static void resetTerminal_(int signum)
{
    // first thing to do when a nuke hits: restore the default signal handler
    signal(signum, SIG_DFL);

    // the following code resets the terminal status and is loosely based on corutils'
    // "stty" utility. The copyright notice for this file is the following:
    //
    // Copyright (C) 1990-2014 Free Software Foundation, Inc.
    //
    // This program is free software: you can redistribute it and/or modify it under the
    // terms of the GNU General Public License as published by the Free Software
    // Foundation, either version 3 of the License, or (at your option) any later
    // version.
    struct termios mode;
    tcgetattr(STDERR_FILENO, &mode);

    mode.c_cc[VINTR] = CINTR;
    mode.c_cc[VQUIT] = CQUIT;
    mode.c_cc[VERASE] = CERASE;
    mode.c_cc[VKILL] = CKILL;
    mode.c_cc[VEOF] = CEOF;
    mode.c_cc[VEOL] = CEOL;
    mode.c_cc[VEOL2] = _POSIX_VDISABLE;
    mode.c_cc[VSWTC] = CSUSP;
    mode.c_cc[VSTART] = CSTART;
    mode.c_cc[VSTOP] = CSTOP;
    mode.c_cc[VSUSP] = CSUSP;
    mode.c_cc[VREPRINT] = CRPRNT;
    mode.c_cc[VWERASE] = CWERASE;
    mode.c_cc[VLNEXT] = CLNEXT;
    mode.c_cc[VDISCARD] = 0x1f & 'o';
    mode.c_cc[VMIN] = 1;
    mode.c_cc[VTIME] = 0;

    // control flags
    mode.c_cflag |= CREAD;

    // input flags
    mode.c_iflag &= ~IGNBRK;
    mode.c_iflag |= BRKINT;
    mode.c_iflag &= ~INLCR;
    mode.c_iflag &= ~IGNCR;
    mode.c_iflag |= ICRNL;
    mode.c_iflag &= ~IXOFF;
    mode.c_iflag &= ~IUCLC;
    mode.c_iflag &= ~IXANY;
    mode.c_iflag |= IMAXBEL;
    mode.c_iflag &= ~IUTF8;

    // output flags
    mode.c_oflag |= OPOST;
    mode.c_oflag &= ~OLCUC;
    mode.c_oflag &= ~OCRNL;
    mode.c_oflag |= ONLCR;
    mode.c_oflag &= ~ONOCR;
    mode.c_oflag &= ~ONLRET;
    mode.c_oflag &= ~OFILL;
    mode.c_oflag &= ~OFDEL;
    mode.c_oflag &= ~NL1;
    mode.c_oflag |= NL0;
    mode.c_oflag &= ~CR3;
    mode.c_oflag &= ~CR2;
    mode.c_oflag &= ~CR1;
    mode.c_oflag |= CR0;
    mode.c_oflag &= ~TAB3;
    mode.c_oflag &= ~TAB2;
    mode.c_oflag &= ~TAB1;
    mode.c_oflag |= TAB0;
    mode.c_oflag &= ~BS1;
    mode.c_oflag |= BS0;
    mode.c_oflag &= ~VT1;
    mode.c_oflag |= VT0;
    mode.c_oflag &= ~FF1;
    mode.c_oflag |= FF0;

    // local flags
    mode.c_lflag |= ISIG;
    mode.c_lflag |= ICANON;
    mode.c_lflag |= IEXTEN;
    mode.c_lflag |= ECHO;
    mode.c_lflag |= ECHOE;
    mode.c_lflag |= ECHOK;
    mode.c_lflag &= ~ECHONL;
    mode.c_lflag &= ~NOFLSH;
    mode.c_lflag &= ~XCASE;
    mode.c_lflag &= ~TOSTOP;
    mode.c_lflag &= ~ECHOPRT;
    mode.c_lflag |= ECHOCTL;
    mode.c_lflag |= ECHOKE;

    // set the control mode of the TTY
    tcsetattr(STDIN_FILENO, TCSADRAIN, &mode);

    const char resetString[] = {
        // switch back to the default character set
        0x0f,

        // reset current attributes
        27, '[', 'm',

        // disable line wrapping
        27, '[', '7', 'l',

        // show text cursor
        27, '[', '2', '5', 'h',

        0
    };

    std::cerr << resetString << std::flush;

    // sleep a bit to give the terminal some time to contemplate about the commands...
    struct timespec sleepTime;
    sleepTime.tv_sec = 0;
    sleepTime.tv_nsec = 100 * 1000 * 1000;
    nanosleep(&sleepTime, NULL);

    // print a new line to decrease the possibility of garbage remaining on the line
    // which shows the command line prompt
    std::cerr << "\r\n" << std::flush;

    // after we did our best to clean the pedestrian way, re-raise the signal
    raise(signum);
}
//! \endcond

/*!
 * \ingroup Startup
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments
 * \param argv The array of the command line arguments
 */
template <class TypeTag>
int start(int argc, char **argv)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;

    // set the signal handlers to reset the TTY to a well defined state on unexpected
    // program aborts
    if (isatty(STDIN_FILENO)) {
        signal(SIGINT, resetTerminal_);
        signal(SIGHUP, resetTerminal_);
        signal(SIGABRT, resetTerminal_);
        signal(SIGFPE, resetTerminal_);
        signal(SIGSEGV, resetTerminal_);
        signal(SIGPIPE, resetTerminal_);
        signal(SIGTERM, resetTerminal_);
    }

    // initialize MPI, finalize is done automatically on exit
    const Ewoms::MpiHelper mpiHelper(argc, argv);

    int myRank = mpiHelper.rank();

    try
    {
        int paramStatus = setupParameters_<TypeTag>(argc, argv);
        if (paramStatus == 1)
            return 1;
        if (paramStatus == 2)
            return 0;

        ThreadManager::init();

        // read the initial time step and the end time
        double endTime;
        double initialTimeStepSize;

        endTime = EWOMS_GET_PARAM(TypeTag, Scalar, EndTime);
        if (endTime < -1e50) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0],
                                                "Mandatory parameter '--end-time' not specified!");
            return 1;
        }

        initialTimeStepSize = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        if (initialTimeStepSize < -1e50) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0],
                                                "Mandatory parameter '--initial-time-step-size' "
                                                "not specified!");
            return 1;
        }


        if (myRank == 0)
            std::cout << "eWoms " << EWOMS_VERSION
#ifdef EWOMS_CODENAME
                      << " (\"" << EWOMS_CODENAME << "\")"
#endif
                      << " will now start the trip. "
                      << "Please sit back, relax and enjoy the ride.\n"
                      << std::flush;

        // print the parameters if requested
        int printParams = EWOMS_GET_PARAM(TypeTag, int, PrintParameters);
        if (myRank == 0) {
            std::string endParametersSeparator("# [end of parameters]\n");
            if (printParams) {
                bool printSeparator = false;
                if (printParams == 1 || !isatty(fileno(stdout))) {
                    Ewoms::Parameters::printValues<TypeTag>();
                    printSeparator = true;
                }
                else
                    // always print the list of specified but unused parameters
                    printSeparator =
                        printSeparator ||
                        Ewoms::Parameters::printUnused<TypeTag>();
                if (printSeparator)
                    std::cout << endParametersSeparator;
            }
            else
                // always print the list of specified but unused parameters
                if (Ewoms::Parameters::printUnused<TypeTag>())
                    std::cout << endParametersSeparator;
        }

        // print the properties if requested
        int printProps = EWOMS_GET_PARAM(TypeTag, int, PrintProperties);
        if (printProps && myRank == 0) {
            if (printProps == 1 || !isatty(fileno(stdout)))
                Opm::Properties::printValues<TypeTag>();
        }

        // instantiate and run the concrete problem. make sure to
        // deallocate the problem and before the time manager and the
        // grid
        Simulator simulator;
        simulator.run();

        if (myRank == 0) {
            std::cout << "eWoms reached the destination. If it is not the one that was intended, "
                      << "change the booking and try again.\n"
                      << std::flush;
        }
        return 0;
    }
    catch (std::exception &e)
    {
        if (myRank == 0)
            std::cout << e.what() << ". Abort!\n" << std::flush;
        return 1;
    }
    catch (Dune::Exception &e)
    {
        if (myRank == 0)
            std::cout << "Dune reported an error: " << e.what() << std::endl  << std::flush;
        return 2;
    }
    catch (...)
    {
        if (myRank == 0)
            std::cout << "Unknown exception thrown!\n" << std::flush;
        return 3;
    }
}

} // namespace Ewoms

#endif
