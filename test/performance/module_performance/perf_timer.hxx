// Copyright (C) 2014 INRIA
// Author(s): Nicolas Claude
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_PERF_TIMER
#define VERDANDI_PERF_TIMER

#include <iostream>
#include <map>

#include "perf_clock.hxx"

//! Some user-friendly macros
#define GET_TIMER(clock_name) Timer::GetClock(clock_name)
#define TIMER(clock_name) Timer::OnOff(clock_name)
#define RESET_TIMER(clock_name) Timer::ResetClock(clock_name)
#define RESET_ALL_TIMER Timer::ResetAll();
#define PRINT_TIMER(clock_name) Timer::PrintClock(clock_name)

using namespace std;


//! This class stores and manages a set of chronos.
/*! This class is used to manage a set of chronos and provides a few
  high-level macros which permit any test user to set marks in his code and to
  compute time lapses or to print them in a gtest-like syntax. */
class Timer
{
public:
    //! Static map of all the existing timers.
    static map<string, Clock> clock_map_;

public:
    static double GetClock(string);
    static void OnOff(string);
    static void ResetClock(string);
    static void ResetAll();
    static void PrintClock(string);
};


#endif
