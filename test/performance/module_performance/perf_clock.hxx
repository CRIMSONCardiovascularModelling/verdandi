// Copyright (C) 2010 INRIA
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


#ifndef VERDANDI_PERF_CLOCK
#define VERDANDI_PERF_CLOCK

#include <iostream>
#include <ctime>

using namespace std;

#ifdef VERDANDI_HAS_CXX11
#include <chrono>
using namespace chrono;
#endif


//! This class is a simple chronometer.
/*! This class is a simple chronometer using high-precision objects found in
  <chrono>. The aim of this class is to be as simple as possible, in order to
  reduce the time spent here and to improve the accuracy of the measure. If
  C++11 is not available, the accuracy of this class is reduced but can still
  be used, thanks to <ctime> library. The directive 'VERDANDI_HAS_CXX11'
  should be set for C++11 to be used. */
class Clock
{
private:
    //! Is the clock running?
    bool is_running_;

#ifdef VERDANDI_HAS_CXX11
    //! Point in time which helps to compute elapsed time.
    high_resolution_clock::time_point first_point_;
    //! Time lapse during which this chrono was running.
    duration<double> elapsed_time_;
#else
    // Point in time which helps to compute elapsed time (C++98).
    time_t first_point_;
    // Time lapse during which this chrono was running (C++98).
    time_t elapsed_time_;
#endif

public:
    Clock();
    ~Clock();
    void Start();
    void Stop();
    void OnOff();
    void Reset();
    double GetCurrentTime();
};


#endif
