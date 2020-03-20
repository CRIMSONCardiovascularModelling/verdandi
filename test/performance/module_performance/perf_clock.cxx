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


#include "perf_clock.hxx"


#ifdef VERDANDI_HAS_CXX11


//////////////
// C++ 2011 //
//////////////


//! Constructor.
Clock::Clock()
{
    // Creates a clock, initially at 0 s.
    bool is_running_ = false;
    first_point_ = high_resolution_clock::now() ;
    elapsed_time_
        = duration_cast<duration<double>>(first_point_ - first_point_);
}


//! Destructor.
Clock::~Clock()
{
}


//! Starts the clock, setting a new initial point in time.
void Clock::Start()
{
    if (is_running_)
        cout << "Clock already running" << endl;
    else
    {
        is_running_ = true;
        first_point_ = high_resolution_clock::now();
    }
}


//! Stops the clock and updates the elapsed time.
void Clock::Stop()
{
    duration<double> time_span;
    if (is_running_)
    {
        // Computes the elapsed time since the start.
        high_resolution_clock::time_point t2  = high_resolution_clock::now();
        time_span = duration_cast<duration<double>>(t2 - first_point_);
        // Computes the new total time.
        elapsed_time_ += time_span;
        is_running_ = false;
    }
    else
        cout << "Clock already stopped" << endl;
}


//! Returns a double corresponding to the elapsed time.
/*!
  \return The elapsed time.
*/
double Clock::GetCurrentTime()
{
    duration<double> time_span;
    if (is_running_)
    {
        high_resolution_clock::time_point t2  = high_resolution_clock::now();
        time_span = duration_cast<duration<double>>(t2 - first_point_);
        return elapsed_time_.count() + time_span.count();
    }
    else
        return elapsed_time_.count();
}


// Turns the clock on or off.
void Clock::OnOff()
{
    if(is_running_)
        this->Stop();
    else
        this->Start();
}


// Resets the clock.
void Clock::Reset()
{
    elapsed_time_
        = duration_cast<duration<double>>(first_point_ - first_point_);
    is_running_ = false;
}


#else


////////////
// C++ 98 //
////////////


//! Constructor.
Clock::Clock()
{
    // Creates a clock, initally stopped at 0s.
    bool is_running_ = false;
    time(&first_point_);
    elapsed_time_ = 0;
}


//! Destructor.
Clock::~Clock()
{
}


//! Starts the clock, setting a new initial point in time.
void Clock::Start()
{
    if (is_running_)
        cout << "Clock already running" << endl;
    else
    {
    is_running_ = true;
    time(&first_point_);
    }
}


//! Stops the clock and updates the elapsed time.
void Clock::Stop()
{
    time_t time_span;
    if (is_running_)
    {
        // Computes the elapsed time since the start.
        time_t t2  = time(NULL);
        time_span = t2 - first_point_;
        // Computes the new total time.
        elapsed_time_ += time_span;
        is_running_ = false;
    }
    else
        cout << "Clock already stopped" << endl;
}


//! Returns a double corresponding to the elapsed time.
/*!
  \return The elapsed time.
*/
double Clock::GetCurrentTime()
{
    time_t time_span;
    if (is_running_)
    {
    time_t t2  = time(NULL);
        time_span = t2 - first_point_;
        return elapsed_time_ + time_span;
    }
    else
        return elapsed_time_;
}


//! Turns the clock on or off.
void Clock::OnOff()
{
    if(is_running_)
        this->Stop();
    else
        this->Start();
}


//! Resets the clock.
void Clock::Reset()
{
    elapsed_time_ = 0;
    is_running_ = false;
}


#endif
