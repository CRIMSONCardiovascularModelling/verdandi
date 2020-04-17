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


#include "perf_timer.hxx"
#include "perf_clock.cxx"


map<string, Clock> Timer::clock_map_;


// Switches the clock  between On and Off.
void Timer::OnOff(string label)
{
    if (Timer::clock_map_.count(label) > 0)
        Timer::clock_map_.at(label).OnOff();
    else
    {
        Clock newClock;
        Timer::clock_map_.insert(pair<string,Clock>(label, newClock));
    }
}


// Returns the actual value of the clock.
double Timer::GetClock(string label)
{
    if (Timer::clock_map_.count(label) > 0)
        return Timer::clock_map_.at(label).GetCurrentTime();
    else
    {
        cout << "No clock named " << label << " found." << endl;
        return 0;
    }
}


// Resets the clock to 0.
void Timer::ResetClock(string label)
{
    if (Timer::clock_map_.count(label) > 0)
        Timer::clock_map_.at(label).Reset();
    else
        cout << "No clock named " << label << " found." << endl;
}


// Resets all the clocks.
void Timer::ResetAll()
{
    for (map<string, Clock>::iterator it = clock_map_.begin();
         it != clock_map_.end(); ++it)
        it->second.Reset();
}


// Prints the clock in a gtest style.
void Timer::PrintClock(string label)
{
    if (Timer::clock_map_.count(label) > 0)
    {
        cout << "[ Perf     ] Clock "<< label << " is at: ";
        cout << Timer::clock_map_.at(label).GetCurrentTime();
        cout << "s" << endl;
    }
    else
        cout << "No clock named " << label << " found." << endl;

}
