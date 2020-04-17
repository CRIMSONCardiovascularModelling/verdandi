// Copyright (C) 2010, INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_SHARE_VERDANDIOPS_CXX
#define VERDANDI_FILE_SHARE_VERDANDIOPS_CXX


#include "VerdandiHeader.hxx"
#include "ops/Ops.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! Nothing is performed. The Lua state is set to NULL.
     */
    VerdandiOps::VerdandiOps(): ::Ops::Ops()
    {
    }


    //! Main constructor.
    /*! The Lua configuration file is loaded and run. An exception may be
      raised during this evaluation.
      \param[in] file_path path to the configuration file.
    */
    VerdandiOps::VerdandiOps(string file_path): ::Ops::Ops(file_path)
    {
    }


    //! Destructor.
    /*! Destroys the Lua state object.
     */
    VerdandiOps::~VerdandiOps()
    {
    }


} // namespace Verdandi.


#endif
