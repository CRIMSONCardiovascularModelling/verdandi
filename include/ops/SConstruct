# Copyright (C) 2010 Vivien Mallet
#
# This file is part of Ops, a library for parsing Lua configuration files.
#
# Ops is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Ops is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Ops. If not, see http://www.gnu.org/licenses/.


import distutils.sysconfig, os, glob


####################
# USEFUL FUNCTIONS #
####################

# This function checks whether a supported argument has been provided, and it
# sets the argument to its default value if the user has not given any value.
def add_argument(name, value_list = None):
    if value_list is None:
        if not ARGUMENTS.has_key(name):
            raise Exception, "The command line argument \"" + name + "\" is" \
                  + " required, but it was not provided."
    else:
        ARGUMENTS[name] = ARGUMENTS.get(name, value_list[0])
        if ARGUMENTS[name] not in value_list:
            raise Exception, "Unsupported option \"" + ARGUMENTS[name] \
                  + "\" for argument \"" + name + "\". Available options " \
                  + "are: " \
                  + ", ".join(["\"" + x + "\"" for x in value_list]) + "."


env = Environment(ENV = os.environ,
                  SWIGFLAGS = ['-Wall', '-c++', '-python' ],
                  CPPPATH = [distutils.sysconfig.get_python_inc(),
                             "/usr/include/lua5.1/"],
                  SHLIBPREFIX = "")

conf = Configure(env)

# Link to the appropriate version of Python.
conf.CheckLib("python" + distutils.sysconfig.get_python_version())
if not conf.CheckLib("lua5.1"):
    conf.CheckLib("lua")

env.Append(CPPFLAGS = " -DOPS_WITH_EXCEPTION")
env.Append(CPPFLAGS = " -fPIC") # for SWIG.

if env['PLATFORM'] == 'win32':
	env.Append(SHLIBSUFFIX = ".pyd")
	env.Replace(LINK = "LINK")
	env.SharedLibrary('_ops', ['Ops.cpp', 'ops.i'])
else:
    env.SharedLibrary('_ops.so', ['Ops.cpp', 'ops.i'])


env_lib = env.Clone()
add_argument("lib", ["no", "yes"])
if ARGUMENTS["lib"] == "yes":
    library_file = glob.glob("*.cxx")
    library_path = "lib/libops.a"
    env_lib.Library(library_path, library_file)
    env_lib.Alias("libops.a", library_path)
