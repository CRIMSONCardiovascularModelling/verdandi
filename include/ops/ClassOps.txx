// Copyright (C) 2010, Vivien Mallet
//
// This file is part of Ops, a library for parsing Lua configuration files.
//
// Ops is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Ops is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Ops. If not, see http://www.gnu.org/licenses/.


#ifndef OPS_FILE_CLASSOPS_TXX


#include "ClassOps.hxx"


namespace Ops
{


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[in] constraint constraint that the entry value must satisfy.
    \param[in] default_value default value for the entry in case it is not
    found in the configuration file.
    \param[out] value value of the entry.
  */
  template<class TD, class T>
  void
  Ops::Set(string name, string constraint, const TD& default_value, T& value)
  {
    SetValue(name, constraint, default_value, true, value);
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[in] constraint constraint that the entry value must satisfy.
    \param[out] value value of the entry.
  */
  template<class T>
  void Ops::Set(string name, string constraint, T& value)
  {
    SetValue(name, constraint, value, false, value);
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[out] value value of the entry.
  */
  template <class T>
  void Ops::Set(string name, T& value)
  {
    SetValue(name, "", value, false, value);
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[in] constraint constraint that the entry value must satisfy.
    \param[in] default_value default value for the entry in case it is not
    found in the configuration file.
    \return The value of the entry.
  */
  template<class T>
  T Ops::Get(string name, string constraint, const T& default_value)
  {
    T value;
    SetValue(name, constraint, default_value, true, value);
    return value;
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \param[in] constraint constraint that the entry value must satisfy.
    \return The value of the entry.
  */
  template<class T>
  T Ops::Get(string name, string constraint)
  {
    T value;
    SetValue(name, constraint, value, false, value);
    return value;
  }


  //! Retrieves a value from the configuration file.
  /*!
    \param[in] name name of the entry.
    \return The value of the entry.
  */
  template <class T>
  T Ops::Get(string name)
  {
    T value;
    SetValue(name, "", value, false, value);
    return value;
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] in parameters of the function.
    \param[out] out outputs of the function.
    \note The prefix is prepended to \a name.
  */
  template<class Tin, class Tout>
  void Ops::Apply(string name, const std::vector<Tin>& in,
                  std::vector<Tout>& out)
  {
    PutOnStack(Name(name));
    PushOnStack(in);

    int n = lua_gettop(state_);

    if (lua_pcall(state_, int(in.size()), LUA_MULTRET, 0) != 0)
      throw Error("Apply(string, vector, vector&)",
                  "While calling " + Function(name) + ":\n  "
                  + lua_tostring(state_, -1));

    n = lua_gettop(state_) - n + 1 + int(in.size());

    out.resize(std::size_t(n));
    for (int i = 0; i < n; i++)
      if (!Convert(i - n, out[std::size_t(i)]))
        {
          std::ostringstream str;
          str << i;
          throw Error("Apply(string, vector, vector&)",
                      "The returned value #" + str.str() + " of \""
                      + Name(name) + "\" is not of correct type.");
        }

    ClearStack();
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \param[in] arg3 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2,
               const T& arg3)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    in.push_back(arg3);
    Apply(name, in, out);
    return out[0];
  }


  //! Applies a Lua function.
  /*!
    \param[in] name name of the function.
    \param[in] arg0 parameter of the function.
    \param[in] arg1 parameter of the function.
    \param[in] arg2 parameter of the function.
    \param[in] arg3 parameter of the function.
    \param[in] arg4 parameter of the function.
    \return First output of the function.
    \note The prefix is prepended to \a name.
  */
  template<class T>
  T Ops::Apply(string name, const T& arg0, const T& arg1, const T& arg2,
               const T& arg3, const T& arg4)
  {
    std::vector<T> in, out;
    in.push_back(arg0);
    in.push_back(arg1);
    in.push_back(arg2);
    in.push_back(arg3);
    in.push_back(arg4);
    Apply(name, in, out);
    return out[0];
  }


  //! Checks whether \a name is of type 'T'.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \return True if the entry is of type 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::Is(string name)
  {
    T value;
    return IsParam(name, value);
  }


  //! Pushes a vector onto the stack.
  /*! Every element of the vector \a v is pushed onto the stack, in the same
    order as in the vector.
    \param[in] v vector to be pushed.
  */
  template<class T>
  void Ops::PushOnStack(const std::vector<T>& v)
  {
    for (std::size_t i = 0; i < v.size(); i++)
      PushOnStack(v[i]);
  }


  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Retrieves a value and checks if it satisfies given constraints.
  /*! If the entry is not found, the default value is returned (if any).
    \param[in] name name of the entry.
    \param[in] constraint constraint to be satisfied.
    \param[in] default_value default value.
    \param[in] with_default is there a default value? If not, \a default_value
    is ignored.
    \param[out] value the value of the entry named \a name.
    \note The default value may not satisfy the constraint.
  */
  template<class TD, class T>
  void Ops::SetValue(string name, string constraint,
                     const TD& default_value, bool with_default,
                     T& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      if (with_default)
        {
          value = default_value;
          ClearStack();
          return;
        }
      else
        throw Error("SetValue",
                    "The " + Entry(name) + " was not found.");

    Convert(-1, value, name);

    if (!CheckConstraint(name, constraint))
      throw Error("SetValue",
                  "The " + Entry(name) + " does not satisfy "
                  + "the constraint:\n" + Constraint(constraint));

    ClearStack();

    Push(Name(name), value);
  }


  //! Retrieves a value and checks if it satisfies given constraints.
  /*! If the entry is not found, the default value is returned (if any).
    \param[in] name name of the entry.
    \param[in] constraint constraint to be satisfied.
    \param[in] default_value default value.
    \param[in] with_default is there a default value? If not, \a default_value
    is ignored.
    \param[out] value the value of the entry named \a name.
    \note The default value may not satisfy the constraint.
  */
  template<class T>
  void Ops::SetValue(string name, string constraint,
                     const std::vector<T>& default_value, bool with_default,
                     std::vector<T>& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      if (with_default)
        {
          value = default_value;
          ClearStack();
          return;
        }
      else
        throw Error("SetValue",
                    "The " + Entry(name) + " was not found.");

    if (!lua_istable(state_, -1))
      throw Error("SetValue",
                  "The " + Entry(name) + " is not a table.");

    std::vector<T> element_list;
    T element;
    std::vector<string> key_list;
    string key;
    // Now loops over all elements of the table.
    lua_pushnil(state_);
    while (lua_next(state_, -2) != 0)
      {
        // Duplicates the key and value so that 'lua_tostring' (applied to
        // them) should not interfere with 'lua_next'.
        lua_pushvalue(state_, -2);
        lua_pushvalue(state_, -2);

        if (!Convert(-2, key))
          throw Error("SetValue",
                      "Unable to read the keys of " + Entry(name) + ".");
        key_list.push_back(key);

        Convert(-1, element, name + "[" + key + "]");
        element_list.push_back(element);

        lua_pop(state_, 3);
      }

    for (std::size_t i = 0; i < key_list.size(); i++)
      if (!CheckConstraint(name + "[" + key_list[i] + "]", constraint))
        throw Error("SetValue",
                    "The " + Entry(name + "[" + key_list[i] + "]")
                    + " does not satisfy the constraint:\n"
                    + Constraint(constraint));

    value = element_list;

    ClearStack();

    Push(Name(name), value);
  }


  //! Checks whether \a name is of type 'T'.
  /*! On exit, the value of the entry (if it exists) is on the stack.
    \param[in] name the name of the entry whose type is checked.
    \param[in] value anything: it is used to determine the type.
    \return True if the entry is of type 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::IsParam(string name, T& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      throw Error("Is(string)",
                  "The " + Entry(name) + " was not found.");

    return Convert(-1, value);
  }


  //! Checks whether \a name is a table of 'T'.
  /*!
    \param[in] name the name of the entry whose type is checked.
    \param[in] value anything: it is used to determine the type.
    \return True if the entry is a table of 'T', false otherwise.
    \note The prefix is prepended to \a name. If \a name does not exist, an
    exception is raised.
  */
  template<class T>
  bool Ops::IsParam(string name, std::vector<T>& value)
  {
    PutOnStack(Name(name));

    if (lua_isnil(state_, -1))
      throw Error("IsParam",
                  "The " + Entry(name) + " was not found.");

    if (!lua_istable(state_, -1))
      return false;

    T element;
    // Now loops over all elements of the table.
    lua_pushnil(state_);
    while (lua_next(state_, -2) != 0)
      {
        // Duplicates the value so that 'lua_tostring' (applied to them)
        // should not interfere with 'lua_next'.
        lua_pushvalue(state_, -1);

        if (!Convert(-1, element))
          return false;
        lua_pop(state_, 2);
      }

    return true;
  }


  //! Pushes all keys of a map into a vector.
  /*!
    \param[in] input the map whose keys should be pushed.
    \param[in,out] vect the vector to which the map keys are pushed back.
  */
  template<class TK, class T>
  void Ops::AppendKey(const std::map<TK, T>& input, std::vector<TK>& vect)
  {
    typename std::map<TK, T>::const_iterator i;
    for (i = input.begin(); i != input.end(); i++)
      vect.push_back(i->first);
  }


}


#define OPS_FILE_CLASSOPS_TXX
#endif
