// Copyright (C) 2008, INRIA
// Author(s): Marc Fragu
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


#ifndef VERDANDI_FILE_SHARE_LOGGER_TXX
#define VERDANDI_FILE_SHARE_LOGGER_TXX


namespace Verdandi
{


    //! Writes a message in the log file.
    /*!
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <int LEVEL, class T, class S>
    void Logger::Log(const T& object, const S& message, int options)
    {
        if (!CheckStatus(options))
            return;
        if (LEVEL >= logging_level_)
            Log(object, to_str(message), options);
    }


    //! Writes a message in the log file.
    /*!
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <class T, class S>
    void Logger::Log(const T& object, const S& message, int options)
    {
        Log(object, to_str(message), options);
    }


    //! Writes a message in the log file.
    /*!
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */
    template <class T>
    void Logger::Log(const T& object, string message, int options)
    {
        if (!CheckStatus(options))
            return;
        WriteMessage(object, message, options);
    }

    //! Writes a message in the standard output and in the log file.
    /*! The message is always sent to the standard output, and it is possibly
      written in a log file if the logging level is lower than or equal to
      'VERDANDI_STDOUT_LOGGING_LEVEL'.
      \tparam S type of the message, which must be convertible to a string
      through 'ostringstream& operator << (ostringstream&, S& message)'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
    */
    template <class T, class S>
    void Logger::StdOut(const T& object, const S& message)
    {
        StdOut(object, to_str(message));
    }


    //! Writes a message in the standard output and in the log file.
    /*! The message is always sent to the standard output, and it is possibly
      written in a log file if the logging level is lower than or equal to
      'VERDANDI_STDOUT_LOGGING_LEVEL'.
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
    */
    template <class T>
    void Logger::StdOut(const T& object, string message)
    {
        WriteMessage(object, message, stdout_);
        if (VERDANDI_STDOUT_LOGGING_LEVEL >= logging_level_)
            WriteMessage(object, message, options_ & ~stdout_);
    }

    //! Writes a message in the log file.
    /*!
      \param[in] object the object that sends the message.
      \param[in] message the message to be written.
      \param[in] options options.
    */


    template <class T>
    void Logger::WriteMessage(const T& object, string message, int options)
    {
        WriteMessage(object.GetName(), message, options);
    }


}


#endif
