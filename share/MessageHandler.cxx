// Copyright (C) 2009, INRIA
// Author(s): Claire Mouton
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


#ifndef VERDANDI_FILE_SHARE_MESSAGEHANDLER_CXX
#define VERDANDI_FILE_SHARE_MESSAGEHANDLER_CXX


#include "VerdandiHeader.hxx"


namespace Verdandi
{
    /////////////////////
    // MESSAGE HANDLER //
    /////////////////////

    //! Adds a new object in the recipient list.
    /*!
      \param[in] recipient the string describing the object to add.
      \param[in] object pointer to the recipient object.
      \param[in] pointer the pointer to the method to add.
    */
    void MessageHandler
    ::AddRecipient(string recipient, void* object,
                   MessageHandler::function_pointer pointer)
    {
        pair<void*, function_pointer> pointer_pair(object, pointer);
        Recipient_map()[recipient].push_back(pointer_pair);
        Recipient_map()["all"].push_back(pointer_pair);
    }


    //! Sends a message to a recipient.
    /*!
      \param[in] recipient the recipient of the message.
      \param[in] message the string containing the message.
    */
    void MessageHandler
    ::Send(string recipient, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        Logger::Log<-10>(MessageHandler::GetName(), string("Message \"")
                         + message + "\" is sent to \"" + recipient + "\".");

        if (Recipient_map().count(recipient) == 0)
            throw ErrorArgument("MessageHandler::Send",
                                string("The object \"") + recipient +
                                "\" is not part of the recipient list.");

        recipient_list my_list = Recipient_map()[recipient];
        SendToList(my_list, message);
#endif
    }


    //! Returns the name of the class, that is, "MessageHandler".
    /*!
      \return The name of the class.
    */
    string MessageHandler::GetName()
    {
        return "MessageHandler";
    }


    //! Sends a message to a list of recipients.
    /*!
      \param[in] my_list the list of recipients.
      \param[in] message the string containing the message.
    */
    void MessageHandler::SendToList(recipient_list& my_list, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        recipient_list::iterator my_iterator;
        for(my_iterator = my_list.begin(); my_iterator != my_list.end();
            ++my_iterator)
            (*my_iterator->second)(my_iterator->first, message);
#endif
    }

    //! Returns a reference to a recipient map
    /*!
      \return The recipient_map
    */
    MessageHandler::recipient_map& MessageHandler::Recipient_map()
    {
        static recipient_map* output = new recipient_map;
        return *output;
    }


} // namespace Verdandi.


#endif
