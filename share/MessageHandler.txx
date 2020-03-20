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


#ifndef VERDANDI_FILE_SHARE_MESSAGEHANDLER_TXX
#define VERDANDI_FILE_SHARE_MESSAGEHANDLER_TXX


namespace Verdandi
{




    /////////////////////
    // MESSAGE HANDLER //
    /////////////////////


    //! Adds a new object in the recipient list.
    /*!
      \tparam R type of the recipient object.
      \param[in] recipient the string describing the object to add.
      \param[in] object reference to the recipient object.
      \param[in] pointer the pointer to the method to add.
    */
    template <class R>
    void MessageHandler
    ::AddRecipient(string recipient, R& object,
                   MessageHandler::function_pointer pointer)
    {
        pair<void*, function_pointer>
            pointer_pair(reinterpret_cast<void*>(&object), pointer);
        Recipient_map()[recipient].push_back(pointer_pair);
        Recipient_map()["all"].push_back(pointer_pair);
    }


    //! Sends a message to a recipient and specifies the sender.
    /*!
      \param[in] sender the object sending the message.
      \param[in] recipient the recipient of the message.
      \param[in] message the string containing the message.
      \tparam Sender the type of the sending object.
    */
    template <class Sender>
    void MessageHandler
    ::Send(const Sender& sender, string recipient, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        Logger::Log<-10>(MessageHandler::GetName(), string("Message \"")
                         + message + "\" is sent to \"" + recipient
                         + "\" by \"" + sender.GetName() + "\".");

        if (Recipient_map().count(recipient) == 0)
            throw ErrorArgument("MessageHandler::Send",
                                string("The object \"") + recipient +
                                "\" is not part of the recipient list.");

        recipient_list my_list = Recipient_map()[recipient];
        SendToList(my_list, string("[") + sender.GetName() + "] "
                   + message);
#endif
    }


    //! Removes an object from the static recipient map.
    /*!
      \tparam R type of the recipient object.
      \param[in] object the object to remove from the list.
    */
    template <class R>
    void MessageHandler::RemoveRecipient(R& object)
    {
        recipient_map::iterator my_map_iterator;
        for (my_map_iterator = Recipient_map().begin();
             my_map_iterator != Recipient_map().end(); )
        {
            recipient_list::iterator my_map_list_iterator;
            recipient_list my_list;
            for (my_map_list_iterator = my_map_iterator->second.begin();
                 my_map_list_iterator != my_map_iterator->second.end();
                 ++my_map_list_iterator)
                if (my_map_list_iterator->first
                    != reinterpret_cast<void*>(&object))
                    my_list.push_back(*my_map_list_iterator);

            // 'my_list' now contains all remaining objects for the current
            // entry in the map. If it is empty, the map entry should be
            // removed.
            if (my_list.empty())
                Recipient_map().erase(my_map_iterator++);
            else
            {
                my_map_iterator->second = my_list;
                ++my_map_iterator;
            }
        }
    }


} // namespace Verdandi.


#endif
