/*
 * ClientSocket.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: firstuser
 */

#include "ClientSocket.h"
#include "SocketException.h"

ClientSocket::ClientSocket() {
	// TODO Auto-generated constructor stub

}

ClientSocket::~ClientSocket() {
	// TODO Auto-generated destructor stub
}

ClientSocket::ClientSocket ( std::string host, int port )
{
  if ( ! Socket::create() )
    {
      throw SocketException ( "Could not create client socket." );
    }

  if ( ! Socket::connect ( host, port ) )
    {
      throw SocketException ( "Could not bind to port." );
    }

}

const ClientSocket& ClientSocket::operator << ( const std::string& s ) const
{
  if ( ! Socket::send ( s ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;

}

const ClientSocket& ClientSocket::operator >> ( std::string& s ) const
{
  if ( ! Socket::recv ( s ) )
    {
      throw SocketException ( "Could not read from socket." );
    }

  return *this;
}

const int ClientSocket::recv ( void* s, const size_t size_s ) const
{
  int size_recv = Socket::recv ( s, size_s );
  if ( ! size_recv )
    {
      throw SocketException ( "Could not read from socket." );
    }

  return size_recv;
}

const ClientSocket& ClientSocket::send ( const void* s, const size_t size_s ) const
{
  if ( ! Socket::send ( s, size_s ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;

}

