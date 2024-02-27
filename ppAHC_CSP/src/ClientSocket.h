/*
 * ClientSocket.h
 *
 *  Created on: Jun 18, 2019
 *      Author: firstuser
 */

#ifndef CLIENTSOCKET_H_
#define CLIENTSOCKET_H_

#include "Socket.h"

class ClientSocket: private Socket {
public:
	ClientSocket();
	virtual ~ClientSocket();
	ClientSocket(std::string host, int port);
  const ClientSocket& operator << ( const std::string& ) const;
  const ClientSocket& operator >> ( std::string& ) const;
  const ClientSocket& send ( const void* s, const size_t size_s ) const;
  const int recv ( void* s, const size_t size_s ) const;

};

#endif /* CLIENTSOCKET_H_ */
