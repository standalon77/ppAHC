/*
 * SocketException.h
 *
 *  Created on: Jun 18, 2019
 *      Author: firstuser
 */

#ifndef SOCKETEXCEPTION_H_
#define SOCKETEXCEPTION_H_

#include <string>

class SocketException
{
 public:
  SocketException ( std::string s ) : m_s ( s ) {};
  ~SocketException (){};

  std::string description() { return m_s; }

 private:

  std::string m_s;

};

#endif /* SOCKETEXCEPTION_H_ */
