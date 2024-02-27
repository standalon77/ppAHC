// Definition of the Socket class

#ifndef Socket_class
#define Socket_class

#define DATA_SQUARE_LENGTH 4					// (check) bits in squared data
#define DATA_NUMBER_LENGTH 7					// (check) bit in number of data
const int DATA_NUM = 11;							// (check)
const int DATA_DIM = 21;							// (check)
const int PARAM_K =  1;							// (check)
const int THREAD1_NUM = 5;						// (check) number of Main threads
const int MOD_SIZE = 512;						// (check) bit

const int KEY_SIZE = MOD_SIZE/8;				// byte: N
const int ENC_SIZE = KEY_SIZE*2;				// byte: N^2
const int GMP_N_SIZE = MOD_SIZE/64;
const int HED_SIZE = 5;
const int HED_LEN  = 3;
const int DTRI_SIZ = DATA_NUM*(DATA_NUM-1)/2;	// upper triangular size of DT


#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <string>
#include <arpa/inet.h>

const int MAXHOSTNAME = 200;
const int MAXCONNECTIONS = 5;
const int MAXRECV = HED_SIZE+DATA_NUM*ENC_SIZE+10;


class Socket
{
 public:
  Socket();
  virtual ~Socket();

  // Server initialization
  bool create();
  bool bind ( const int port );
  bool listen() const;
  bool accept ( Socket& ) const;

  // Client initialization
  bool connect ( const std::string host, const int port );

  // Data Transimission
  bool send ( const std::string ) const;
  bool send ( const void* s, const size_t size_s ) const;
  int recv ( std::string& ) const;
  int recv ( void* s, const size_t size_s ) const;


  void set_non_blocking ( const bool );

  bool is_valid() const { return m_sock != -1; }

 private:

  int m_sock;
  sockaddr_in m_addr;


};


#endif
