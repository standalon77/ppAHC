/*
 * PaillierCrypto.h
 *
 *  Created on: Jun 18, 2019
 *      Author: PJS
 */

#ifndef PAILLIERCRYPTO_H_
#define PAILLIERCRYPTO_H_

#define _DEBUG_INIT_1
//#define _DEBUG_Assert
//#define _DEBUG_THREAD
//#define _DEBUG_Initialization
//#define _DEBUG_SM
//#define _DEBUG_LSB
//#define _DEBUG_SVR
//#define _DEBUG_SingleOutput
//#define _DEBUG_SEQ1
//#define _DEBUG_SEQ2
//#define _DEBUG_STEP2
//#define _DEBUG_TERMINATE_PGM
//#define _DEBUG_Communication

#include <iostream>
#include <string>
#include <gmp.h>
extern "C"{
	#include "paillier.h"
}
#include "ClientSocket.h"
#include <assert.h>
#include <queue>
#include <mutex>
#include <condition_variable>

typedef paillier_ciphertext_t cTxt_t;
typedef paillier_plaintext_t  pTxt_t;

// Command Tag
const unsigned char COM_MUL1 	= 0x01;
const unsigned char COM_MUL2 	= 0x02;
const unsigned char COM_LSB 	= 0x03;
const unsigned char COM_SVR 	= 0x04;
const unsigned char COM_SEQ1 	= 0x05;
const unsigned char COM_SEQ2 	= 0x06;
const unsigned char COM_SING 	= 0x07;
const unsigned char COM_STP2 	= 0x08;
const unsigned char COM_ST2D1 	= 0x09;
const unsigned char COM_ST2D2 	= 0x0A;
const unsigned char COM_ST2V 	= 0x0B;
const unsigned char COM_TERM 	= 0xFF;

// structure for recv thread
typedef struct {
	std::mutex m;
	std::condition_variable cv;
	unsigned char* pa[THREAD_NUM];
} recv_t;

class PaillierCrypto {
private:
	ClientSocket 			*mCSocket;
	paillier_pubkey_t		*mPubKey;
	paillier_prvkey_t		*mPrvKey;
	int mSt2Dst, mSt2Src;

	inline void SetSendMsg(unsigned char* ucSendPtr, const unsigned char* ucRecvPtr, unsigned char* Data, const unsigned short Len);
	inline void SetSendMsg(unsigned char* ucSendPtr, const unsigned char* ucRecvPtr, unsigned char bData);
	inline void SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data, const unsigned short Idx, const unsigned char Tag, const unsigned short Len);
	inline unsigned short Byte2Short(unsigned char* In);
	inline void DebugOut(const char* pMsg, const mpz_t pData, 			const short idx, const short DHidx);
	inline void DebugOut(const char* pMsg, short pData, 					const short idx, const short DHidx);
	inline void DebugDec(const char* pMsg, cTxt_t* cTxt, const short idx, const short DHidx);
	inline void DebugCom(const char* pMsg, unsigned char *data, int len, const short idx, const short DHidx);
	inline void paillier_ciphertext_from_bytes(cTxt_t* ct, void* c, int len );
	inline void paillier_ciphertext_to_bytes(unsigned char* buf, int len, cTxt_t* ct );

public:
	PaillierCrypto();
	PaillierCrypto(int pModSize);
	virtual ~PaillierCrypto();

	bool SecMul1	 	(unsigned short idx, unsigned char* ucRecvPtr);
	bool SecMul2	 	(unsigned short idx, unsigned char* ucRecvPtr);
	void EncryptedLSB	(unsigned short idx, unsigned char* ucRecvPtr);
	void SVR		 	(unsigned short idx, unsigned char* ucRecvPtr);
	bool SingleOutput	(unsigned short idx, unsigned char* ucRecvPtr);
	bool SEQ1	 		(unsigned short idx, unsigned char* ucRecvPtr);
	int  SEQ2			(unsigned short idx, unsigned char* ucRecvPtr);
	int  Step_2		(unsigned short idx, unsigned char* ucRecvPtr, recv_t* tRecv);
	int  Step_2_Delta	(unsigned short idx, unsigned char* ucRecvPtr);
	int  Step_2_v		(unsigned short idx, unsigned char* ucRecvPtr);
	void TerminatePgm	(unsigned short idx, unsigned char* ucRecvPtr);

	int Receiver(recv_t* tRecv);
	bool distributePubKey();
	bool distributePrvKey();
	void DebugComMain(const char* pMsg, unsigned char *data, int len, const short idx, const short DHidx);
};

#endif /* PAILLIERCRYPTO_H_ */
