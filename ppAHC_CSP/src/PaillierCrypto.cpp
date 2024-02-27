/*
 * PaillierCrypto.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: PJS
 */


#include "PaillierCrypto.h"
#include <iostream>
#include <string.h>
#include "SocketException.h"
#include <cstring>

PaillierCrypto::PaillierCrypto(): mCSocket(NULL), mPubKey(NULL), mPrvKey(NULL)
{
	mSt2Dst = 0;
	mSt2Src = 0;
}

PaillierCrypto::PaillierCrypto(int pModSize)
{
    // Generate public and secret keys
    paillier_keygen(pModSize, &mPubKey, &mPrvKey, paillier_get_rand_devurandom);

	#ifdef _DEBUG_Initialization
    gmp_printf("[CSP] Public key(n)\t\t: %ZX\n", mPubKey->n);
    gmp_printf("[CSP] Private key(lambda)\t: %ZX\n", mPrvKey->lambda);
	#endif

	try {
		// connect to DH
		//mCSocket = new ClientSocket("192.168.25.2", 3000);
		mCSocket = new ClientSocket("localhost", 3000);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	mSt2Dst = 0;
	mSt2Src = 0;

	return;
}

PaillierCrypto::~PaillierCrypto()
{
	paillier_freepubkey(mPubKey);
	paillier_freeprvkey(mPrvKey);
	delete(mCSocket);
}

// send public key to DH
bool PaillierCrypto::distributePubKey()
{
	// distribute public key
	char* cpPubKey = paillier_pubkey_to_hex(mPubKey);
	#ifdef _DEBUG_Initialization
    std::cout << "[CSP] Public Key (Hex)\t\t: " << cpPubKey << std::endl;
	#endif

	try {
		mCSocket->send(cpPubKey, KEY_SIZE*2);		// string으로 전송
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	delete[] cpPubKey;

	return true;
}

// for debugging
bool PaillierCrypto::distributePrvKey()
{
	// distribute private key
	char* cpPrvKey = paillier_prvkey_to_hex(mPrvKey);
	#ifdef _DEBUG_Initialization
    std::cout << "[CSP] Private Key (Hex)\t\t: " << cpPrvKey << std::endl;
	#endif

	try {
		mCSocket->send(cpPrvKey, KEY_SIZE*2);		// string으로 전송
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	delete[] cpPrvKey;

	return true;
}

// secure multiplication protocol
bool PaillierCrypto::SecMul1(unsigned short idx, unsigned char* ucRecvPtr)
{
	cTxt_t cEh, cEa;
	pTxt_t  pHa, pH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cEh.c, cEa.c, pHa.m, pH.m, NULL);
	#else
	mpz_init2(cEh.c, 2*GMP_N_SIZE*2);
	mpz_init2(cEa.c, 2*GMP_N_SIZE);
	mpz_init2(pHa.m, 2*GMP_N_SIZE+1);
	mpz_init2(pH.m, 2*GMP_N_SIZE);
	#endif

	unsigned char bEh[ENC_SIZE]={0,};
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SM
	DebugCom("Encrypted a' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == ENC_SIZE);
	#endif

	// a' = E(a+r_a)
	paillier_ciphertext_from_bytes(&cEa, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
	DebugOut("Encrypted a' (Copied)\t\t\t", cEa.c, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// h_a = D(a')
	paillier_dec(&pHa, mPubKey, mPrvKey, &cEa);
	#ifdef _DEBUG_SM
    DebugOut("Decrypted ha\t\t\t", pHa.m, idx, sDHidx);
	#endif

    // h = h_a^2 mod N
    mpz_mul(pH.m, pHa.m, pHa.m);
    mpz_mod(pH.m, pH.m, mPubKey->n);
	#ifdef _DEBUG_SM
    DebugOut("h = ha^2 mod N\t\t\t", pH.m, idx, sDHidx);
	#endif

    // h' = E(h)
	paillier_enc(&cEh, mPubKey, &pH, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
    DebugOut("h' = E(h)\t\t\t", cEh.c, idx, sDHidx);
	#endif

	//------------------------------------------------------------------------------
    // send h'
	paillier_ciphertext_to_bytes(bEh, ENC_SIZE, &cEh);
	#ifdef _DEBUG_SM
    DebugCom("h' = E(h) (Hex)\t\t\t", bEh, ENC_SIZE, idx, sDHidx);
	#endif

    sLen = ENC_SIZE;
	SetSendMsg(ucSendPtr, ucRecvPtr, bEh, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_SM
    DebugCom("h' = E(h) (copied)\t\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pHa.m->_mp_alloc == 2*GMP_N_SIZE+1);
	assert( pH.m->_mp_alloc == 2*GMP_N_SIZE);
	assert(cEh.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(cEa.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+2*ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return true;
}

// secure multiplication protocol (square)
bool PaillierCrypto::SecMul2(unsigned short idx, unsigned char* ucRecvPtr)
{
	cTxt_t cEh, cEa, cEb;
	pTxt_t  pHa, pHb, pH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cEh.c, cEa.c, cEb.c, pHa.m, pHb.m, pH.m, NULL);
	#else
	mpz_init2(cEh.c, 2*GMP_N_SIZE*2);
	mpz_init2(cEa.c, 2*GMP_N_SIZE);
	mpz_init2(cEb.c, 2*GMP_N_SIZE);
	mpz_init2(pHa.m, 2*GMP_N_SIZE+1);
	mpz_init2(pHb.m, 2*GMP_N_SIZE+1);
	mpz_init2(pH.m, 2*GMP_N_SIZE);
	#endif

	unsigned char bEh[ENC_SIZE]={0, };
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SM
	DebugCom("Encrypted a' and b' (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == 2*ENC_SIZE);
	#endif

	// a' = E(a+r_a)
	paillier_ciphertext_from_bytes(&cEa, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
	DebugOut("Encrypted a' (Copied)\t\t", cEa.c, idx, sDHidx);
	#endif

	// b' = E(b+r_b)
	paillier_ciphertext_from_bytes(&cEb, ucRecvPtr+HED_SIZE+ENC_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
	DebugOut("Encrypted b' (Copied)\t\t", cEb.c, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// h_a = D(a'), h_b = D(b')
	paillier_dec(&pHa, mPubKey, mPrvKey, &cEa);
	#ifdef _DEBUG_SM
    DebugOut("Decrypted ha\t\t\t", pHa.m, idx, sDHidx);
	#endif

    paillier_dec(&pHb, mPubKey, mPrvKey, &cEb);
    #ifdef _DEBUG_SM
    DebugOut("Decrypted hb\t\t\t", pHb.m, idx, sDHidx);
	#endif

    // h = h_a * h_b  mod N
    mpz_mul(pH.m, pHa.m, pHb.m);
    mpz_mod(pH.m, pH.m, mPubKey->n);
	#ifdef _DEBUG_SM
    DebugOut("h = ha * hb mod N\t\t\t", pH.m, idx, sDHidx);
	#endif

    // h' = E(h)
	paillier_enc(&cEh, mPubKey, &pH, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
    DebugOut("h' = E(h)\t\t\t\t", cEh.c, idx, sDHidx);
	#endif

	//------------------------------------------------------------------------------
    // send h'
	paillier_ciphertext_to_bytes(bEh, ENC_SIZE, &cEh);
	#ifdef _DEBUG_SM
    DebugCom("h' = E(h) (Hex)\t\t\t", bEh, ENC_SIZE, idx, sDHidx);
	#endif

    sLen = ENC_SIZE;
	SetSendMsg(ucSendPtr, ucRecvPtr, bEh, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_SM
    DebugCom("h' = E(h) (copied)\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pHa.m->_mp_alloc == 2*GMP_N_SIZE+1);
	assert(pHb.m->_mp_alloc == 2*GMP_N_SIZE+1);
	assert( pH.m->_mp_alloc == 2*GMP_N_SIZE);
	assert(cEh.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(cEa.c->_mp_alloc == 2*GMP_N_SIZE);
	assert(cEb.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return true;
}

// subprotocol in secure big-decomposition protocol
void PaillierCrypto::EncryptedLSB(unsigned short idx, unsigned char* ucRecvPtr)
{
	pTxt_t pY, p;
	cTxt_t c, cY;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pY.m, p.m, c.c, cY.c, NULL);
	#else
	mpz_init2(pY.m, 2*GMP_N_SIZE+1);
	mpz_init2(p.m, 1);
	mpz_init2(c.c, 3*GMP_N_SIZE);
	mpz_init2(cY.c, 2*GMP_N_SIZE);
	#endif

	int iEven;
	unsigned char bAlpha[ENC_SIZE] = {0, };
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE] = {0, };
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	// Y = E(x+r)
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_LSB
    DebugCom("Y = E(x+r) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == ENC_SIZE);
	#endif

	paillier_ciphertext_from_bytes(&cY, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_LSB
	DebugOut("Y = E(x+r) (Copied)\t\t", cY.c, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// y = D(Y) = x+r
	paillier_dec(&pY, mPubKey, mPrvKey, &cY);
	#ifdef _DEBUG_LSB
    DebugOut("Decrypted y = x+r\t\t\t", pY.m, idx, sDHidx);
	#endif

	// determine whether random y is even
	iEven = mpz_even_p(pY.m);
	#ifdef _DEBUG_LSB
    DebugOut("y is even?\t\t\t", iEven, idx, sDHidx);
	#endif

    // send alpha = E(0) or E(1)
	if (iEven) {
		mpz_set_ui(p.m, 0);
		paillier_enc(&c, mPubKey, &p, paillier_get_rand_devurandom);
		//------------------------------------------------------------------------------
		paillier_ciphertext_to_bytes(bAlpha, ENC_SIZE, &c);
		#ifdef _DEBUG_LSB
	    DebugCom("(y:even) alpha = E(0) (Hex)\t", bAlpha, ENC_SIZE, idx, sDHidx);
		#endif
	}
	else {
		mpz_set_ui(p.m, 1);
		paillier_enc(&c, mPubKey, &p, paillier_get_rand_devurandom);
		//------------------------------------------------------------------------------
		paillier_ciphertext_to_bytes(bAlpha, ENC_SIZE, &c);
		#ifdef _DEBUG_LSB
	    DebugCom("(y:odd) alpha = E(1) (Hex)\t", bAlpha, ENC_SIZE, idx, sDHidx);
		#endif
	}

	//------------------------------------------------------------------------------
    sLen = ENC_SIZE;
	SetSendMsg(ucSendPtr, ucRecvPtr, bAlpha, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_LSB
    DebugCom("h' = E(h) (copied)\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(pY.m->_mp_alloc == 2*GMP_N_SIZE+1);
	assert( p.m->_mp_alloc == 1);
	assert( c.c->_mp_alloc <= 3*GMP_N_SIZE);
	assert(cY.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+2*ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

// subprotocol in secure big-decomposition protocol
void PaillierCrypto::SVR(unsigned short idx, unsigned char* ucRecvPtr)
{
	pTxt_t pW;
	cTxt_t cW;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pW.m, cW.c, NULL);
	#else
	mpz_init2(pW.m, GMP_N_SIZE);
	mpz_init2(cW.c, 2*GMP_N_SIZE);
	#endif

	unsigned char ucSendPtr[HED_SIZE+1];
	int iZero;
	char bGamma;
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	// W
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SVR
    DebugCom("W (Hex)\t\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == ENC_SIZE);
	#endif

	paillier_ciphertext_from_bytes(&cW, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SVR
	DebugOut("W (Copied)\t\t\t", cW.c, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// w = D(W)
	paillier_dec(&pW, mPubKey, mPrvKey, &cW);
	#ifdef _DEBUG_SVR
    DebugOut("Decrypted w\t\t\t", pW.m, idx, sDHidx);
	#endif

    // send gamma (result)
    iZero = mpz_cmp_ui(pW.m, 0);
    if (iZero==0)  	{ bGamma = 1; }	// Success
    else	    	{ bGamma = 0; }		// Fail

	//------------------------------------------------------------------------------
	SetSendMsg(ucSendPtr, ucRecvPtr, bGamma);
    sLen = HED_SIZE+1;
	#ifdef _DEBUG_SVR
	DebugCom("Gamma (copied)\t\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	#ifdef _DEBUG_Assert
//	assert(pW.m->_mp_alloc ==   GMP_N_SIZE);
//	assert(cW.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+2*ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

// SMin protocol
bool PaillierCrypto::SingleOutput(unsigned short idx, unsigned char* ucRecvPtr)
{
	pTxt_t  pX, pY;
	cTxt_t cX, cY;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pX.m, pY.m, cX.c, cY.c, NULL);
	#else
	mpz_init2(pX.m, 2*GMP_N_SIZE+1);
	mpz_init2(pY.m, 1);
	mpz_init2(cX.c, 2*GMP_N_SIZE);
	mpz_init2(cY.c, 3*GMP_N_SIZE);
	#endif

	unsigned char bY[ENC_SIZE] = {0, };
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE] = {0, };
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	// X'
    sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SingleOutput
	DebugCom("(buf) X' (Hex)\t\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
    assert(sLen == ENC_SIZE);
	#endif

	paillier_ciphertext_from_bytes(&cX, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SingleOutput
	DebugDec("X' \t\t\t", &cX, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// x = D(x')
	paillier_dec(&pX, mPubKey, mPrvKey, &cX);
	#ifdef _DEBUG_SingleOutput
	DebugOut("Decrypted x = D(x')\t\t", pX.m, idx, sDHidx);
	#endif

	// (x==0) y=E(1), (x!=0) y=E(0)
	if (mpz_cmp_ui(pX.m, 0) == 0)
		mpz_set_ui(pY.m, 1);
	else
		mpz_set_ui(pY.m, 0);

	paillier_enc(&cY, mPubKey, &pY, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SingleOutput
	//DebugOut("(x==0) y=E(1), (x!=0) y=E(0)", cY.c, idx, sDHidx);
	DebugDec("(x==0) y=E(1), (x!=0) y=E(0)\t", &cY, idx, sDHidx);
	#endif

	//------------------------------------------------------------------------------
	// send Y'
	paillier_ciphertext_to_bytes(bY, ENC_SIZE, &cY);
	#ifdef _DEBUG_SingleOutput
	DebugCom("(buf) Y' (Hex)\t\t\t", bY, ENC_SIZE, idx, sDHidx);
	#endif

    sLen = ENC_SIZE;
	SetSendMsg(ucSendPtr, ucRecvPtr, bY, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_SingleOutput
    DebugCom("(buf) y' = E(y) (copied)\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	#ifdef _DEBUG_Assert
//	assert(pX.m->_mp_alloc <= 2*GMP_N_SIZE+1);
//	assert(pV.m->_mp_alloc == 1);
//	assert(cV.c->_mp_alloc <= 3*GMP_N_SIZE);
//	assert(cX.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return true;
}

// secure equality protocol
bool PaillierCrypto::SEQ1(unsigned short idx, unsigned char* ucRecvPtr)
{
	cTxt_t	cX, cXi[DATA_NUMBER_LENGTH];
	pTxt_t	pX, pXi[DATA_NUMBER_LENGTH];

	#ifdef _DEBUG_INIT_1
	mpz_inits(cX.c, pX.m, NULL);
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
		mpz_inits(cXi[i].c, pXi[i].m, NULL);
	#else
	mpz_init2(cX.c, 2*GMP_N_SIZE*2);
	mpz_init2(pX.m, 2*GMP_N_SIZE);
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		mpz_init2(cXi[i].c, 2*GMP_N_SIZE*2);
		mpz_init2(pXi[i].m, 2*GMP_N_SIZE);
	}
	#endif

	unsigned char bX[ENC_SIZE*DATA_NUMBER_LENGTH]={0,};
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE*DATA_NUMBER_LENGTH]={0,};
	short sLen, sDHidx;

	//------------------------------------------------------------------------------
	// X'
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SEQ1
	DebugCom("(buf) X' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == ENC_SIZE);
	#endif

	// receive x' = E(x)
	paillier_ciphertext_from_bytes(&cX, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SEQ1
	DebugOut("x' = E(x) \t\t\t", cX.c, idx, sDHidx);
	#endif
	//------------------------------------------------------------------------------

	// x = D(x')
	paillier_dec(&pX, mPubKey, mPrvKey, &cX);
	#ifdef _DEBUG_SEQ1
    DebugOut("Decrypted x\t\t\t", pX.m, idx, sDHidx);
	#endif

	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		mpz_fdiv_r_ui(pXi[i].m, pX.m, 2);
		paillier_enc(&cXi[i], mPubKey, &pXi[i], paillier_get_rand_devurandom);

		mpz_sub(pX.m, pX.m, pXi[i].m);
		mpz_tdiv_q_2exp(pX.m, pX.m, 1);
	}

	// (for debugging) decomposed bit x[i]
	#ifdef _DEBUG_SEQ1
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
	    DebugOut("pXi[i] \t\t\t", cXi[i].c, idx, sDHidx);
	    DebugDec("pXi[i] \t\t\t", cXi+i, idx, sDHidx);
	}
	#endif

	//------------------------------------------------------------------------------
    // send E(x[])
    for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		paillier_ciphertext_to_bytes(bX+(ENC_SIZE*i), ENC_SIZE, cXi+i);
		#ifdef _DEBUG_SEQ1
		DebugCom("(buf) x' = E(x) (Hex)\t\t\t", bX, ENC_SIZE*DATA_NUMBER_LENGTH, idx, sDHidx);
		#endif
    }

    sLen = ENC_SIZE*DATA_NUMBER_LENGTH;
	SetSendMsg(ucSendPtr, ucRecvPtr, bX, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_SEQ1
    DebugCom("(buf) E(x[]) (copied)\t\t\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	#ifdef _DEBUG_Assert
//	assert(pHa.m->_mp_alloc == 2*GMP_N_SIZE+1);
//	assert( pH.m->_mp_alloc == 2*GMP_N_SIZE);
//	assert(cEh.c->_mp_alloc == 2*GMP_N_SIZE*2);
//	assert(cEa.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return true;
}

// secure equality protocol
int PaillierCrypto::SEQ2(unsigned short idx, unsigned char* ucRecvPtr)
{
	pTxt_t  pD, pDelta;
	cTxt_t cD, cDelta;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pD.m, pDelta.m, cD.c, cDelta.c, NULL);
	#else
	mpz_init2(pD.m, 2*GMP_N_SIZE+1);
	mpz_init2(pDelta.m, 1);
	mpz_init2(cD.c, 2*GMP_N_SIZE);
	mpz_init2(cDelta.c, 2*GMP_N_SIZE);
	#endif

	unsigned char bDelta[ENC_SIZE] = {0, };
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE] = {0, };
	short sLen, sDHidx;
	int iZero=0;

	//------------------------------------------------------------------------------
   sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SEQ2
	DebugCom("(buf) cDi[size]' (Hex)\t\t\t\t", ucRecvPtr, sLen+HED_SIZE*DATA_NUMBER_LENGTH, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert((*(ucRecvPtr+2)==COM_SEQ2)&&(sLen == ENC_SIZE*DATA_NUMBER_LENGTH));
	#endif

	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		// receive D' = E(D)
		paillier_ciphertext_from_bytes(&cD, ucRecvPtr+HED_SIZE+i*ENC_SIZE, ENC_SIZE);
		#ifdef _DEBUG_SEQ2
		DebugOut("Encrypted cDi[i] (Copied)\t\t", cD.c, idx, sDHidx);
		#endif
		//------------------------------------------------------------------------------

		// cDi[i]' = D(cDi[i])
		paillier_dec(&pD, mPubKey, mPrvKey, &cD);
		#ifdef _DEBUG_SEQ2
		DebugOut("Decrypted cDi[i]\t\t\t", pD.m, idx, sDHidx);
		#endif

		// compare with 0,
		if (mpz_cmp_d(pD.m, 0) == 0)
			iZero++;
	}

	// check whether there exist 0
	if (iZero > 0) 		mpz_set_ui(pDelta.m, 1);
	else 					mpz_set_ui(pDelta.m, 0);

	// Delta'=E(0) or E(1)
	paillier_enc(&cDelta, mPubKey, &pDelta, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SEQ2
	DebugDec("Delta'=E(0) or E(1)\t\t", &cDelta, idx, sDHidx);
	#endif

	//------------------------------------------------------------------------------
	// send E(delta)
    paillier_ciphertext_to_bytes(bDelta, ENC_SIZE, &cDelta);
	#ifdef _DEBUG_SEQ2
    DebugCom("(buf) Delta'=E(0) or E(1) (Hex)\t\t", bDelta, ENC_SIZE, idx, sDHidx);
	#endif

    sLen = ENC_SIZE;
	SetSendMsg(ucSendPtr, ucRecvPtr, bDelta, sLen);
    sLen += HED_SIZE;
	#ifdef _DEBUG_SEQ2
    DebugCom("(buf) Delta'=E(0) or E(1) (copied)\t", ucSendPtr, sLen, idx, sDHidx);
	#endif

	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}
	//------------------------------------------------------------------------------

	#ifdef _DEBUG_Assert
//	assert(	 pV.m->_mp_alloc >= GMP_N_SIZE);
//	assert(   pT.m->_mp_alloc == 1);
//	assert(	 cV.c->_mp_alloc == 2*GMP_N_SIZE);
//	assert(cB.c->_mp_alloc >= 1);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return 0;
}

// message processing in step 2
int PaillierCrypto::Step_2(unsigned short idx, unsigned char* ucRecvPtr, recv_t* tRecv)
{
	cTxt_t cT;
	pTxt_t pT, pV[DATA_NUM];

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, pT.m, NULL);
	for (int i=0 ; i<DATA_NUM ; i++)
		mpz_inits(pV[i].m, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE+1);
	mpz_init2(pT.m, 2*GMP_N_SIZE+1);
	for (int i=0 ; i<DATA_NUM ; i++) {
		mpz_init2(pV[i].m, 2*GMP_N_SIZE+1);
	}
	#endif

	unsigned char bW[ENC_SIZE]={0,};
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned short sLen, sDHidx;

	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	int iClusterSize, iDelta, src, dst;
	short sNum;

	//------------------------------------------------------------------------------
   sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_STEP2
	DebugCom("(buf) header of step 2 (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_STP2)&&(sLen==4));
	#endif

	iClusterSize = (int)Byte2Short(ucRecvPtr+HED_SIZE);
	sNum = Byte2Short(ucRecvPtr+HED_SIZE+2);

	SetSendMsg(ucSendPtr, ucRecvPtr, 0);			// ack message

	// memory initialization
	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;
	//-----------------------------------------------------

	// receive delta[i][j]
	for (int j=0 ; j<2 ; j++) {
		for (int i=0 ; i<iClusterSize ; i++) {
			// send ack (header + 00: 6 bytes)
			#ifdef _DEBUG_STEP2
			DebugCom("(buf) ack (Hex)\t\t", ucSendPtr, HED_SIZE+1, idx, sDHidx);
			#endif

		    // DH <-- CSP (ack)
			try {
				mCSocket->send(ucSendPtr, HED_SIZE+1);
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}
			//-----------------------------------------------------

			// DH --> CSP
			ulRecv.lock();
			tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
			ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
			ulRecv.unlock();

			// receive delta[i][j]
			sDHidx = Byte2Short(ucRecvPtr);
			sLen = Byte2Short(ucRecvPtr+HED_LEN);
			#ifdef _DEBUG_STEP2
		    DebugCom("(buf) delta[i][j] (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
			#endif

			#ifdef _DEBUG_Assert
		    assert((*(ucRecvPtr+2)==COM_STP2)&&(sLen==ENC_SIZE));
			#endif

			paillier_ciphertext_from_bytes(&cT, ucRecvPtr+HED_SIZE, ENC_SIZE);
			#ifdef _DEBUG_STEP2
			DebugOut("delta[i][j] (Copied)\t\t", cT.c, idx, sDHidx);
			DebugDec("delta[i][j]' \t\t", &cT, idx, sDHidx);
			#endif

			// ack message
			SetSendMsg(ucSendPtr, ucRecvPtr, 0);

			// memory initialization
			memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
			ucRecvPtr[1] = 0xff;
			//-----------------------------------------------------

			paillier_dec(&pT, mPubKey, mPrvKey, &cT);

			#ifdef _DEBUG_STEP2
			DebugOut("Encrypted cDelta[i,j] \t\t", cT.c, idx, sDHidx);
			DebugOut("Decrypted cDelta[i,j] \t\t", pT.m, idx, sDHidx);
			#endif

			// compare with 1
			iDelta = mpz_cmp_d(pT.m, 1);
			if ((iDelta==0)&&(j==0))
				dst = i;
			else if ((iDelta==0)&&(j==1))
				src = i;
		}
	}

	for (int j=0 ; j<sNum ; j++) {

	    // receive v[i][iSize]
		for (int i=0 ; i<iClusterSize ; i++) {
			// send ack (header + 00: 6 bytes)
			#ifdef _DEBUG_STEP2
			DebugCom("(buf) ack (Hex)\t\t", ucSendPtr, HED_SIZE+1, idx, sDHidx);
			#endif

		    // DH <-- CSP (ack)
			try {
				mCSocket->send(ucSendPtr, HED_SIZE+1);
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}
			//-----------------------------------------------------

			// DH --> CSP : v[i][iSize]
			ulRecv.lock();
			tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
			ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
			ulRecv.unlock();

			sDHidx = Byte2Short(ucRecvPtr);
			sLen = Byte2Short(ucRecvPtr+HED_LEN);
			#ifdef _DEBUG_STEP2
		    DebugCom("(buf) v[i][iSize] (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
			#endif

			#ifdef _DEBUG_Assert
		    assert((*(ucRecvPtr+2)==COM_STP2)&&(sLen==ENC_SIZE));
			#endif

			paillier_ciphertext_from_bytes(&cT, ucRecvPtr+HED_SIZE, ENC_SIZE);
			#ifdef _DEBUG_STEP2
			DebugOut("v[i][iSize] (Copied)\t\t", cT.c, idx, sDHidx);
			DebugDec("v[i][iSize] \t\t", &cT, idx, sDHidx);
			#endif

			// ack message
			SetSendMsg(ucSendPtr, ucRecvPtr, 0);

			// memory initialization
			memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
			ucRecvPtr[1] = 0xff;
			//-----------------------------------------------------

			paillier_dec(pV+i, mPubKey, mPrvKey, &cT);

			#ifdef _DEBUG_STEP2
			DebugOut("v[i][iSize] (ciphertext)\t\t", cT.c, idx, sDHidx);
			DebugOut("v[i][iSize] (plaintext)\t\t", pV[i].m, idx, sDHidx);
			#endif
		}

		//-----------------------------------------------------
		// clustering
		mpz_add(pV[dst].m, pV[dst].m, pV[src].m);
		//-----------------------------------------------------

		// send E(w[])
		for (int i=0 ; i<iClusterSize ; i++) {
			paillier_enc(&cT, mPubKey, pV+i, paillier_get_rand_devurandom);
			#ifdef _DEBUG_STEP2
			DebugOut("w[i][iSize] (ciphertext)\t", cT.c, idx, sDHidx);
			DebugDec("w[i][iSize] \t\t\t", &cT, idx, sDHidx);
			#endif

			//-----------------------------------------------------
			bW[0]=0;
			paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cT);
			#ifdef _DEBUG_STEP2
			DebugCom("(buf) w[i][iSize] (Hex)\t", bW, ENC_SIZE, idx, sDHidx);
			#endif

			sLen = ENC_SIZE;
			SetSendMsg(ucSendPtr, bW, idx, COM_STP2, ENC_SIZE);
			sLen += HED_SIZE;
			#ifdef _DEBUG_STEP2
			DebugCom("w[i][iSize] (copied)\t\t", ucSendPtr, sLen, idx, sDHidx);
			#endif

			try {
				mCSocket->send(ucSendPtr, sLen);
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}
			//-----------------------------------------------------

			// DH --> CSP (ack)
			ulRecv.lock();
			tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
			ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
			ulRecv.unlock();

			// ack (header + 00: 6 bytes)
			#ifdef _DEBUG_STEP2
		    DebugCom("(buf) ack (Hex)\t\t", ucRecvPtr, HED_SIZE+1, idx, sDHidx);
			#endif
			#ifdef _DEBUG_Assert
		    SetSendMsg(ucSendPtr, ucRecvPtr, 0);		// ack 확인용 메시지
		    assert(memcmp(ucRecvPtr, ucSendPtr, HED_SIZE+1)==0);
			#endif

			// memory initialization
			memset(ucRecvPtr, 0, HED_SIZE);
			ucRecvPtr[1] = 0xff;
		}
	}

	#ifdef _DEBUG_Assert
//	assert(	 pV.m->_mp_alloc >= GMP_N_SIZE);
//	assert(   pT.m->_mp_alloc == 1);
//	assert(	 cV.c->_mp_alloc == 2*GMP_N_SIZE);
//	assert(cB.c->_mp_alloc >= 1);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return 0;
}

// message processing in step 2
int PaillierCrypto::Step_2_Delta(unsigned short idx, unsigned char* ucRecvPtr)
{
	cTxt_t cT;
	pTxt_t pT;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, pT.m, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE+1);
	mpz_init2(pT.m, 2*GMP_N_SIZE+1);
	#endif

	unsigned short sLen, sDHidx;
	int iClusterSize, iDelta;
	unsigned char ucComDelta12;

	//------------------------------------------------------------------------------
	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_STEP2
	DebugCom("(buf) delta[n][j] (Hex)\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen%ENC_SIZE==0);
	#endif

	iClusterSize = (int)(sLen/ENC_SIZE);
	ucComDelta12 = ucRecvPtr[2];

	for (int i=0 ; i<iClusterSize ; i++) {
		paillier_ciphertext_from_bytes(&cT, ucRecvPtr+HED_SIZE+i*ENC_SIZE, ENC_SIZE);
		#ifdef _DEBUG_STEP2
		DebugOut("delta[i][j] (Copied)\t\t", cT.c, idx, sDHidx);
		DebugDec("delta[i][j]' \t\t\t", &cT, idx, sDHidx);
		#endif

		paillier_dec(&pT, mPubKey, mPrvKey, &cT);

		#ifdef _DEBUG_STEP2
		DebugOut("Decrypted cDelta[i,j] \t", pT.m, idx, sDHidx);
		#endif

		// compare with 1
		iDelta = mpz_cmp_d(pT.m, 1);
		if ((iDelta==0)&&(ucComDelta12==COM_ST2D1))
			mSt2Dst = i;
		else if ((iDelta==0)&&(ucComDelta12==COM_ST2D2))
			mSt2Src = i;
	}

	#ifdef _DEBUG_Assert
//	assert(pT.m->_mp_alloc == 1);
//	assert(cV.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+iClusterSize*ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return 0;
}

// message processing in step 2
int PaillierCrypto::Step_2_v(unsigned short idx, unsigned char* ucRecvPtr)
{
	cTxt_t cT;
	pTxt_t pV[DATA_NUM];

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, NULL);
	for (int i=0 ; i<DATA_NUM ; i++)
		mpz_inits(pV[i].m, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE+1);
	for (int i=0 ; i<DATA_NUM ; i++) {
		mpz_init2(pV[i].m, 2*GMP_N_SIZE+1);
	}
	#endif

	unsigned char bW[ENC_SIZE]={0,};
	unsigned char ucSendPtr[HED_SIZE+DATA_NUM*ENC_SIZE]={0,};
	unsigned short sLen, sDHidx;
	int iClusterSize;

	//------------------------------------------------------------------------------
   sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_STEP2
	DebugCom("(buf) header of step 2 (Hex)\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen%ENC_SIZE==0);
	#endif

	iClusterSize = (int)(sLen/ENC_SIZE);

	for (int i=0 ; i<iClusterSize ; i++) {
		paillier_ciphertext_from_bytes(&cT, ucRecvPtr+HED_SIZE+i*ENC_SIZE, ENC_SIZE);
		#ifdef _DEBUG_STEP2
		DebugOut("v[i][iSize] (Copied)\t\t", cT.c, idx, sDHidx);
		//DebugDec("v[i][iSize] \t\t", &cT, idx, sDHidx);
		#endif

		paillier_dec(pV+i, mPubKey, mPrvKey, &cT);

		#ifdef _DEBUG_STEP2
		DebugOut("v[i][iSize] (plaintext)\t", pV[i].m, idx, sDHidx);
		#endif
	}

	//-----------------------------------------------------
	// clustering
	#ifdef _DEBUG_STEP2
	DebugOut("Source:\t\t\t", mSt2Src, idx, sDHidx);
	DebugOut("Source:\t\t\t", pV[mSt2Src].m, idx, sDHidx);
	DebugOut("Destination:\t\t\t", mSt2Dst, idx, sDHidx);
	DebugOut("Destination:\t\t\t", pV[mSt2Dst].m, idx, sDHidx);
	#endif

	mpz_add(pV[mSt2Dst].m, pV[mSt2Dst].m, pV[mSt2Src].m);
	//-----------------------------------------------------

	// send E(w[])
	SetSendMsg(ucSendPtr, bW, idx, COM_ST2V, 0);
	for (int i=0 ; i<iClusterSize ; i++) {
		paillier_enc(&cT, mPubKey, pV+i, paillier_get_rand_devurandom);
		#ifdef _DEBUG_STEP2
		DebugOut("w[i][iSize] (ciphertext)\t", cT.c, idx, sDHidx);
		DebugDec("w[i][iSize] \t\t\t", &cT, idx, sDHidx);
		#endif

		//-----------------------------------------------------
		bW[0]=0;
		paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cT);
		memcpy(ucSendPtr+HED_SIZE+(i*ENC_SIZE), bW, ENC_SIZE);
	}
	ucSendPtr[3] = (unsigned char)((iClusterSize*ENC_SIZE) & 0x000000ff);
	ucSendPtr[4] = (unsigned char)(((iClusterSize*ENC_SIZE) & 0x0000ff00) >> 8);
	#ifdef _DEBUG_Step2
	DebugCom("(buf) v[n][j] (Hex)\t\t", ucSendPtr, HED_SIZE+(iClusterSize*ENC_SIZE), idx, sDHidx);
	#endif

	sLen = HED_SIZE+(iClusterSize*ENC_SIZE);
	try {
		mCSocket->send(ucSendPtr, sLen);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	#ifdef _DEBUG_Assert
//	assert(cV.c->_mp_alloc == 2*GMP_N_SIZE);
//	assert(pV.m->_mp_alloc >= GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+iClusterSize*ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return 0;
}

// termination message processing
void PaillierCrypto::TerminatePgm(unsigned short idx, unsigned char* ucRecvPtr)
{
	short sLen, sDHidx;

	sDHidx = Byte2Short(ucRecvPtr);
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_TERMINATE_PGM
    DebugCom("Terminate Program (received)\t", ucRecvPtr, sLen+HED_SIZE, idx, sDHidx);
	#endif

	#ifdef _DEBUG_Assert
	assert(sLen == 3);
	#endif

	delete[] ucRecvPtr;

	return;
}

// construct sending message
inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, const unsigned char* ucRecvPtr, unsigned char* Data, const unsigned short Len)
{
    memcpy(ucSendPtr, ucRecvPtr, HED_LEN);
	ucSendPtr[3] = (unsigned char)(Len & 0x000000ff);
	ucSendPtr[4] = (unsigned char)((Len & 0x0000ff00) >> 8);
    memcpy(ucSendPtr+HED_SIZE, Data, Len);

	return;
}

// construct sending message
inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, const unsigned char* ucRecvPtr, unsigned char bData)
{
    memcpy(ucSendPtr, ucRecvPtr, HED_LEN);
    ucSendPtr[HED_LEN]  = 1;
    ucSendPtr[HED_LEN+1]= 0;
    ucSendPtr[HED_SIZE] = bData;

	return;
}

// construct sending message
inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data,
		const unsigned short Idx, const unsigned char Tag, const unsigned short Len)
{
	ucSendPtr[0] = (unsigned char)(Idx & 0x000000ff);
	ucSendPtr[1] = (unsigned char)((Idx & 0x0000ff00) >> 8);
	ucSendPtr[2] = Tag;
	ucSendPtr[3] = (unsigned char)(Len & 0x000000ff);
	ucSendPtr[4] = (unsigned char)((Len & 0x0000ff00) >> 8);
	memcpy(ucSendPtr+HED_SIZE, Data, Len);

	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const mpz_t pData, const short idx, const short DHidx)
{
	gmp_printf("[%03d-CSP-%03d] %s : %ZX\n", idx, DHidx, pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, short pData, const short idx, const short DHidx)
{
	printf("[%03d-CSP-%03d] %s : %d\n", idx, DHidx, pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugDec(const char* pMsg, cTxt_t* cTxt, const short idx, const short DHidx)
{
	pTxt_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("[%03d-CSP-%03d] %s : %ZX\n", idx, DHidx, pMsg, &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 2*GMP_N_SIZE+1);
	#endif

    return;
}

inline void PaillierCrypto::DebugCom(const char* pMsg, unsigned char *data, int len, const short idx, const short DHidx)
{
	constexpr char hexmap[] = {'0', '1', '2', '3', '4', '5', '6', '7',
							   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	char s[len*2+1]={0,};
	for (int i = 0; i < len; ++i) {
		s[2*i]	 = hexmap[(data[i] & 0xF0) >> 4];
		s[2*i+1] = hexmap[data[i] & 0x0F];
	}
	printf("[%03d-CSP-%03d] %s : %s\n", idx, DHidx, pMsg, s);

	return;
}

void PaillierCrypto::DebugComMain(const char* pMsg, unsigned char *data, int len, const short idx, const short DHidx)
{
	constexpr char hexmap[] = {'0', '1', '2', '3', '4', '5', '6', '7',
							   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	char s[len*2+1]={0,};
	for (int i = 0; i < len; ++i) {
		s[2*i]	 = hexmap[(data[i] & 0xF0) >> 4];
		s[2*i+1] = hexmap[data[i] & 0x0F];
	}
	printf("[%03d-CSP-%03d] %s : %s\n", idx, DHidx, pMsg, s);

	return;
}

inline unsigned short PaillierCrypto::Byte2Short(unsigned char* In)
{
	return (In[1]<<8) + In[0];
}

inline void PaillierCrypto::paillier_ciphertext_from_bytes(cTxt_t* ct, void* c, int len )
{
	mpz_import(ct->c, len, 1, 1, 0, 0, c);
	return;
}

inline void PaillierCrypto::paillier_ciphertext_to_bytes(unsigned char* buf, int len, cTxt_t* ct )
{
	int cur_len;
	cur_len = mpz_sizeinbase(ct->c, 2);
	cur_len = PAILLIER_BITS_TO_BYTES(cur_len);
	mpz_export(buf + (len - cur_len), 0, 1, 1, 0, 0, ct->c);
	return;
}

int PaillierCrypto::Receiver(recv_t* tRecv)
{
	unsigned char* ucRecvBuf[THREAD_NUM+1];
	unsigned short sLen;
	short sIdx, sbx=0;
	int iReceivedLen;

	for (int i=0 ; i<THREAD_NUM+1 ; i++) {
		ucRecvBuf[i] = new unsigned char[MAXRECV];
		memset(ucRecvBuf[i], 0, MAXRECV);
		ucRecvBuf[i][1] = 0xff;
	}

	while (1) {
		while (ucRecvBuf[sbx][1]!=0xff) {
			sbx++;
			if (sbx>=THREAD_NUM+1)
				sbx=0;
		}
		iReceivedLen = 0;

		try {
			iReceivedLen = mCSocket->recv(ucRecvBuf[sbx], HED_SIZE);
			sLen = Byte2Short(ucRecvBuf[sbx]+HED_LEN);
			while (iReceivedLen != HED_SIZE+sLen)
				iReceivedLen += mCSocket->recv(ucRecvBuf[sbx]+iReceivedLen, HED_SIZE+sLen-iReceivedLen);
	    }
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		sIdx = Byte2Short(ucRecvBuf[sbx]);

		#ifdef _DEBUG_Communication
	    DebugCom("(Receiver Thread) RecvBuf (Hex)\t", ucRecvBuf[sbx], HED_SIZE+sLen, 999, sIdx);
		#endif

		//-----------------------------------------------------
		if (*(ucRecvBuf[sbx]+2) == COM_TERM) {
			assert(*(ucRecvBuf[sbx]+HED_SIZE) == 0xff);
			unsigned char ucKillMsg[HED_SIZE+3] = {0, };
			memcpy(ucKillMsg, ucRecvBuf[sbx], HED_LEN);
			ucKillMsg[HED_LEN] = 3;
			ucKillMsg[HED_SIZE+0] = (unsigned char) (MOD_SIZE & 0x000000ff);
			ucKillMsg[HED_SIZE+1] = (unsigned char)((MOD_SIZE & 0x0000ff00) >> 8);
			ucKillMsg[HED_SIZE+2] = THREAD_NUM;

			#ifdef _DEBUG_Communication
			DebugCom("(buf) KillingMsg (copied)\t", ucKillMsg, HED_SIZE+3, 999, sIdx);
			#endif

			try {
				mCSocket->send(ucKillMsg, sizeof(ucKillMsg));
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}

			for (int i=0 ; i<THREAD_NUM ; i++) {
				memcpy(ucRecvBuf[i], ucKillMsg, HED_SIZE+1);
				tRecv->m.lock();
				tRecv->pa[i] = ucRecvBuf[i];
				tRecv->m.unlock();
			}
			tRecv->cv.notify_all();

			return 0;
		}
		//-----------------------------------------------------

		tRecv->m.lock();
		tRecv->pa[sIdx] = ucRecvBuf[sbx];
		tRecv->m.unlock();
		tRecv->cv.notify_all();
	}

	return 0;
}
