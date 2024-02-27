/*
 * PaillierCrypto.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: PJS
 */


#include "PaillierCrypto.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <wchar.h>
#include <cstring>
#include "SocketException.h"


PaillierCrypto::PaillierCrypto() : mPubKey(NULL)
{
	try {
		// set Server port
		ServerSocket Server(3000);
		mSSocket = new ServerSocket();
		Server.accept(*mSSocket);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	gmp_randinit_default(state);
	mPrvKey = NULL;
}

PaillierCrypto::~PaillierCrypto()
{
	paillier_freepubkey(mPubKey);
	paillier_freeprvkey(mPrvKey);
	delete(mSSocket);
}

void PaillierCrypto::SetPubKey()
{
	char cRecvBuf[KEY_SIZE*2+1]={0,};
	int iRecvSize;

	// get public key from CSP
	try {
		iRecvSize = mSSocket->recv(cRecvBuf, KEY_SIZE*2);
    }
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	#ifdef _DEBUG_Assert
	assert(iRecvSize == KEY_SIZE*2);
	#endif

	#ifdef _DEBUG_Setting
    std::cout << "[DH] Public Key (Hex)\t\t : " << cRecvBuf << std::endl;
	#endif

	mPubKey = paillier_pubkey_from_hex(cRecvBuf);
	#ifdef _DEBUG_Setting
    DebugOut("Public Key (Copied)\t", mPubKey->n, 999);
	#endif

	return;
}

// for debugging
void PaillierCrypto::SetPrvKey()
{
	char cRecvBuf[KEY_SIZE*2+1]={0,};
	int iRecvSize;

	// get private key from CSP for debugging
	try {
		iRecvSize = mSSocket->recv(cRecvBuf, KEY_SIZE*2);
    }
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	#ifdef _DEBUG_Assert
	assert(iRecvSize == KEY_SIZE*2);
	#endif

	#ifdef _DEBUG_Setting
    std::cout << "[DH] Private Key (Hex)\t\t : " << cRecvBuf << std::endl;
	#endif

	mPrvKey = paillier_prvkey_from_hex(cRecvBuf, mPubKey);
	#ifdef _DEBUG_Setting
    DebugOut("Private Key (Copied)\t", mPrvKey->lambda, 999);
	#endif

	return ;
}

// initialize shared data
void PaillierCrypto::InitializeShdDat(sync_t *tSync, recv_t *tRecv, shdTbl_t *tTbl, shdDat_t* tDat, shdVal_t* tVal)
{
	pTxt_t 				pT;

	// thread의 sync를 맞추기 위한 변수 초기화
	// initialization the sync threads
	tSync->c1 = tSync->c2 = 0;
	for (int i=0 ; i<THREAD1_NUM ; i++) {
		tSync->d1[i] = tSync->d2[i] = 0;
		tRecv->pa[i] = NULL;
	}

	#ifdef _DEBUG_INIT_1
	mpz_init(pT.m);
	#else
	mpz_init2(pT.m, GMP_N_SIZE+1);
	#endif

	// ----- initialize Shared Tables
	#ifdef _DEBUG_INIT_1
	for (int i=0 ; i<DATA_NUM ; i++) {
		for (int j=0 ; j<DATA_NUM ; j++) {
			mpz_inits(tTbl->cCT[i][j].c, tTbl->cBeta[i][j].c, NULL);
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
				mpz_inits(tTbl->cDT3[i][j][k].c, tTbl->cDT4[i][j][k].c, NULL);
		}
		mpz_inits(tTbl->cGam2[i][0].c, tTbl->cGam2[i][1].c, tTbl->cGam1[i].c, NULL);
	}
	#else
	for (int i=0 ; i<DATA_NUM ; i++) {
		for (int j=0 ; j<DATA_NUM ; j++) {
			mpz_init2(tDat->cCT[i][j].c, 3*GMP_N_SIZE+1);
			mpz_init2(tDat->cBeta[i][j].c, 3*GMP_N_SIZE+1);
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				mpz_init2(tDat->cDT3[i][j][k].c, 3*GMP_N_SIZE+1);
				mpz_init2(tDat->cDT4[i][j][k].c, 3*GMP_N_SIZE+1);
			}
		}
		mpz_init2(tDat->cGam2[i][0].c, 3*GMP_N_SIZE+1);
		mpz_init2(tDat->cGam2[i][1].c, 3*GMP_N_SIZE+1);
		mpz_init2(tDat->cGam1[i].c, 3*GMP_N_SIZE+1);
	}
	#endif

	// ----- initialize Shared Data

	// initialization
	#ifdef _DEBUG_INIT_1
	for (int i=0 ; i<DATA_NUM ; i++)
		for (int j=0 ; j<DATA_DIM ; j++)
			mpz_inits(tDat->cD[i][j].c, tDat->_cD[i][j].c, NULL);

	// Step 2
	for (int i=0 ; i<DATA_NUM-1 ; i++)
		mpz_init(tDat->cGamS[i].c);

	// Step 4
	for (int i=0 ; i<DATA_NUM-1 ; i++) {
		for (int j=0 ; j<DATA_NUM-1 ; j++)
			mpz_init(tDat->cBs[i][j].c);
	}
	for (int i=0 ; i<DATA_NUM-2 ; i++)
		mpz_init(tDat->cGs[i].c);

	// SMIN
	for (int i=0 ; i<THREAD1_NUM ; i++) {
		for (int j=0 ; j<DTRI_SIZ ; j++)
			mpz_init(tDat->cC[i][j].c);
		for (int j=0 ; j<THREAD1_NUM ; j++)
			mpz_init(tDat->cS[i][j].c);
		mpz_inits(tDat->cQ[i].c, tDat->cSum[i].c, NULL);
	}

	// SMIN^mode
	for (int i=0 ; i<THREAD1_NUM ; i++) {
		for (int j=0 ; j<DTRI_SIZ ; j++)
			mpz_init(tDat->cU[i][j].c);
		for (int j=0 ; j<DATA_NUM-1 ; j++)
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
				mpz_init(tDat->cDmin_ij[i][j][k].c);
	}
	for (int i=0 ; i<DATA_NUM ; i++) {
		for (int j=0 ; j<DATA_NUM ; j++)
			mpz_init(tDat->cMDT[i][j].c);
		for (int j=0 ; j<DATA_SQUARE_LENGTH ; j++)
			mpz_init(tDat->cDmin_i[i][j].c);
	}

	#else
	for (int i=0 ; i<DATA_NUM ; i++)
		for (int j=0 ; j<DATA_DIM ; j++) {
			mpz_init2(tDat->cD[i][j].c, 3*GMP_N_SIZE+1);
			mpz_init2(tDat->_cD[i][j].c, 3*GMP_N_SIZE+1);
			mpz_init2(tDat->cMDT[i][j].c, 3*GMP_N_SIZE+1);
		}
	for (int i=0 ; i<THREAD1_NUM ; i++) {
		for (int j=0 ; j<THREAD1_NUM ; j++)
			mpz_init2(tDat->cS[i][j].c, 3*GMP_N_SIZE+1);
		mpz_init2(tDat->cSum[i].c, 3*GMP_N_SIZE+1);
		mpz_init2(tDat->cQ[i].c, 3*GMP_N_SIZE+1);
	}
	for (int i=0 ; i<DTRI_SIZ ; i++)
		mpz_init2(tDat->cC[i].c, 3*GMP_N_SIZE+1);
	for (int j=0 ; j<DATA_SQUARE_LENGTH ; j++) {
		for (int i=0 ; i<DATA_NUM-1 ; i++)
			mpz_init2(tDat->cDmin_ij[i][j].c, 3*GMP_N_SIZE+1);
		for (int i=0 ; i<DATA_NUM-1 ; i++)
			mpz_init2(tDat->cDmin_i[i][j].c, 3*GMP_N_SIZE+1);
	}
	for (int i=0 ; i<DATA_NUM-1 ; i++) {
		mpz_init2(tDat->cGamS[i].c, 3*GMP_N_SIZE+1);
		for (int j=0 ; j<DATA_NUM-1 ; j++)
			mpz_init2(tDat->cBs[i][j].c, 3*GMP_N_SIZE+1);
	}
	for (int i=0 ; i<DATA_NUM-2 ; i++) {
		mpz_init2(tDat->cGs[i].c, 3*GMP_N_SIZE+1);
	}
	#endif

	tDat->iNumActThds = THREAD1_NUM;

	// ----- initialize Shared Values
	#ifdef _DEBUG_INIT_1
	mpz_inits(tVal->pN1.m, tVal->c0.c, tVal->c1.c, tVal->cN1.c, NULL);
	#else
	mpz_init2(tVal->pN1.m, GMP_N_SIZE+1);
	mpz_init2(tVal->c0.c, 3*GMP_N_SIZE);
	mpz_init2(tVal->c1.c, 3*GMP_N_SIZE);
	mpz_init2(tVal->cN1.c, 3*GMP_N_SIZE);
	#endif


	// ----- compute Shared Value

	// N-1
	mpz_sub_ui(tVal->pN1.m, mPubKey->n, 1);

	#ifdef _DEBUG_Setting
	int idx=0;
	DebugOut("N\t\t\t", mPubKey->n, idx);
	DebugOut("N-1\t\t\t", tVal->pN1.m, idx);
	#endif

	// c0 = E(0)
	mpz_set_ui(pT.m, 0);													// pT=0
	paillier_enc(&(tVal->c0), mPubKey, &pT, paillier_get_rand_devurandom);

	// c1 = E(1)
	mpz_set_ui(pT.m, 1);													// pT=1
	paillier_enc(&(tVal->c1), mPubKey, &pT, paillier_get_rand_devurandom);

	#ifdef _DEBUG_Setting
	DebugOut("c1 = E(1) (ciphertext)\t\t", tVal->c1.c, idx);
	DebugDec("c1 = E(1)\t\t", &(tVal->c1), idx);
	#endif

	// cN1 = E(N-1)
	paillier_enc(&(tVal->cN1), mPubKey, &(tVal->pN1), paillier_get_rand_devurandom);


	// ----- (step 0: initialization)

	// initialize CT
	for (int i=0 ; i<DATA_NUM ; i++) {
		for (int j=0 ; j<DATA_NUM ; j++) {
			if (i==j)
				mpz_set(tTbl->cCT[i][j].c, tVal->c1.c);	// CT[i,i] = E(1)
			else
				mpz_set_ui(tTbl->cCT[i][j].c, 1);			// CT[i,j] = E(0)
		}
	}

	// initialize DT : DT[i,j,k] = E(0) (complement)
	for (int i=0 ; i<DATA_NUM ; i++)
		for (int j=0 ; j<DATA_NUM ; j++)
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				mpz_set(tTbl->cDT3[i][j][k].c, tVal->c1.c);
				mpz_set(tTbl->cDT4[i][j][k].c, tVal->c1.c);
			}

	return;
}

bool PaillierCrypto::readDataset(cTxt_t cD[DATA_NUM][DATA_DIM])
{
	std::fstream DatasetFile("dataset.dat", std::fstream::in);

	#ifdef _DEBUG_Assert
	assert(DatasetFile.is_open());
	#endif

	std::string 	stTmp, stToken;
	pTxt_t			pTmp;

	#ifdef _DEBUG_INIT_1
	mpz_init(pTmp.m);
	#else
	mpz_init2(pTmp.m, 1);
	#endif

	printf("\n************   read input dataset   ************\n\n");

	// input dataset
	for (int i=0 ; i<DATA_NUM ; i++) {
		std::getline(DatasetFile, stTmp);
		#ifdef _DEBUG_Setting
		std::cout << "[DH] Dataset:\t\t" << stTmp << std::endl;
		#endif

		std::stringstream stStream(stTmp);
		for (int j=0 ; j<DATA_DIM ; j++) {
			stStream >> stToken ;
			mpz_set_ui(pTmp.m, std::atoi(stToken.c_str()));
			paillier_enc(&cD[i][j], mPubKey, &pTmp, paillier_get_rand_devurandom);

			#ifdef _DEBUG_Setting
			std::cout << "[DH] Dataset:\t\t" << stToken << std::endl;
			gmp_printf("[DH] Dataset:\t\t%ZX\n", &pTmp);
			gmp_printf("[DH] Encrypted Dataset: %ZX\n", cD[i]+j);
			#endif
		}

		#ifdef _DEBUG_Assert
		assert(i<=DATA_NUM);
		#endif
	}

	DatasetFile.close();

	#ifdef _DEBUG_Assert
	assert(pTmp.m->_mp_alloc == 1);
	#endif

	return true;
}

// secure multiplication protocol
void PaillierCrypto::SecMul(cTxt_t* cRes, cTxt_t* cEa, cTxt_t* cEb, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+2*ENC_SIZE]={0,};
	unsigned char bA[ENC_SIZE]={0,}, bB[ENC_SIZE]={0,};

	pTxt_t  pRa, pRb, pNt;
	cTxt_t cA, cB, cT, cH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pRa.m, pRb.m, pNt.m, cA.c, cB.c, cT.c, cH.c, NULL);
	#else
	mpz_init2(pRa.m, GMP_N_SIZE);
	mpz_init2(pRb.m, GMP_N_SIZE);
	mpz_init2(pNt.m, 2*GMP_N_SIZE);
	mpz_init2(cA.c, 2*GMP_N_SIZE*2);
	mpz_init2(cB.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cH.c, 2*GMP_N_SIZE);
	#endif

	// parameter : E(a), E(b)
	#ifdef _DEBUG_SM
	DebugDec("input a\t\t\t", cEa, idx);
	DebugDec("input b\t\t\t", cEb, idx);
	#endif

	// 0 <= r_a, r_b < N  (Pick	two random numbers)
	mpz_urandomm(pRa.m, state, mPubKey->n);
	mpz_urandomm(pRb.m, state, mPubKey->n);
	#ifdef _DEBUG_SM
	DebugOut("r_a (random number 1)\t\t", pRa.m, idx);
	DebugOut("r_b (random number 2)\t\t", pRb.m, idx);
	#endif

	// E(r_a), E(r_b)
	paillier_enc(&cA, mPubKey, &pRa, paillier_get_rand_devurandom);
	paillier_enc(&cB, mPubKey, &pRb, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
	DebugOut("E(r_a)\t\t\t\t", cA.c, idx);
	DebugOut("E(r_b)\t\t\t\t", cB.c, idx);
	#endif

	// a' = E(a)*E(r_a) = E(a+r_a)
	// b' = E(b)*E(r_b) = E(b+r_b)
	paillier_mul(mPubKey, &cA, &cA, cEa);
	paillier_mul(mPubKey, &cB, &cB, cEb);
	#ifdef _DEBUG_SM
	DebugOut("a'=E(a+r_a)\t\t\t", cA.c, idx);
	DebugOut("b'=E(b+r_b)\t\t\t", cB.c, idx);
	DebugDec("a + r_a\t\t\t", &cA, idx);
	DebugDec("b + r_b\t\t\t", &cB, idx);
	#endif

	//------------------------------------------------------------------------------
	// send a' and b'
	paillier_ciphertext_to_bytes(bA, ENC_SIZE, &cA);
	paillier_ciphertext_to_bytes(bB, ENC_SIZE, &cB);
	SetSendMsg(ucSendPtr, bA, bB, idx, COM_MUL2, ENC_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
    DebugCom("(buf) a'=E(a+r_a), b'=E(b+r_b) (Hex)\t", ucSendPtr, 2*ENC_SIZE+HED_SIZE, idx);
	#endif

    // DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
		CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	//-----------------------------------------------------
	// t1 = E(a)^(N-r_b) = -r_b * E(a)
	// ==> Nt = N-r_b, A = E(a)^Nt
	mpz_sub(pNt.m, mPubKey->n, pRb.m);
	paillier_exp(mPubKey, &cA, cEa, &pNt);
	#ifdef _DEBUG_SM
	DebugOut("t1 = E(a)^(N-r_b) = -r_b * E(a)", cA.c, idx);
	DebugDec("t1 = E(a)^(N-r_b) = -r_b * E(a)", &cA, idx);
	#endif

	// t2 = E(b)^(N-r_a) = -r_a * E(b)
	// ==> Nt = N-r_a, B = E(b)^Nt
	mpz_sub(pNt.m, mPubKey->n, pRa.m);
	paillier_exp(mPubKey, &cB, cEb, &pNt);
	#ifdef _DEBUG_SM
	DebugOut("t2 = E(b)^(N-r_a) = -r_a * E(b)", cB.c, idx);
	DebugDec("t2 = E(b)^(N-r_a) = -r_a * E(b)", &cB, idx);
	#endif

	// tmp_Res1 = E(-r_b*a) * E(-r_a*b) = E(- r_b*a - r_a*b)
	// ==> tmp_Res1 = t1 * t2
    paillier_mul(mPubKey, cRes, &cA, &cB);
	#ifdef _DEBUG_SM
	DebugOut("tmp_Res1 = t1 * t2\t\t", cRes->c, idx);
	DebugDec("tmp_Res1 = t1 * t2\t\t", cRes, idx);
	#endif

	// E(r_a*r_b)^(N-1) = E(- r_a*r_b) = E( (N-r_a)*r_b )
	// ==> Nt = N-r_a, Nt = Nt * r_b,
	// ==> T = E(Nt)
	mpz_sub(pNt.m, mPubKey->n, pRa.m);
	mpz_mul(pNt.m, pNt.m, pRb.m);
	mpz_mod(pNt.m, pNt.m, mPubKey->n);
	paillier_enc(&cT, mPubKey, &pNt, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
	DebugOut("Tmp = E(r_a*r_b)^(N-1)\t\t", cT.c, idx);
	DebugDec("Tmp = E(r_a*r_b)^(N-1)\t\t", &cT, idx);
	#endif

	// tmp_Res2 =tmp_Res1 * T
	paillier_mul(mPubKey, cRes, cRes, &cT);
	#ifdef _DEBUG_SM
	DebugOut("E(-r_b*a)*E(-r_a*b)*E(-r_a*r_b)", cRes->c, idx);
	DebugDec("E(-r_b*a)*E(-r_a*b)*E(-r_a*r_b)", cRes, idx);
	#endif
	//-----------------------------------------------------

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive E(h) = E( (a+r_a)*(b+r_b) ) = E( ab + a*r_b + b*r_a + r_a*r_b )
	short sLen ;
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_SM
    DebugCom("(buf) Encrypted h' (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
    sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_MUL2)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(&cH, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
	DebugOut("E(h) (Copied)\t\t", cH.c, idx);
	#endif
	//------------------------------------------------------------------------------

	// E(ab) = E( ab + a*r_b + b*r_a + r_a*r_b ) * E( - a*r_b - b*r_a - r_a*r_b )
    paillier_mul(mPubKey, cRes, cRes, &cH);
	#ifdef _DEBUG_SM
	//DebugOut("E(a*b)", cRes->c, idx);
	DebugDec("E(a*b)\t\t\t\t", cRes, idx);
	#endif

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pRa.m->_mp_alloc ==   GMP_N_SIZE);
	assert(pRb.m->_mp_alloc ==   GMP_N_SIZE);
	assert(pNt.m->_mp_alloc == 2*GMP_N_SIZE);
	assert( cA.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cB.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cT.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cH.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

// secure multiplication protocol (square)
void PaillierCrypto::SecMul(cTxt_t* cRes, cTxt_t* cEa, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bA[ENC_SIZE]={0,};

	pTxt_t  pRa, pNt;
	cTxt_t cA, cT, cH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pRa.m, pNt.m, cA.c, cT.c, cH.c, NULL);
	#else
	mpz_init2(pRa.m, GMP_N_SIZE);
	mpz_init2(pNt.m, GMP_N_SIZE+1);
	mpz_init2(cA.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cH.c, 2*GMP_N_SIZE);
	#endif

	// parameter : E(a)
	#ifdef _DEBUG_SM
	DebugDec("input a\t\t\t", cEa, idx);
	#endif

	// 0 <= r_a < N  (Pick	a random number)
	mpz_urandomm(pRa.m, state, mPubKey->n);
	#ifdef _DEBUG_SM
	DebugOut("r (random number)\t\t\t", pRa.m, idx);
	#endif

	// E(r_a)
	paillier_enc(&cA, mPubKey, &pRa, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
	DebugOut("E(r_a)\t\t\t\t\t", cA.c, idx);
	#endif

	// a' = E(a)*E(r_a) = E(a+r_a)
	paillier_mul(mPubKey, &cA, &cA, cEa);
	#ifdef _DEBUG_SM
	DebugOut("a'=E(a+r_a)\t\t\t\t", cA.c, idx);
	DebugDec("a + r_a\t\t\t\t", &cA, idx);
	#endif

	//------------------------------------------------------------------------------
	// send a'
	paillier_ciphertext_to_bytes(bA, ENC_SIZE, &cA);
	SetSendMsg(ucSendPtr, bA, idx, COM_MUL1, ENC_SIZE);
	#ifdef _DEBUG_SM
    DebugCom("(buf) a'=E(a+r_a) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

    // DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
		CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	//-----------------------------------------------------
	// t1 = E(a)^(N-2r_a) = -2r_a * E(a)
	// ==> Nt = 2r_a (mod N), Nt = N - 2r_a, A ;= E(a)^Nt
	mpz_mul_2exp(pNt.m, pRa.m, 1);
	mpz_mod(pNt.m, pNt.m, mPubKey->n);
	mpz_sub(pNt.m, mPubKey->n, pNt.m);
	paillier_exp(mPubKey, &cA, cEa, &pNt);
	#ifdef _DEBUG_SM
	DebugOut("t1 = E(a)^(N-2r_a) = -2r_a * E(a)\t", cA.c, idx);
	DebugDec("t1 = E(a)^(N-2r_a) = -2r_a * E(a)\t", &cA, idx);
	#endif

	// t2 = E(-r_a^2) = E(N - r_a^2)
	// ==> Nt = r_a^2 (mod N), Nt = N - r_a^2
	// ==> T = E(Nt)
	mpz_powm_ui(pNt.m, pRa.m, 2, mPubKey->n);
	mpz_sub(pNt.m, mPubKey->n, pNt.m);
	paillier_enc(&cT, mPubKey, &pNt, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SM
	DebugOut("t2 = E(-r_a^2)\t\t\t\t", cT.c, idx);
	DebugDec("t2 = E(-r_a^2)\t\t\t\t", &cT, idx);
	#endif

	// tmp_Res2 = T * A
	paillier_mul(mPubKey, cRes, &cT, &cA);
	#ifdef _DEBUG_SM
	DebugOut("E(-2r_a) * E(-r_a^2)\t\t\t", cRes->c, idx);
	DebugDec("E(-2r_a) * E(-r_a^2)\t\t\t", cRes, idx);
	#endif
	//-----------------------------------------------------

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive E(h) = E( (a+r_a)^2 ) = E( a^2 + 2*a*r_a + r_a^2 )
	short sLen;
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_SM
    DebugCom("(buf) Encrypted h' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_MUL1)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(&cH, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SM
	DebugOut("Encrypted h' (Copied)\t\t\t", cH.c, idx);
	#endif
	//------------------------------------------------------------------------------

	// E(ab) = E( a^2 + 2*a*r_a + r_a^2 ) * E(-2r_a*a - r_a^2)
    paillier_mul(mPubKey, cRes, cRes, &cH);
	#ifdef _DEBUG_SM
	//DebugOut("E(a^2)", cRes->c, idx);
	DebugDec("E(a^2)\t\t\t\t\t", cRes, idx);
	#endif

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pRa.m->_mp_alloc == GMP_N_SIZE);
	assert(pNt.m->_mp_alloc == GMP_N_SIZE+1);
	assert( cA.c->_mp_alloc ==2*GMP_N_SIZE*2);
	assert( cT.c->_mp_alloc ==2*GMP_N_SIZE*2);
	assert( cH.c->_mp_alloc ==2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

// secure big-decomposition protocol
void PaillierCrypto::SBD(cTxt_t* cXb, cTxt_t* cXc, cTxt_t* cX, int iLen, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	cTxt_t cT, cZ, cT1;
	pTxt_t pL;
	int iGamma;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, cZ.c, cT1.c, pL.m, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT1.c, 2*GMP_N_SIZE*2);
	mpz_init2(pL.m, 2*GMP_N_SIZE*2);
	#endif

	// l = 2^(-1) mod N
	mpz_set_ui(pL.m, 2);
	assert(mpz_invert(pL.m, pL.m, mPubKey->n));

	while (1) {
		// T = E(x)
		mpz_set(cT.c, cX->c);
		#ifdef _DEBUG_SBD
		//DebugOut("T (copied)", cT.c, idx);
		DebugDec("T = E(x) \t", &cT, idx); // problem
		#endif

		for (int j=0 ; j<iLen ; j++) {
			#ifdef _DEBUG_SBD
			printf("\n<<<  %03d - th bit (Bit-Decomposition) >>>\n\n", j);
			#endif

			// compute the j-th LSB of E(x)
			EncryptedLSB(cXb+j, &cT, tVal, idx, tRecv, CommAmt);
			#ifdef _DEBUG_SBD
			//DebugOut(j+"-th bit \t\t\t", cXb+j->c, idx);
			DebugDec("j-th bit \t\t\t", cXb+j, idx);
			#endif

			// E(-x_i) = E(x_i)^(N-1)
			paillier_exp(mPubKey, &cT1, cXb+j, &(tVal->pN1));
			#ifdef _DEBUG_SBD
			//DebugOut("E(-x_i)\t\t\t", cT1.c, idx);
			DebugDec("E(-x_i)\t\t\t", &cT1, idx);
			#endif

			// E(~x_i) = E(1) * E(-x_i) : 추가 계산
			paillier_mul(mPubKey, cXc+j, &(tVal->c1), &cT1);
			#ifdef _DEBUG_SBD
			//DebugOut("E(~x_i) = E(1) * E(-x_i) \t\t", cXc[j].c, idx);
			DebugDec("E(~x_i) = E(1) * E(-x_i) \t\t", cXc+j, idx);
			#endif

			// Z = T * E(-x_i)
			paillier_mul(mPubKey, &cZ, &cT, &cT1);
			#ifdef _DEBUG_SBD
			//DebugOut("Z = T * E(-x_i)\t\t", cZ.c, idx);
			DebugDec("Z = T * E(-x_i)\t\t", &cZ, idx);
			#endif

			// T = Z^l (right shift)
			paillier_exp(mPubKey, &cT, &cZ, &pL);
			#ifdef _DEBUG_SBD
			//DebugOut("T = Z^l (right shift)\t\t", T.c, idx);
			DebugDec("T = Z^l (right shift)\t\t", &cT, idx);
			#endif
		}

		#ifdef _DEBUG_SBD
		// parameter : E(x), < E(x_0), E(x_1), ... , E(x_m-1) >
		DebugDec("E(X) \t\t\t\t", cX, idx);
		for (int j=iLen-1 ; j>=0 ; j--) {
			DebugDecBit(cXb+j);
			if (j%4 == 0)
				printf(" ");
		}
		printf("\n");
		#endif

		iGamma = SVR(cX, cXb, iLen, tVal, idx, tRecv, CommAmt);

		if (iGamma == 1)
			break;

		printf("\n<<<  ***  FAIL to SVR !!!  *** >>>\n\n");
	}
}

// subprotocol in secure big-decomposition protocol
void PaillierCrypto::EncryptedLSB(cTxt_t* cRes, cTxt_t* cT, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bY[ENC_SIZE]={0,};

	pTxt_t  pT;
	cTxt_t cRY;
	int iOdd;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pT.m, cRY.c, NULL);
	#else
	mpz_init2(pT.m, 2*GMP_N_SIZE);
	mpz_init2(cRY.c, 2*GMP_N_SIZE*2);
	#endif

	// parameter : Enrypted T
	#ifdef _DEBUG_LSB
	//DebugOut("input T\t\t\t", cT->c, idx);
	DebugDec("input T\t\t\t", cT, idx);
	#endif

	// 0 <= r < N  (Pick random number)
	mpz_urandomm(pT.m, state, mPubKey->n); 				// pT = pR
	#ifdef _DEBUG_LSB
	DebugOut("r (random number)\t\t", pT.m, idx);
	#endif

	// E(r)
	paillier_enc(&cRY, mPubKey, &pT, paillier_get_rand_devurandom);
	#ifdef _DEBUG_LSB
	//DebugOut("E(r)\t\t\t\t", cRY.c, idx);
	DebugDec("r\t\t\t\t", &cRY, idx);
	#endif

	// Y = T * E(r)
	paillier_mul(mPubKey, &cRY, cT, &cRY);
	#ifdef _DEBUG_LSB
	//DebugOut("Y = T * E(r)\t\t\t", cRY.c, idx);
	DebugDec("y = x + r\t\t\t", &cRY, idx);
	#endif

	//------------------------------------------------------------------------------
	// send Y
	paillier_ciphertext_to_bytes(bY, ENC_SIZE, &cRY);
	SetSendMsg(ucSendPtr, bY, idx, COM_LSB, ENC_SIZE);
	#ifdef _DEBUG_LSB
    DebugCom("(buf) Y = T * E(r) (Hex):\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

    // DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
		CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive E(alpha)
	short sLen ;
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_LSB
	DebugCom("(buf) E(alpha) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_LSB)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(cRes, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_LSB
	//DebugOut("E(alpha) (Copied)", cAlpha.c, idx);
	DebugDec("E(alpha) \t\t", cRes, idx);
	#endif
	//------------------------------------------------------------------------------

	iOdd = mpz_odd_p(pT.m);
	#ifdef _DEBUG_LSB
	DebugOut("r is odd?\t\t\t", iOdd, idx);
	#endif

	// (r:even) E(x_i) = alpha
	// (r:odd)  E(x_i) = E(1) * alpha^(N-1)
	if (iOdd) {
		// alpha^(N-1)
		paillier_exp(mPubKey, &cRY, cRes, &(tVal->pN1));
		#ifdef _DEBUG_LSB
		//DebugOut("alpha^(N-1)", cRY.c, idx);
		DebugDec("alpha^(N-1)\t\t\t", &cRY, idx);
		#endif

		// E(x_i) = E(1) * alpha^(N-1)
		paillier_mul(mPubKey, cRes, &(tVal->c1), &cRY);
		#ifdef _DEBUG_LSB
		//DebugOut("E(x_i) = E(1) * alpha^(N-1)", cRes->c, idx);
		DebugDec("E(x_i) = E(1) * alpha^(N-1)\t", cRes, idx);
		#endif
	}

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert( pT.m->_mp_alloc <= 2*GMP_N_SIZE);
	assert(cRY.c->_mp_alloc == 2*GMP_N_SIZE*2);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;						// 2N*2
}

// subprotocol in secure big-decomposition protocol
int PaillierCrypto::SVR(cTxt_t* cEX, cTxt_t* cEXi, int len, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bW[ENC_SIZE]={0,};

	pTxt_t pT;
	cTxt_t cUVW, cT;
	int iGamma;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pT.m, cUVW.c, cT.c, NULL);
	#else
	mpz_init2(pT.m, GMP_N_SIZE+1);
	mpz_init2(cUVW.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE);Ma
	#endif

	// parameter : E(x)
	#ifdef _DEBUG_SVR
	DebugOut("input E(x)\t\t\t", cEX->c, idx);
	DebugDec("input x\t\t\t", cEX, idx);
	#endif

	#ifdef _DEBUG_SVR
	// parameter : < E(x_0), E(x_1), ... , E(x_m-1) >
    printf("input < x_0, x_1, ..., x_m-1 >\n");
    for (int j=len-1 ; j>=0 ; j--) {
		DebugDecBit(&cEXi[j]);
		if (j%8==0)		printf("\t");
		if (j%64==0)	printf("\n");
	}
	#endif

	// U = SUM_i=0~m-1 (E(x_i))^(2^i)
    mpz_set_ui(pT.m, 1);
    mpz_set_ui(cUVW.c, 1);
	for (int i=0 ; i<len ; i++) {
		paillier_exp(mPubKey, &cT, &cEXi[i], &pT);
		paillier_mul(mPubKey, &cUVW, &cUVW, &cT);
		#ifdef _DEBUG_SVR
		//DebugOut("U = SUM_i=0~m-1 (E(x_i))^(2^i)", cUVW.c, idx);
		DebugDec("U = SUM_i=0~m-1 (E(x_i))^(2^i)", &cUVW, idx);
		#endif

		mpz_mul_2exp(pT.m, pT.m, 1);
	}

	// N-1
	#ifdef _DEBUG_SVR
	DebugOut("N\t\t\t\t", mPubKey->n, idx);	DebugOut("N-1\t\t\t\t", tVal->pN1.m, idx);
	#endif

	// V = U * E(x)^(N-1)
	paillier_exp(mPubKey, &cT, cEX, &(tVal->pN1));
	paillier_mul(mPubKey, &cUVW, &cUVW, &cT);				// 여기부터 cUVW = cV
	#ifdef _DEBUG_SVR
	//DebugOut("V = U * E(-x)", cUVW.c, idx);
	DebugDec("V = U * E(-x)\t\t\t", &cUVW, idx);
	#endif

	// 0 <= r < N  (Pick random number)
	mpz_urandomm(pT.m, state, mPubKey->n);					// 여기부터 pT = r
	#ifdef _DEBUG_SVR
	DebugOut("r (random number)\t\t", pT.m, idx);
	#endif

	// W = V^r
	paillier_exp(mPubKey, &cUVW, &cUVW, &pT);				// 여기부터 cUVW = cW
	#ifdef _DEBUG_SVR
	//DebugOut("W = V^r", cUVW.c, idx);
	DebugDec("W = V^r\t\t\t", &cUVW, idx);
	#endif

	//------------------------------------------------------------------------------
    // send E(w)
	paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cUVW);
	SetSendMsg(ucSendPtr, bW, idx, COM_SVR, ENC_SIZE);
	#ifdef _DEBUG_SVR
    DebugCom("(buf) W (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

    // DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
		CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive ooo
	short sLen ;
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_SVR
	DebugCom("(buf) Gamma' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	// a' = E(a+r_a)
	#ifdef _DEBUG_Assert
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_SVR)&&(sLen==1));
	#endif
	//------------------------------------------------------------------------------

    iGamma = *(ucRecvPtr+HED_SIZE);

	// E( (a+r_a)^2 ) = E( a^2 + 2*a*r_a + r_a^2 )
	#ifdef _DEBUG_SVR
	printf("[DH-%03d] Gamma \t\t\t\t : %d\n", idx, iGamma);
	#endif

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(  pT.m->_mp_alloc ==   GMP_N_SIZE);
	assert(cUVW.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(  cT.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return iGamma;
}

// compute the task amount of each thread (subprotocol of assignTask)
inline void PaillierCrypto::ComputeTaskAmount(int *SizeEachTask, int SizeTotalTask)
{
	int tmp, rem;

	rem = SizeTotalTask % THREAD1_NUM ;

	if (rem == 0)	{
		tmp = SizeTotalTask / THREAD1_NUM ;
		for (int i=0 ; i<THREAD1_NUM ; i++)
			SizeEachTask[i] = tmp ;						// assign tasks to threads equally
	}
	else {
		tmp = SizeTotalTask / THREAD1_NUM ;
		for (int i=0 ; i<rem ; i++)
			SizeEachTask[i] = tmp+1 ;					// task amount of thread : T+1
		for (int i=rem ; i<THREAD1_NUM ; i++)
			SizeEachTask[i] = tmp ;						// task amount of thread : T
	}

	return;
}

// compute the task amount of each thread (subprotocol of assignTask)
inline void PaillierCrypto::ComputeTaskAmount(int *SizeEachTask, int thSize, int SizeTotalTask)
{
	int tmp, rem;

	rem = SizeTotalTask % thSize ;

	if (rem == 0)	{
		tmp = SizeTotalTask / thSize ;
		for (int i=0 ; i<thSize ; i++)
			SizeEachTask[i] = tmp ;				// assign tasks to threads equally
	}
	else {
		tmp = SizeTotalTask / thSize ;
		for (int i=0 ; i<rem ; i++)
			SizeEachTask[i] = tmp+1 ;			// task amount of thread : T+1
		for (int i=rem ; i<thSize ; i++)
			SizeEachTask[i] = tmp ;				// task amount of thread : T
	}

	return;
}

// assign tasks to each thread
void PaillierCrypto::assignTask(shdDat_t *tDat)
{
	int SizeTotalTask, tmp, rem, CntTask, CntThd, thNum ;
	int SizeEachTask[THREAD1_NUM];

	memset(tDat->thSrtPos, -1, sizeof(int)*(NUM_THD_TASK*(DATA_NUM+1)*(THREAD1_NUM+1)*2));

	// ----------------------------------------------------------------------------------------
	// [POS_TRI][size][#th][i,j]	upper triangular 	task position

	for (int size=DATA_NUM ; size>0 ; size--) {
		SizeTotalTask = size * (size-1) / 2 ;								// size of total tasks

		if (SizeTotalTask >= THREAD1_NUM) {

			ComputeTaskAmount(SizeEachTask, SizeTotalTask);

			// compute the starting position of each thread's task
			tDat->thSrtPos[POS_TRI][size][0][0] = 0 ;						// the starting position of the first thread
			tDat->thSrtPos[POS_TRI][size][0][1] = 1 ;
			CntTask=0; 	CntThd=0;
			for (int i=0 ; i<size-1 ; i++) {
				for (int j=i+1 ; j<size ; j++) {
					if (CntTask >= SizeEachTask[CntThd]) {
						CntThd++;
						tDat->thSrtPos[POS_TRI][size][CntThd][0] = i ;
						tDat->thSrtPos[POS_TRI][size][CntThd][1] = j ;
						CntTask=0;
					}
					CntTask++;
				}
			}

			tDat->thSrtPos[POS_TRI][size][THREAD1_NUM][0] = size ;		// the ending position of the last thread
			tDat->thSrtPos[POS_TRI][size][THREAD1_NUM][1] = size ;
		}
		else {
			CntThd=0;
			for (int i=0 ; i<size-1 ; i++) {
				for (int j=i+1 ; j<size ; j++) {
					tDat->thSrtPos[POS_TRI][size][CntThd][0] = i ;
					tDat->thSrtPos[POS_TRI][size][CntThd][1] = j ;
					CntThd++;
				}
			}

			tDat->thSrtPos[POS_TRI][size][SizeTotalTask][0] = size ;	// the ending position of the last thread
			tDat->thSrtPos[POS_TRI][size][SizeTotalTask][1] = size ;
		}
	}

	// ----------------------------------------------------------------------------------------
	// [POS_SQU][size][#th][i,j]	square task position

	for (int size=DATA_NUM ; size>0 ; size--) {
		SizeTotalTask = size * size ;										// size of total tasks

		if (SizeTotalTask >= THREAD1_NUM) {

			ComputeTaskAmount(SizeEachTask, SizeTotalTask);

			// compute the starting position of each thread's task
			tDat->thSrtPos[POS_SQU][size][0][0] = 0 ;						// the starting position of the first thread
			tDat->thSrtPos[POS_SQU][size][0][1] = 0 ;
			tmp = 0;
			for (int i=1 ; i<THREAD1_NUM ; i++) {
				tmp += SizeEachTask[i-1];
				tDat->thSrtPos[POS_SQU][size][i][0] = tmp / size;
				tDat->thSrtPos[POS_SQU][size][i][1] = tmp % size;
			}
			tDat->thSrtPos[POS_SQU][size][THREAD1_NUM][0] = size ;		// the ending position of the last thread
			tDat->thSrtPos[POS_SQU][size][THREAD1_NUM][1] = size ;
		}
		else {
			CntThd=0;
			for (int i=0 ; i<size ; i++) {
				for (int j=0 ; j<size ; j++) {
					tDat->thSrtPos[POS_SQU][size][CntThd][0] = i ;
					tDat->thSrtPos[POS_SQU][size][CntThd][1] = j ;
					CntThd++;
				}
			}
			tDat->thSrtPos[POS_SQU][size][SizeTotalTask][0] = size ;	// the ending position of the last thread
			tDat->thSrtPos[POS_SQU][size][SizeTotalTask][1] = size ;
		}
	}

	// ----------------------------------------------------------------------------------------
	// [LIN][size][#th][0,j]	(linear) n*(n-1)/2 task position

	for (int size=DATA_NUM ; size>0 ; size--) {
		SizeTotalTask = size * (size-1) / 2  ;									// size of total tasks

		if (SizeTotalTask >= THREAD1_NUM) {
			ComputeTaskAmount(SizeEachTask, SizeTotalTask);

			// compute the starting position of each thread's task
			tDat->thSrtPos[POS_LIN][size][0][0] = 0 ;							// the starting position of the first thread
			for (int i=1 ; i<THREAD1_NUM ; i++) {
				tDat->thSrtPos[POS_LIN][size][i][0] = tDat->thSrtPos[POS_LIN][size][i-1][0] + SizeEachTask[i-1];
			}
			tDat->thSrtPos[POS_LIN][size][THREAD1_NUM][0] = SizeTotalTask ;// the ending position of the last thread
		}
		else {
			for (int i=0 ; i<SizeTotalTask ; i++) {
				tDat->thSrtPos[POS_LIN][size][i][0] = i;
			}
			tDat->thSrtPos[POS_LIN][size][SizeTotalTask][0] = SizeTotalTask ;// the ending position of the last thread
		}
	}

	// ----------------------------------------------------------------------------------------
	// [LIN][size][#th][1,j]	(linear) n task position

	for (int size=DATA_NUM ; size>0 ; size--) {
		SizeTotalTask = size  ;													// size of total tasks

		if (SizeTotalTask >= THREAD1_NUM) {
			ComputeTaskAmount(SizeEachTask, SizeTotalTask);

			// compute the starting position of each thread's task
			tDat->thSrtPos[POS_LIN][size][0][1] = 0 ;							// the starting position of the first thread
			for (int i=1 ; i<THREAD1_NUM ; i++) {
				tDat->thSrtPos[POS_LIN][size][i][1] = tDat->thSrtPos[POS_LIN][size][i-1][1] + SizeEachTask[i-1];
			}
			tDat->thSrtPos[POS_LIN][size][THREAD1_NUM][1] = SizeTotalTask ;// the ending position of the last thread
		}
		else {
			for (int i=0 ; i<SizeTotalTask ; i++) {
				tDat->thSrtPos[POS_LIN][size][i][1] = i;
			}
			tDat->thSrtPos[POS_LIN][size][SizeTotalTask][1] = SizeTotalTask ;// the ending position of the last thread
		}
	}

	// ----------------------------------------------------------------------------------------
	// [ST3][size][#th][i,j]	(Step 3) n row's task position

	memset(tDat->thNumS3RRow, -1, sizeof(int)*(DATA_NUM+1)*THREAD1_NUM);

	for (int size=DATA_NUM ; size>0 ; size--) {
		SizeTotalTask = size  ;											// size of total tasks

		tmp = SizeTotalTask / THREAD1_NUM ;
		rem = SizeTotalTask % THREAD1_NUM ;

		// compute the row (task) number in SMIN^PSV execution (single-thread)
		if (SizeTotalTask >= THREAD1_NUM) {
			tDat->thSrtPos[POS_ST3][size][0][0] = 0 ;					// the starting position of the first thread
			for (int i=1 ; i<=THREAD1_NUM ; i++)
				tDat->thSrtPos[POS_ST3][size][i][0] = tDat->thSrtPos[POS_ST3][size][i-1][0] + tmp;
		}

		// compute the row (task) number in SMIN^PSV execution (multi-thread)
		if (rem != 0) {
			thNum = size>THREAD1_NUM ? THREAD1_NUM : size ;
			for (int i=0 ; i<rem ; i++)
				tDat->thNumS3RRow[size][i] = 0;
			for (int i=0 ; i<thNum ; i++) {
				tDat->thNumS3RRow[size][i%rem]++;
			}
			for (int i=0 ; i<THREAD1_NUM ; i++)
				tDat->thSrtPos[POS_ST3][size][i][1] = (tmp*THREAD1_NUM) + (i%rem);
		}
	}

	// ----------------------------------------------------------------------------------------
	// [size][thSize][#th]	(step 3) n task position

	memset(tDat->thPosS3RRow, -1, sizeof(int)*DATA_NUM*(THREAD1_NUM+1)*(THREAD1_NUM+1));

	for (int size=DATA_NUM-1 ; size>0 ; size--) {
		SizeTotalTask = size  ;											// size of total tasks

		for (int thSize=1 ; thSize<=THREAD1_NUM ; thSize++) {
			if (SizeTotalTask >= thSize) {
				ComputeTaskAmount(SizeEachTask, thSize, SizeTotalTask);

				// compute the starting position of each thread's task
				tDat->thPosS3RRow[size][thSize][0] = 0 ;				// the starting position of the first thread
				for (int i=1 ; i<thSize ; i++) {
					tDat->thPosS3RRow[size][thSize][i] = tDat->thPosS3RRow[size][thSize][i-1] + SizeEachTask[i-1];
				}
				tDat->thPosS3RRow[size][thSize][thSize] = SizeTotalTask ;	// the ending position of the last thread
			}
			else {
				for (int i=0 ; i<SizeTotalTask ; i++) {
					tDat->thPosS3RRow[size][thSize][i] = i;
				}
			}
		}
	}
	return;
}


void PaillierCrypto::computeDTpos(shdDat_t *tDat)
{
	int cnt;

	memset(tDat->DTposST1, -1, sizeof(int)*(DATA_NUM+1)*(DTRI_SIZ+1)*2);
	// ----------------------------------------------------------------------------------------
	// (step 1: SMIN^ASB) data index of SMIN -> coordinate of DT

	for (int size = 2 ; size <= DATA_NUM ; size++) {
		cnt = 0;
		for (int i = 0 ; i < size ; i++) {
			for (int j=i+1 ; j < size ; j++) {
				tDat->DTposST1[size][cnt][0] = i;
				tDat->DTposST1[size][cnt][1] = j;
				cnt++;
			}
		}
	}

	memset(tDat->DTposST3, -1, sizeof(int)*(DATA_NUM+1)*DATA_NUM*DATA_NUM);
	// ----------------------------------------------------------------------------------------
	// (step 3: SMIN^PSV) data coordinate of SMIN -> coordinate of DT

	for (int size = 2 ; size <= DATA_NUM ; size++) {
		for (int row = 0 ; row < size ; row++) {
			cnt = 0;
			for (int j=0 ; j < size ; j++) {
				if (j == row)				continue;
				tDat->DTposST3[size][row][cnt] = j;
				cnt++;
			}
		}
	}

	return;
}

// initialization step
void PaillierCrypto::initialization(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned long uiCommAmt[2]={0, 0};

	cTxt_t cW, cX, cDist, cT[DATA_SQUARE_LENGTH];
	bool	bStart=true, bExit = false;
	int Start_i, Start_j, End_i, End_j;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cW.c, cX.c, cDist.c, NULL);
	for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
		mpz_init(cT[k].c);
	#else
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(cX.c, 2*GMP_N_SIZE*2);
	mpz_init2(cDist.c, 2*GMP_N_SIZE*2);
	for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
		mpz_init2(cT[k].c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);
	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	(*c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c1, idx);
	#endif
	if (*c1 < THREAD1_NUM)
		cv->wait(ulSync, [&] {return *c1>=THREAD1_NUM;});
	else {
		// --------------------------------------------------------------------
		printf("\n************   [initialization] compute DT   ************\n\n");
		tDat->start = time(NULL);

		// --------------------------------------------------------------------
		tTbl->n = DATA_NUM;		// (initialization) cluster size = data size

		// M-1
		#ifdef _DEBUG_Initialization
		DebugOut("N   \t", mPubKey->n, idx);
		DebugOut("N-1 \t", tVal->pN1.m, idx);
		#endif

		// p[j][k]^(M-1)  for k=1~m
		for (int j=1 ; j<DATA_NUM ; j++) {
			for (int k=0 ; k<DATA_DIM ; k++) {
				#ifdef _DEBUG_Initialization
				DebugDec("p[j][k] \t", tDat->cD[j]+k, idx);
				#endif

				paillier_exp(mPubKey, tDat->_cD[j]+k, tDat->cD[j]+k, &(tVal->pN1));

				#ifdef _DEBUG_Initialization
				//DebugOut("p[j][k]^(M-1) \t", tDat->_cD[j]+k, idx);
				DebugDec("p[j][k]^(M-1) \t", tDat->_cD[j]+k, idx);
				#endif
			}
		}
		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
	}
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	Start_i = tDat->thSrtPos[POS_TRI][DATA_NUM][idx][0] ;	Start_j = tDat->thSrtPos[POS_TRI][DATA_NUM][idx][1] ;
	End_i 	 = tDat->thSrtPos[POS_TRI][DATA_NUM][idx+1][0] ;	End_j 	 = tDat->thSrtPos[POS_TRI][DATA_NUM][idx+1][1] ;
	for (int i=0 ; i<DATA_NUM-1 ; i++) {
		if (bExit)			break;
		for (int j=i+1 ; j<DATA_NUM ; j++) {
			if (bStart) {
				i = Start_i ;		j = Start_j ;
				bStart = false;
			}
			if ((i>=End_i) && (j>=End_j)) {
				bExit = true;
				break;
			}

			// initialize Dist = E(0)
			mpz_set_ui(cDist.c, 1);
			for (int k=0 ; k<DATA_DIM ; k++) {

				// w[i,j,k] = p[i,k] * -p[j,k]
				paillier_mul(mPubKey, &cW, tDat->cD[i]+k, tDat->_cD[j]+k);
				#ifdef _DEBUG_Initialization
				//DebugOut("w[i,j,k] = p[i,k] * -p[j,k] \t\t", cW.c, idx);
				DebugDec("w[i,j,k] = p[i,k] * -p[j,k] \t\t", &cW, idx);
				#endif

				// x[i,j,k] = w[i,j,k]^2
				SecMul(&cX, &cW, idx, tRecv, uiCommAmt);
				#ifdef _DEBUG_Initialization
				//DebugOut("x[i,j,k] = w[i,j,k]^^2 \t\t\t", cX.c, idx);
				DebugDec("x[i,j,k] = w[i,j,k]^2 \t\t\t", &cX, idx);
				#endif

				// Dist[i,j] = Sum_(k=1~m) x[i,j,k]
				paillier_mul(mPubKey, &cDist, &cDist, &cX);
				#ifdef _DEBUG_Initialization
				//DebugOut("Dist[i,j] = Sum_(k=1~m) x[i,j,k] \t\t", cDist.c, idx);
				DebugDec("Dist[i,j] = Sum_(k=1~m) x[i,j,k] \t\t", &cDist, idx);
				#endif
			}

			SBD(cT, tTbl->cDT4[i][j], &cDist, DATA_SQUARE_LENGTH, tVal, idx, tRecv, uiCommAmt);

			#ifdef _DEBUG_Initialization
			for (int k=DATA_SQUARE_LENGTH-1 ; k>=0 ; k--) {
				DebugDecBit(tTbl->cDT4[i][j]+k);
				if (k%8==0)		printf("\t");
				if (k%64==0)		printf("\n");
			}
			#endif
		}
	}

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	tDat->uiCommAmtIni[0] += uiCommAmt[0];
	tDat->uiCommAmtIni[1] += uiCommAmt[1];
    (*c2)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c2, idx);
	#endif
    if (*c2 < THREAD1_NUM)
    	cv->wait(ulSync, [&] {return *c2>=THREAD1_NUM;});
    else {
    	tDat->end = time(NULL);
    	tDat->dTimeEtc[TIME_INIT][0] = (double) (tDat->end - tDat->start);
    	printf("[DH-%03d] Time \t : %.2f (seconds) \n", idx, tDat->dTimeEtc[TIME_INIT][0]);
    	printf("[DH-%03d] Comm Amount  : (DH) %lu \t(CSP) %lu \t(bytes) \n", idx, tDat->uiCommAmtIni[0], tDat->uiCommAmtIni[1]);
    	tDat->dTimeEtc[TIME_INIT][1] = -1;
		// --------------------------------------------------------------------
		#ifdef _DEBUG_Detail
		for (int i=0 ; i<DATA_NUM-1 ; i++) {
			for (int j=i+1 ; j<DATA_NUM ; j++) {
				DebugDec(tTbl->cDT4[i][j], idx);
			    for (int k=DATA_SQUARE_LENGTH-1 ; k>=0 ; k--) {
					DebugDecBit(tTbl->cDT4[i][j]+k);
					if (k%8==0)	printf("\t");
					if (k%64==0)	printf("\n");
				}
			}
		}
		#endif

		//DebugDT("Initialization", tTbl, 0, 4, idx);
		// --------------------------------------------------------------------
		*c1=0;		cv->notify_all();
    }
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	return;
}

// step 1
void PaillierCrypto::step_1(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned long uiCommAmt[2]={0, 0};
	cTxt_t *cMDT_i, (*cMDT)[DATA_NUM], *cBeta_i, *cGam2_i, *cGam1, *cGam2 ;
	int iNumActThds = tDat->iNumActThds;

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);
	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	(*c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c1, idx);
	#endif
	if (*c1 < iNumActThds)
		cv->wait(ulSync, [&] {return *c1>=iNumActThds;});
	else {
		// --------------------------------------------------------------------
		printf("\n************   [Step 1] compute MT   ************\n\n");
		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
		tDat->start = time(NULL);
	}
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------
	iNumActThds = tDat->iNumActThds;

	SmaxSmin_ASB(tDat->cMDT, tTbl, tDat, tVal, idx, tRecv, tSync, uiCommAmt);

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	tDat->uiCommAmt[DATA_NUM - tTbl->n][0][0] += uiCommAmt[0];
	tDat->uiCommAmt[DATA_NUM - tTbl->n][0][1] += uiCommAmt[1];
    (*c2)++;
    if (*c2 < iNumActThds)
    	cv->wait(ulSync, [&] {return *c2>=iNumActThds;});
    else {
    	int iClusterSize = tTbl->n;
		// --------------------------------------------------------------------
    	// Results of SMIN_ASB
		#ifdef _DEBUG_Detail
		for (int i=0 ; i<iClusterSize-1 ; i++) {
			for (int j=i+1 ; j<iClusterSize ; j++)
				DebugDec("cMDT[i][j]: \t\t", tDat->cMDT[i]+j, idx);
		}
		#endif
		// --------------------------------------------------------------------

		// single-thread process due to small task
    	cTxt_t cR[iClusterSize], cC[iClusterSize] ;

    	for (int i=0 ; i<iClusterSize ; i++)
    		mpz_inits(cR[i].c, cC[i].c, NULL);

    	// initialization
    	for (int i=0 ; i<iClusterSize ; i++) {
    		mpz_set_ui(cR[i].c, 1);
    		mpz_set_ui(cC[i].c, 1);
    	}

    	// compute row and column
    	for (int i=0 ; i<iClusterSize-1 ; i++) {
    		cMDT_i = tDat->cMDT[i];
    		for (int k=i+1 ; k<iClusterSize ; k++) {
    			paillier_mul(mPubKey, cR+i, cR+i, cMDT_i+k);
    		}
    	}

    	// r[n] = E(0), c[0] = E(0)
    	mpz_set(cR[iClusterSize-1].c, tVal->c0.c);
    	mpz_set(cC[0].c, tVal->c0.c);

    	cMDT = tDat->cMDT;
    	for (int j=1 ; j<iClusterSize ; j++) {
    		for (int k=0 ; k<j ; k++) {
    			paillier_mul(mPubKey, cC+j, cC+j, cMDT[k]+j);
    		}
    	}

    	#ifdef _DEBUG_Step1
    	for (int i=0 ; i<iClusterSize ; i++) {
    		//DebugOut("cR[i] (ciphertext): \t\t", cR[i].c, idx);
    		DebugDec("cR[i]: \t\t", cR+i, idx);
    	}
    	for (int j=0 ; j<iClusterSize ; j++) {
    		//DebugOut("cC[j] (ciphertext): \t\t", cC[j].c, idx);
    		DebugDec("cC[j]: \t\t", cC+j, idx);
    	}
    	#endif

		cGam1 = tTbl->cGam1;
    	// compute beta-table, gamma-table, and gamma-array
    	for (int i=0 ; i<iClusterSize ; i++) {
    		cBeta_i = tTbl->cBeta[i];
    		cGam2_i = tTbl->cGam2[i];
    		for (int j=0 ; j<iClusterSize ; j++) {
    			if (i == j) {
    				mpz_set(cBeta_i[i].c, cC[i].c);
    				mpz_set(cGam2_i[0].c, cR[i].c);
    				mpz_set(cGam2_i[1].c, cC[i].c);
    				paillier_mul(mPubKey, cGam1+i, cR+i, cC+i);
    			}
    			else {
    				paillier_mul(mPubKey, cBeta_i+j, cC+i, cC+j);
    			}
    		}
    	}

		// --------------------------------------------------------------------
    	tDat->end = time(NULL);
    	tDat->dTimeRnd[DATA_NUM - tTbl->n][0] = (double) (tDat->end - tDat->start);
    	printf("[DH-%03d] Time \t     : %.2f (seconds) \n", idx, tDat->dTimeRnd[DATA_NUM - tTbl->n][0]);
    	printf("[DH-%03d] Comm Amount : (DH) %lu \t(CSP) %lu \t(bytes) \n", idx, tDat->uiCommAmt[DATA_NUM - tTbl->n][0][0], tDat->uiCommAmt[DATA_NUM - tTbl->n][0][1]);

    	// print computation results
    	#ifdef _DEBUG_Step1
    	for (int i=0 ; i<iClusterSize ; i++) {
    		cBeta_i = tTbl->cBeta[i];
    		for (int j=0 ; j<iClusterSize ; j++) {
    			DebugOut("cBeta[i][j] (ciphertext): \t", cBeta_i[j].c, idx);
    			DebugDec("cBeta[i][j]: \t\t\t", cBeta_i+j, idx);
    		}
    	}
		//cGam2 = &(tTbl->cGam2[0][0]);
    	for (int j=0 ; j<2 ; j++) {
    		for (int i=0 ; i<iClusterSize ; i++) {
    			//DebugOut("cGam2[i][2] (ciphertext): \t", cGam2[i][j]->c, idx);
    			//DebugDec("cGam2[i][2]: \t\t\t", cGam2[i]+j, idx);
    			DebugOut("cGam2[i][2] (ciphertext): \t", tTbl->cGam2[i][j].c, idx);
    			DebugDec("cGam2[i][2]: \t\t\t", tTbl->cGam2[i]+j, idx);
    		}
    	}
		cGam1 = tTbl->cGam1;
    	for (int i=0 ; i<iClusterSize ; i++) {
    		DebugOut("cGam1[i] (ciphertext): \t\t", cGam1[i].c, idx);
    		DebugDec("cGam1[i]: \t\t\t\t", cGam1+i, idx);
    	}
    	#endif

    	//DebugMT("Step 1", tTbl, idx);
		// --------------------------------------------------------------------
    	*c1=0;    	cv->notify_all();
    }
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	return ;
}

// step 2
void PaillierCrypto::step_2(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync)
{
	int iClusterSize = tTbl->n;
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE*iClusterSize]={0,};
	unsigned char bV[ENC_SIZE]={0,};
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned long uiCommAmt[2]={0, 0};
	int iNumActThds = tDat->iNumActThds;

	int Start = tDat->thSrtPos[POS_LIN][DATA_NUM][idx][1] ;
	int iSize = tDat->thSrtPos[POS_LIN][DATA_NUM][idx+1][1] - Start;
	unsigned int uiR[iClusterSize]={0,}, uiT, tmp;
	cTxt_t cT, cY, cZ, c_Rb, *cGamS, cR[iClusterSize][iSize], cU[iClusterSize][iSize], cX[iClusterSize][iSize];
	pTxt_t pR;
	cTxt_t (*cCT)[DATA_NUM];

	#ifdef _DEBUG_INIT_1
	mpz_inits(pR.m, cT.c, cY.c, cZ.c, c_Rb.c, NULL);
	for (int i=0 ; i<iClusterSize ; i++) {
		for (int j=0 ; j<iSize ; j++)
			mpz_inits(cR[i][j].c, cU[i][j].c, cX[i][j].c, NULL);
	}

	#else
	mpz_init2(pR.m, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cY.c, 2*GMP_N_SIZE*2);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	mpz_init2(c_Rb.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iClusterSize ; i++) {
		for (int j=0 ; j<iSize ; j++) {
			mpz_init2(cR[i][j].c, 2*GMP_N_SIZE*2);
			mpz_init2(cU[i][j].c, 2*GMP_N_SIZE*2);
			mpz_init2(cX[i][j].c, 2*GMP_N_SIZE*2);
		}
		mpz_init2(cX[i].c, 2*GMP_N_SIZE*2);
	}
	#endif

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);

	cGamS = tDat->cGamS;
	cCT = tTbl->cCT;
	srand((unsigned int)time(0));
	// ----------------------------------------------------------------------------------------
	// compute permutation phi
	for (int i=0 ; i<iClusterSize ; i++)
		uiR[i] = i;

	#ifdef _DEBUG_Step2
	printf("[DH-%03d] phi permutation : ", idx);
	for (int i=0 ; i<iClusterSize ; i++)
		printf(" %d ", uiR[i]);
	std::cout << std::endl;
	#endif

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	(*c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c1, idx);
	#endif
	if (*c1 < iNumActThds)
		cv->wait(ulSync, [&] {return *c1>=iNumActThds;});
	else {
		// --------------------------------------------------------------------
		printf("\n************   [Step 2] compute CT   ************\n\n");
		tDat->start = time(NULL);

		// --------------------------------------------------------------------
		// computing summation array of gamma(1)-table
		mpz_set(cGamS[0].c, tTbl->cGam2[0][1].c);		// initialization
		for (int i=1 ; i<iClusterSize-1 ; i++) {
			paillier_mul(mPubKey, cGamS+i, cGamS+i-1, tTbl->cGam2[i]+1);
		}
		// --------------------------------------------------------------------
		// send delta[n][j]
		for (int j=0 ; j<2 ; j++) {
			if (j==0)
				SetSendMsg(ucSendPtr, bV, idx, COM_ST2D1, 0);
			else
				SetSendMsg(ucSendPtr, bV, idx, COM_ST2D2, 0);
			for (int i=0 ; i<iClusterSize ; i++) {
				bV[0]=0;
				paillier_ciphertext_to_bytes(bV, ENC_SIZE, tTbl->cGam2[uiR[i]]+j);
				memcpy(ucSendPtr+HED_SIZE+(i*ENC_SIZE), bV, ENC_SIZE);
			}
			ucSendPtr[3] = (unsigned char)((iClusterSize*ENC_SIZE) & 0x000000ff);
			ucSendPtr[4] = (unsigned char)(((iClusterSize*ENC_SIZE) & 0x0000ff00) >> 8);
			#ifdef _DEBUG_Step2
			DebugCom("(buf) delta[i][j] (Hex)\t", ucSendPtr, HED_SIZE+(iClusterSize*ENC_SIZE), idx);
			#endif

			// DH --> CSP
			try {
				mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
				uiCommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}
		}
		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
	}
	ulSync.unlock();

	// ----------------------------------------------------------------------------------------

	// randomize CT
	for (int j=0 ; j<iSize ; j++) {
		for (int i=0 ; i<iClusterSize ; i++) {
			mpz_urandomm(pR.m, state, mPubKey->n);
			mpz_set_ui(pR.m, 0x100);
			paillier_enc(cR[i]+j, mPubKey, &pR, paillier_get_rand_devurandom);
			paillier_mul(mPubKey, cU[i]+j, cCT[i]+(Start+j), cR[i]+j);
		}
	}

	// print gamma[i][j] and u[i][]
	#ifdef _DEBUG_Step2
	for (int j=0 ; j<2 ; j++) {
		for (int i=0 ; i<iClusterSize ; i++) {
			DebugOut("cGamma[i][j] (ciphertext) \t", tTbl->cGam2[i][j].c, idx);
			DebugDec("cGamma[i][j] \t", tTbl->cGam2[i]+j, idx);
		}
	}
	for (int j=0 ; j<iSize ; j++) {
		for (int i=0 ; i<iClusterSize ; i++) {
			DebugOut("cU[i][j] (ciphertext) \t", cU[i][j].c, idx);
			DebugDec("cU[i][j] \t", cU[i]+j, idx);
		}
	}
	#endif

	// send v[n][j]
	for (int j=0 ; j<iSize ; j++) {
		SetSendMsg(ucSendPtr, bV, idx, COM_ST2V, 0);
		for (int i=0 ; i<iClusterSize ; i++) {
			bV[0]=0;
			paillier_ciphertext_to_bytes(bV, ENC_SIZE, cU[uiR[i]]+j);
			memcpy(ucSendPtr+HED_SIZE+(i*ENC_SIZE), bV, ENC_SIZE);
		}
		ucSendPtr[3] = (unsigned char)((iClusterSize*ENC_SIZE) & 0x000000ff);
		ucSendPtr[4] = (unsigned char)(((iClusterSize*ENC_SIZE) & 0x0000ff00) >> 8);
		#ifdef _DEBUG_Step2
		DebugCom("(buf) v[n][j] (Hex):\t", ucSendPtr, HED_SIZE+(iClusterSize*ENC_SIZE), idx);
		#endif

		// DH --> CSP
		try {
			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
			uiCommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
		}
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		// DH <-- CSP
		ulRecv.lock();
		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
		ulRecv.unlock();

		// receive E(w[i])
		short sLen ;
		sLen = Byte2Short(ucRecvPtr+HED_LEN);
		uiCommAmt[1] += sLen;
		#ifdef _DEBUG_Step2
		DebugCom("(buf) E(w[i]) (Hex)\t", ucRecvPtr, sLen+HED_SIZE, idx);
		#endif

		#ifdef _DEBUG_Assert
		sLen = Byte2Short(ucRecvPtr+HED_LEN);
		assert((*(ucRecvPtr+2)==COM_ST2V)&&(sLen==ENC_SIZE*iClusterSize));
		#endif

		for (int i=0 ; i<iClusterSize ; i++) {
			paillier_ciphertext_from_bytes(cX[uiR[i]]+j, ucRecvPtr+HED_SIZE+(ENC_SIZE*i), ENC_SIZE);
			#ifdef _DEBUG_Step2
			//DebugOut("cX[i][j]' (ciphertext)\t", cX[uiR[i]][j].c, idx);
			DebugDec("cX[i][j]' \t\t", cX[uiR[i]]+j, idx);
			#endif
		}

		// memory initialization
		memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE*iClusterSize);
		ucRecvPtr[1] = 0xff;
	}

	for (int j=0 ; j<iSize ; j++) {

		// E(-r_(beta,j))
		mpz_set_ui(c_Rb.c, 1);		// initialize T
		for (int i=0 ; i<iClusterSize ; i++) {
			paillier_exp(mPubKey, cR[i]+j, cR[i]+j, &(tVal->pN1));
			SecMul(&cT, tTbl->cGam2[i]+1, cR[i]+j, idx, tRecv, uiCommAmt);
			paillier_mul(mPubKey, &c_Rb, &c_Rb, &cT);
		}

		// remove random value
		for (int i=0 ; i<iClusterSize ; i++) {
			SecMul(&cY, tTbl->cGam2[i]+0, &c_Rb, idx, tRecv, uiCommAmt);
			paillier_mul(mPubKey, &cY, cR[i]+j, &cY);
			paillier_mul(mPubKey, cCT[i]+(Start+j), cX[i]+j, &cY);
		}

		#ifdef _DEBUG_Step2
		for (int i=0 ; i<iClusterSize ; i++) {
			DebugDec("cCT[i,j]' \t\t", cCT[i]+(Start+j), idx);
		}
		#endif

		// reduce CT
		for (int i=0 ; i<iClusterSize-1 ; i++) {
			paillier_exp(mPubKey, &cT, cCT[i]+(Start+j), &(tVal->pN1));
			paillier_mul(mPubKey, &cT, cCT[i+1]+(Start+j), &cT);
			SecMul(&cZ, cGamS+i, &cT, idx, tRecv, uiCommAmt);
			paillier_mul(mPubKey, cCT[i]+(Start+j), &cZ, cCT[i]+(Start+j));
		}
	}

	// print computation results
	#ifdef _DEBUG_Step2
	for (int i=0 ; i<iClusterSize-1 ; i++) {
		printf("-- cluster %d -- \n", i+1);
		for (int j=0 ; j<iSize ; j++)
			DebugDec("cCT[i,j]' \t\t", cCT[i]+(Start+j), idx);
	}
	#endif

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	tDat->uiCommAmt[DATA_NUM - tTbl->n][1][0] += uiCommAmt[0];
	tDat->uiCommAmt[DATA_NUM - tTbl->n][1][1] += uiCommAmt[1];
    (*c2)++;
    if (*c2 < iNumActThds)
    	cv->wait(ulSync, [&] {return *c2>=iNumActThds;});
    else {
    	tDat->end = time(NULL);
    	tDat->dTimeRnd[DATA_NUM - tTbl->n][1] = (double) (tDat->end - tDat->start);
    	printf("[DH-%03d] Time \t     : %.2f (seconds) \n", idx, tDat->dTimeRnd[DATA_NUM - tTbl->n][1]);
    	printf("[DH-%03d] Comm Amount : (DH) %lu \t(CSP) %lu \t(bytes) \n", idx, tDat->uiCommAmt[DATA_NUM - tTbl->n][1][0], tDat->uiCommAmt[DATA_NUM - tTbl->n][1][1]);

		// --------------------------------------------------------------------
		//DebugCT("Step 2", tTbl, -1, idx);
		// --------------------------------------------------------------------
		*c1=0;		cv->notify_all();
    }
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	return;
}

// step 3
void PaillierCrypto::step_3(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync)
{
	int iClusterSize = tTbl->n;
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned long uiCommAmt[2]={0, 0};
	int iNumActThds = tDat->iNumActThds;

	cTxt_t cLocalDmin_ij[iClusterSize-1][DATA_SQUARE_LENGTH];
	cTxt_t cT, cW, cX, cY, cZ;
	cTxt_t *cDT4_ij, *cDT4_ji, *cDT3_ij, *cDmin_i ;
	bool	bStart, bExit;
	int Start, End, Start_i, Start_j, End_i, End_j;
	int RemRow, RemTask;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, cW.c, cX.c, cY.c, cZ.c, NULL);
	for (int i=0 ; i<iClusterSize-1 ; i++)
		for (int j=0 ; j<DATA_SQUARE_LENGTH ; j++)
			mpz_init(cLocalDmin_ij[i][j].c);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(cX.c, 2*GMP_N_SIZE*2);
	mpz_init2(cY.c, 2*GMP_N_SIZE*2);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	for (int i=1 ; i<iClusterSize-1 ; i++)
		for (int j=0 ; j<DATA_SQUARE_LENGTH ; j++)
			mpz_init2(cLocalDmin_ij[i][j].c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);
	// ----------------------------------------------------------------------------------------

	ulSync.lock();
	(*c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c1, idx);
	#endif
	if (*c1 < iNumActThds)
		cv->wait(ulSync, [&] {return *c1>=iNumActThds;});
	else {
		// --------------------------------------------------------------------
		printf("\n************   [Step 3] compute CDT   ************\n\n");
		tDat->start = time(NULL);

		// --------------------------------------------------------------------

		#ifdef _DEBUG_Step3
		DebugDT("before copying upper triangular part", tTbl, 0, 4, idx);
		#endif

		// copy upper triangular part
		for (int i=1 ; i<iClusterSize ; i++) {
			for (int j=0 ; j<i ; j++) {
				cDT4_ij = tTbl->cDT4[i][j];
				cDT4_ji = tTbl->cDT4[j][i];
				for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
					mpz_set(cDT4_ij[k].c, cDT4_ji[k].c);
				}
			}
		}

		#ifdef _DEBUG_Step3
		for (int i=0 ; i<iClusterSize ; i++) {
			for (int j=0 ; j<iClusterSize ; j++) {
				cDT4_ij = tTbl->cDT4[i][j];
				for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
					printf("(%02d, %02d, %02d) ", i,j,k);
					//DebugOut("cDT[i,j,k] (ciphertext)\t", cDT4_ij[k].c, idx);
					DebugDec("cDT[i,j,k] \t\t", cDT4_ij+k, idx);
				}
			}
		}
		#endif

		#ifdef _DEBUG_Step3
		DebugDT("after copying upper triangular part", tTbl, 0, 4, idx);
		#endif

		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
	}
	ulSync.unlock();

	// ----------------------------------------------------------------------------------------

	Start = tDat->thSrtPos[POS_LIN][iClusterSize][idx][1] ;
	End   = tDat->thSrtPos[POS_LIN][iClusterSize][idx+1][1] ;

	// initialize cDmin_i[i][k]
	if (Start >= 0) {
		for (int i=Start ; i<End ; i++) {
			cDmin_i = tDat->cDmin_i[i];
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
				mpz_set_ui(cDmin_i[k].c, 1);
		}
	}

	Start = tDat->thSrtPos[POS_ST3][iClusterSize][idx][0] ;
	End   = tDat->thSrtPos[POS_ST3][iClusterSize][idx+1][0] ;

	// compute distances between new merged cluster and other clusters
	if (Start >= 0) {
		for (int i=Start ; i<End ; i++) {
			SmaxSmin_PSV(cLocalDmin_ij, i, tTbl, tDat, tVal, idx, tRecv, uiCommAmt);
			cDmin_i = tDat->cDmin_i[i];
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				for (int j=0 ; j<iClusterSize-1 ; j++)
					paillier_mul(mPubKey, cDmin_i+k, cDmin_i+k, cLocalDmin_ij[j]+k);
				#ifdef _DEBUG_Step3
				DebugDec("cDmin[i,k]: \t\t", cDmin_i+k, idx);
				#endif
			}
		}
	}

	// threads compute the remaining tasks (rows)
	RemRow = tDat->thSrtPos[POS_ST3][iClusterSize][idx][1];
	if (RemRow >= 0) {
		RemTask = RemRow % iNumActThds;
		SmaxSmin_PSV(tDat->cDmin_ij[RemTask], tDat->cDmin_i, RemRow, tTbl, tDat, tVal, idx, tRecv, tSync, uiCommAmt);

		#ifdef _DEBUG_Step3
		cDmin_i = tDat->cDmin_i[RemRow];
		for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++)
			DebugDec("cDmin[i,k]: \t\t", cDmin_i+k, idx);
		#endif
	}

	bStart=true, bExit = false;
	Start_i = tDat->thSrtPos[POS_SQU][iClusterSize][idx][0] ;	Start_j = tDat->thSrtPos[POS_SQU][iClusterSize][idx][1] ;
	End_i 	 = tDat->thSrtPos[POS_SQU][iClusterSize][idx+1][0];	End_j 	 = tDat->thSrtPos[POS_SQU][iClusterSize][idx+1][1] ;

	// set distances between new merged cluster and other clusters
	for (int i=0 ; i<iClusterSize ; i++) {
		if (bExit)			break;
		for (int j=0 ; j<iClusterSize ; j++) {
			if (bStart) {
				i = Start_i ;		j = Start_j ;
				bStart = false;
			}
			if ((i>=End_i) && (j>=End_j)) {
				bExit = true;
				break;
			}

			if (i == j)	continue;
			cDT4_ij = tTbl->cDT4[i][j];
			cDmin_i = tDat->cDmin_i[i];
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				paillier_exp(mPubKey, &cT, cDT4_ij+k, &(tVal->pN1));
				paillier_mul(mPubKey, &cW, cDmin_i+k, &cT);
				SecMul(&cX, tTbl->cGam2[j]+0, &cW, idx, tRecv, uiCommAmt);
				paillier_mul(mPubKey, cDT4_ij+k, cDT4_ij+k, &cX);

				#ifdef _DEBUG_Step3
				printf("(%02d, %02d, %02d) ", i,j,k);
				//DebugOut("cDT[i,j,k] (ciphertext)\t", cDT4_ij[k].c, idx);
				DebugDec("cDT[i,j,k] \t\t", cDT4_ij+k, idx);
				#endif
			}
		}
	}

	bStart=true, bExit = false;
	Start_i = tDat->thSrtPos[POS_TRI][iClusterSize][idx][0] ;	Start_j = tDat->thSrtPos[POS_TRI][iClusterSize][idx][1] ;
	End_i 	 = tDat->thSrtPos[POS_TRI][iClusterSize][idx+1][0];	End_j 	 = tDat->thSrtPos[POS_TRI][iClusterSize][idx+1][1] ;

	// copy lower triangular part
	for (int i=0 ; i<iClusterSize-1 ; i++) {
		if (bExit)			break;
		for (int j=i+1 ; j<iClusterSize ; j++) {
			if (bStart) {
				i = Start_i ;		j = Start_j ;
				bStart = false;
			}
			if ((i>=End_i) && (j>=End_j)) {
				bExit = true;
				break;
			}
			cDT4_ij = tTbl->cDT4[i][j];
			cDT4_ji = tTbl->cDT4[j][i];
			cDT3_ij = tTbl->cDT3[i][j];
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				paillier_exp(mPubKey, &cT, cDT4_ji+k, &(tVal->pN1));
				paillier_mul(mPubKey, &cY, cDT4_ij+k, &cT);
				SecMul(&cZ, tTbl->cGam1+j, &cY, idx, tRecv, uiCommAmt);
				paillier_mul(mPubKey, cDT3_ij+k, cDT4_ji+k, &cZ);

				#ifdef _DEBUG_Step3
				printf("(%02d, %02d, %02d) ", i,j,k);
				//DebugOut("cDT[i,j,k] (ciphertext)\t", cDT3_ij[k].c, idx);
				DebugDec("cDT[i,j,k] \t\t", cDT3_ij+k, idx);
				#endif
			}
		}
	}

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	tDat->uiCommAmt[DATA_NUM - tTbl->n][2][0] += uiCommAmt[0];
	tDat->uiCommAmt[DATA_NUM - tTbl->n][2][1] += uiCommAmt[1];
    (*c2)++;
    if (*c2 < iNumActThds)
    	cv->wait(ulSync, [&] {return *c2>=iNumActThds;});
    else {
    	tDat->end = time(NULL);
    	tDat->dTimeRnd[DATA_NUM - tTbl->n][2] = (double) (tDat->end - tDat->start);
    	printf("[DH-%03d] Time \t     : %.2f (seconds) \n", idx, tDat->dTimeRnd[DATA_NUM - tTbl->n][2]);
    	printf("[DH-%03d] Comm Amount : (DH) %lu \t(CSP) %lu \t(bytes) \n", idx, tDat->uiCommAmt[DATA_NUM - tTbl->n][2][0], tDat->uiCommAmt[DATA_NUM - tTbl->n][2][1]);
		// --------------------------------------------------------------------
		#ifdef _DEBUG_Detail
    	for (int i=0 ; i<iClusterSize ; i++) {
    		for (int j=i+1 ; j<iClusterSize ; j++) {
			    for (int k=DATA_SQUARE_LENGTH-1 ; k>=0 ; k--) {
					DebugDecBit(tTbl->cDT3[i][j]+k);
					if (k%8==0)	printf("\t");
					if (k%64==0)	printf("\n");
				}
			}
		}
		#endif

    	//DebugDT("Step 3", tTbl, 0, 3, idx);
		// --------------------------------------------------------------------
    	*c1=0;		cv->notify_all();
    }
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	return;
}

// step 4
void PaillierCrypto::step_4(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned long uiCommAmt[2]={0, 0};
	int iNumActThds = tDat->iNumActThds;

	int iReducedSize = tTbl->n - 1 ;
	cTxt_t cT, cU, cV, cW, cX, cY, cZ;
	cTxt_t (*cBs)[DATA_NUM-1], *cGs, *cBeta_i, *cDT3_11, *cDT3_01, *cDT3_00, *cDT4_00;
	bool	bStart=true, bExit = false;
	int Start_i, Start_j, End_i, End_j;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, cU.c, cV.c, cW.c, cX.c, cY.c, cZ.c, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cU.c, 2*GMP_N_SIZE*2);
	mpz_init2(cV.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(cX.c, 2*GMP_N_SIZE*2);
	mpz_init2(cY.c, 2*GMP_N_SIZE*2);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);

	cBs = tDat->cBs;
	cGs = tDat->cGs;
	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	(*c1)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c1, idx);
	#endif
	if (*c1 < iNumActThds)
		cv->wait(ulSync, [&] {return *c1>=iNumActThds;});
	else {
		// --------------------------------------------------------------------
		printf("\n************   [Step 4] compute RDT   ************\n\n");
		tDat->start = time(NULL);

		// --------------------------------------------------------------------
		cTxt_t cB[iReducedSize][iReducedSize];

		for (int i=0 ; i<iReducedSize ; i++)
			for (int j=0 ; j<iReducedSize ; j++)
				mpz_init(cB[i][j].c);

		// computing b-table and g-array
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			cBeta_i = tTbl->cBeta[i];
			paillier_exp(mPubKey, &cT, cBeta_i+i, &(tVal->pN1));
			for (int j=i+1 ; j<iReducedSize ; j++) {
				paillier_mul(mPubKey, cB[i]+j, cBeta_i+j, &cT);		// b-table
			}
		}

		#ifdef _DEBUG_Step4
		for (int i=0 ; i<iReducedSize ; i++) {
			for (int j=i+1 ; j<iReducedSize ; j++)
				DebugDec("cB-table[i,j]: \t\t", cB[i]+j, idx);
		}
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			DebugDec("cG-array[i]: \t\t", tTbl->cBeta[i]+i, idx);
		}
		#endif

		// initializing b^s-table and g^s-array
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			mpz_set_ui(cGs[i].c, 1);									//	g^s-array
			for (int j=i+1 ; j<iReducedSize ; j++) {
				mpz_set_ui(cBs[i][j].c, 1);							//	b^s-table
			}
		}

		// computing b^s-table and g^s-array
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			mpz_set(cBs[i][i+1].c, cB[i][i+1].c);
			for (int j=i+2 ; j<iReducedSize ; j++) {
				paillier_mul(mPubKey, cBs[i]+j, cBs[i]+j-1, cB[i]+j);	//	b^s-table
			}
		}

		mpz_set(cGs[0].c, tTbl->cBeta[0][0].c);
		for (int i=1 ; i<iReducedSize-1 ; i++) {
			paillier_mul(mPubKey, cGs+i, cGs+i-1, tTbl->cBeta[i]+i);	//	g^s-array
		}

		#ifdef _DEBUG_Step4
		for (int i=0 ; i<iReducedSize ; i++) {
			for (int j=i+1 ; j<iReducedSize ; j++) {
				printf("(%d, %d) ", i, j);
				DebugDec("cB^s-table[i,j]: \t\t", cBs[i]+j, idx);
			}
		}
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			printf("(%d) ", i);
			DebugDec("cG^s-array[i]: \t\t", cGs+i, idx);
		}
		#endif

		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
	}
	ulSync.unlock();
	// ----------------------------------------------------------------------------------------

	Start_i = tDat->thSrtPos[POS_TRI][iReducedSize][idx][0];	Start_j = tDat->thSrtPos[POS_TRI][iReducedSize][idx][1] ;
	End_i 	 = tDat->thSrtPos[POS_TRI][iReducedSize][idx+1][0];	End_j 	 = tDat->thSrtPos[POS_TRI][iReducedSize][idx+1][1] ;

	// reducing DT (RDT)
	for (int i=0 ; i<iReducedSize-1 ; i++) {
		if (bExit)			break;
		paillier_mul(mPubKey, &cU, cGs+i, &(tVal->cN1));
		for (int j=i+1 ; j<iReducedSize ; j++) {
			if (bStart) {
				i = Start_i ;		j = Start_j ;
				bStart = false;
				paillier_mul(mPubKey, &cU, cGs+i, &(tVal->cN1));
			}
			if ((i>=End_i) && (j>=End_j)) {
				bExit = true;
				break;
			}
			paillier_mul(mPubKey, &cV, cBs[i]+j, &(tVal->cN1));
			SecMul(&cW, &cU, &cV, idx, tRecv, uiCommAmt);

			cDT3_11 = tTbl->cDT3[i+1][j+1];
			cDT3_01 = tTbl->cDT3[i]  [j+1];
			cDT3_00 = tTbl->cDT3[i]  [j];
			cDT4_00 = tTbl->cDT4[i]  [j];
			for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
				SecMul(&cX, cGs+i, cDT3_11+k, idx, tRecv, uiCommAmt);
				SecMul(&cY, cBs[i]+j, cDT3_01+k, idx, tRecv, uiCommAmt);
				SecMul(&cZ, &cW, cDT3_00+k, idx, tRecv, uiCommAmt);
				paillier_mul(mPubKey, &cT, &cX, &cY);
				paillier_mul(mPubKey, cDT4_00+k, &cT, &cZ);

				// print computation results
				#ifdef _DEBUG_Step4
				printf("(%02d, %02d, %02d) ", i,j,k);
				//DebugOut("cDT[i,j,k] (ciphertext)\t", cDT4_00[k].c, idx);
				DebugDec("cDT[i,j,k] \t\t", cDT4_00+k, idx);
				#endif
			}
		}
	}

	// ----------------------------------------------------------------------------------------
	ulSync.lock();
	tDat->uiCommAmt[DATA_NUM - tTbl->n][3][0] += uiCommAmt[0];
	tDat->uiCommAmt[DATA_NUM - tTbl->n][3][1] += uiCommAmt[1];
    (*c2)++;
	#ifdef _DEBUG_THREAD
	printf("(Count: %02d) <<<  %03d - th Thread : continuing  >>>\n", *c2, idx);
	#endif
    if (*c2 < iNumActThds)
    	cv->wait(ulSync, [&] {return *c2>=iNumActThds;});
    else {
		// --------------------------------------------------------------------
    	tDat->end = time(NULL);
    	int iRnd = DATA_NUM - tTbl->n;
    	tDat->dTimeRnd[iRnd][3] = (double) (tDat->end - tDat->start);
    	printf("[DH-%03d] Time \t     : %.2f (seconds) \n", idx, tDat->dTimeRnd[iRnd][3]);
    	printf("[DH-%03d] Comm Amount : (DH) %lu \t(CSP) %lu \t(bytes) \n", idx, tDat->uiCommAmt[iRnd][3][0], tDat->uiCommAmt[iRnd][3][1]);

		#ifdef _DEBUG_Detail
		for (int i=0 ; i<iReducedSize-1 ; i++) {
			for (int j=i+1 ; j<iReducedSize ; j++) {
			    for (int k=DATA_SQUARE_LENGTH-1 ; k>=0 ; k--) {
					DebugDecBit(tTbl->cDT4[i][j]+k);
					if (k%8==0)	printf("\t");
					if (k%64==0)	printf("\n");
				}
			}
		}
		#endif

		//DebugDT("step 4", tTbl, -1, 4, idx);
		// --------------------------------------------------------------------
		// print the running time and communication amount of a round (steps 1-4)
		// --------------------------------------------------------------------
		tTbl->n = tTbl->n - 1;
		printf("-----------------------------------------------------------------------------------------------------------------\n");
    	printf("[DH-%03d] %02d rounds (cluster size: %02d): %d (seconds)\t ( %d, \t %d, \t %d, \t %d ) \n",
    			idx, iRnd+1, tTbl->n, (int)(tDat->dTimeRnd[iRnd][0]+tDat->dTimeRnd[iRnd][1]+tDat->dTimeRnd[iRnd][2]+tDat->dTimeRnd[iRnd][3]),
				(int)tDat->dTimeRnd[iRnd][0], (int)tDat->dTimeRnd[iRnd][1], (int)tDat->dTimeRnd[iRnd][2], (int)tDat->dTimeRnd[iRnd][3]);
    	printf("[DH-%03d] (DH/CSP) %lu / %lu (bytes)\t( %lu / %lu ,\t%lu / %lu ,\t%lu / %lu ,\t%lu / %lu ) \n",
    			idx, tDat->uiCommAmt[iRnd][0][0]+tDat->uiCommAmt[iRnd][1][0]+tDat->uiCommAmt[iRnd][2][0]+tDat->uiCommAmt[iRnd][3][0],
				tDat->uiCommAmt[iRnd][0][1]+tDat->uiCommAmt[iRnd][1][1]+tDat->uiCommAmt[iRnd][2][1]+tDat->uiCommAmt[iRnd][3][1],
				tDat->uiCommAmt[iRnd][0][0], tDat->uiCommAmt[iRnd][0][1], tDat->uiCommAmt[iRnd][1][0], tDat->uiCommAmt[iRnd][1][1],
				tDat->uiCommAmt[iRnd][2][0], tDat->uiCommAmt[iRnd][2][1], tDat->uiCommAmt[iRnd][3][0], tDat->uiCommAmt[iRnd][3][1]);
		printf("-----------------------------------------------------------------------------------------------------------------\n\n");

		// --------------------------------------------------------------------

		if (iNumActThds > tTbl->n) {
			tDat->iNumActThds--;
		}

		*c1=0;		cv->notify_all();
    }
	ulSync.unlock();

	// ----------------------------------------------------------------------------------------

	return;
}

// SMin^ASB protocol
void PaillierCrypto::SmaxSmin_ASB(cTxt_t cMDT[][DATA_NUM],
		shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, sync_t* tSync, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bW[ENC_SIZE]={0,};
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	int iNumActThds = tDat->iNumActThds;

	int iClusterSize = tTbl->n;
	int iSizeTotal = iClusterSize*(iClusterSize-1)/2;
	int iSize = tDat->thSrtPos[POS_LIN][iClusterSize][idx+1][0] - tDat->thSrtPos[POS_LIN][iClusterSize][idx][0];
	unsigned int uiR[iSize], uiT, tmp;
	cTxt_t cS, cW, *cC, *cU;
	pTxt_t pR;
	int Start, End;
	int perPos;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cS.c, cW.c, pR.m, NULL);
	#else
	mpz_init2(cS.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(pR.c, 2*GMP_N_SIZE*2);
	#endif

	srand((unsigned int)time(0));

	// ----------------------------------------------------------------------------------------
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);

	cC = tDat->cC[0] ;
	cU = tDat->cU[0] ;
	// ----------------------------------------------------------------------------------------

	Start = tDat->thSrtPos[POS_LIN][iClusterSize][idx][0] ;
	End   = tDat->thSrtPos[POS_LIN][iClusterSize][idx+1][0] ;

	// all candidate mode: C[i,j] = E(1)
	for (int i=Start ; i<End ; i++)
		mpz_set(cC[i].c, tVal->c1.c);

	// SMin protocol
	SmaxSmin(cC+Start, iSize, tDat->DTposST1[iClusterSize], Start, tDat->cS[0], &(tDat->cSum[0]), tTbl, tDat, tVal, idx, tRecv, tSync, CommAmt);

	// ------------------------------------------------------------------------------------
	ulSync.lock();
    (*c1)++;
    if (*c1 < iNumActThds)
    	cv->wait(ulSync, [&] {return *c1>=iNumActThds;});
    else {
		// --------------------------------------------------------------------
		// print the result of SMIN
		#ifdef _DEBUG_Detail
		for (int i=0 ; i<iSizeTotal ; i++) {
			printf("(%02d) ", i);
			DebugDec("C[i]: \t\t", cC+i, idx);
		}
		#endif
		// --------------------------------------------------------------------

		// E(u[i]) = sum_j=1~i-1 E(C[j])
		mpz_set_ui(cU[0].c, 1);
		for (int i=1 ; i<iSizeTotal ; i++)
			paillier_mul(mPubKey, cU+i, cU+(i-1), cC+(i-1));

		#ifdef _DEBUG_SMIN_ASB
		for (int i=0 ; i<iSizeTotal ; i++)
			DebugDec("cU: \t\t", cU+i, idx);
		#endif

		// --------------------------------------------------------------------
		*c2=0;		cv->notify_all();
    }
	ulSync.unlock();
	// ------------------------------------------------------------------------------------

	// phi permutation idx / depermutation idx
	for (int i=0 ; i<iSize ; i++)
		uiR[i] = i;
	for (int i=iSize-1 ; i>0 ; i--) {
		tmp = rand() % (i+1);
		uiT = uiR[i];
		uiR[i] = uiR[tmp];
		uiR[tmp] = uiT;
	}

	#ifdef _DEBUG_SMIN_ASB
	printf("[DH-%03d] phi permutation : ", idx);
	for (int i=0 ; i<iSize ; i++)
		printf(" %d ", uiR[i]);
	std::cout << std::endl;
	#endif

	for (int i=0 ; i<iSize ; i++) {
		perPos = uiR[i];
		#ifdef _DEBUG_SMIN_ASB
			DebugDec("(permuted) cC: \t\t", cC+Start+perPos, idx);
			DebugDec("(permuted) cU: \t\t", cU+Start+perPos, idx);
		#endif

		paillier_mul(mPubKey, &cS, cU+Start+perPos, cU+Start+perPos);
		paillier_mul(mPubKey, &cS, cC+Start+perPos, &cS);
		paillier_mul(mPubKey, &cS, &(tVal->cN1), &cS);

		#ifdef _DEBUG_SMIN_ASB
			DebugDec("(permuted) cS: \t\t", &cS, idx);
		#endif

		// E(w[i]) = E(s[i])^r_i
		mpz_urandomm(pR.m, state, mPubKey->n);
		paillier_exp(mPubKey, &cW, &cS, &pR);

		#ifdef _DEBUG_SMIN_ASB
			DebugDec("cW: \t\t", &cW, idx);
		#endif

		//------------------------------------------------------------------------------
		// send w[i]'
		bW[0]=0;
		paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cW);
		SetSendMsg(ucSendPtr, bW, idx, COM_SING, ENC_SIZE);
		#ifdef _DEBUG_SMIN_ASB
	    DebugCom("(buf) E(x_i) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
		#endif

	    // DH --> CSP
		try {
			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
			CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
		}
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		// DH <-- CSP
		ulRecv.lock();
		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
		ulRecv.unlock();

		// receive E(y[i])
		short sLen ;
		sLen = Byte2Short(ucRecvPtr+HED_LEN);
		CommAmt[1] += sLen;
		#ifdef _DEBUG_SMIN_ASB
		//sLen = Byte2Short(ucRecvPtr+HED_LEN);
	    DebugCom("(buf) E(y[i]) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
		#endif

		#ifdef _DEBUG_Assert
	    sLen = Byte2Short(ucRecvPtr+HED_LEN);
	    assert((*(ucRecvPtr+2)==COM_SING)&&(sLen==ENC_SIZE));
		#endif

	    int row = tDat->DTposST1[iClusterSize][Start+perPos][0];
	    int col = tDat->DTposST1[iClusterSize][Start+perPos][1];
	    paillier_ciphertext_from_bytes(&cMDT[row][col], ucRecvPtr+HED_SIZE, ENC_SIZE);
	    #ifdef _DEBUG_SMIN_ASB
		DebugDec("cMDT[row][col]: \t\t", &cMDT[row][col], idx);
		#endif
		//------------------------------------------------------------------------------

		memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
		ucRecvPtr[1] = 0xff;
	}

	return;
}

// SMin protocol for SMin^ASB
void PaillierCrypto::SmaxSmin(cTxt_t *cC, const int iSize, const int DTpos[][2], int Start, cTxt_t *cLocalS_th, cTxt_t *cGlobalS,
		 shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *c1, *c2;
	std::condition_variable *cv;
	unsigned int *pSyncCnt2, *pSyncCnt1;
	int iNumActThds = tDat->iNumActThds;

	cTxt_t *cD, *cLocalS, *cQ, cW, cP[iSize] ;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cW.c, NULL);
	for (int i=0 ; i<iSize ; i++)
		mpz_inits(cP[i].c, NULL);
	#else
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iSize ; i++)
		mpz_init2(cP[i].c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	cLocalS = cLocalS_th+idx;
	cQ = &(tDat->cQ[0]) ;
	c1 = &(tSync->c1);
	c2 = &(tSync->c2);
	cv = &(tSync->cv);
	// ----------------------------------------------------------------------------------------

	for (int j=DATA_SQUARE_LENGTH-1 ; j>=0 ; j--) {

		#ifdef _DEBUG_Detail
		printf("\n---------  Round %02d  ---------\n", j);
		#endif

		// --------------------------------------------------------------------
		// Step 1: compute P[i] = SM(C[i], d[i,j])
		mpz_set_ui(cLocalS->c, 1);		// initialize local S

		for (int i=0 ; i<iSize ; i++) {
			cD = tTbl->cDT4[DTpos[Start+i][0]][DTpos[Start+i][1]]+j ;

			#ifdef _DEBUG_SMIN
			DebugDec("d[i,j] \t\t\t", cD, idx);
			DebugDec("C[i] \t\t\t", cC+i, idx);
			#endif

			// P[i] = SM(C[i], d[i,j])
			SecMul(cP+i, cC+i, cD, idx, tRecv, CommAmt);

			// step 2: local S = Sum_(i=1~n) P[i]
			paillier_mul(mPubKey, cLocalS, cLocalS, cP+i);

			#ifdef _DEBUG_SMIN
			DebugDec("P[i] = SM(C[i], d[i,j]) \t", cP+i, idx);
			DebugDec("local S = Sum_(i=1~n) P[i] \t", cLocalS, idx);
			#endif
		}
		// --------------------------------------------------------------------
		ulSync.lock();
		#if (DATA_SQUARE_LENGTH-1) % 2 == 0
		if (j%2 == 0) 		{	pSyncCnt1 = c1;	pSyncCnt2 = c2;	}
		else 					{	pSyncCnt1 = c2;	pSyncCnt2 = c1;	}
		#else
		if (j%2 == 0) 		{	pSyncCnt1 = c2;	pSyncCnt2 = c1;	}
		else 					{	pSyncCnt1 = c1;	pSyncCnt2 = c2;	}
		#endif
		(*pSyncCnt2)++;
		if (*pSyncCnt2 < iNumActThds)
			cv->wait(ulSync, [&] {return *pSyncCnt2>=iNumActThds;});
		else {
			// -----------------------------------------------------------------
			// Step 2: compute global S
			mpz_set_ui(cGlobalS->c, 1);		// initialize global S
			for (int i=0 ; i<iNumActThds ; i++) {
				#ifdef _DEBUG_SMIN
				DebugDec("local S = Sum_(i=1~n) P[i] \t", cLocalS_th+i, idx);
				#endif
				paillier_mul(mPubKey, cGlobalS, cGlobalS, cLocalS_th+i);
			}

			#ifdef _DEBUG_SMIN
			DebugDec("compute global S \t", cGlobalS, idx);
			#endif
			// -----------------------------------------------------------------
			// Step 3: compare q = SEQ(S, 0)
			SEQ(cQ, cGlobalS, &(tVal->c0), tVal, idx, tRecv, CommAmt);

			#ifdef _DEBUG_SMIN
			DebugDec("Q = \t\t\t\t", cQ, idx);
			#endif
			// -----------------------------------------------------------------
			*pSyncCnt1=0;		cv->notify_all();
		}
		ulSync.unlock();
		// --------------------------------------------------------------------
		// Step 4: compute C[i] = SM(q, C[i]) * P[i]

		for (int i=0 ; i<iSize ; i++) {
			// w = SM(q, C[i])
			SecMul(&cW, cQ, cC+i, idx, tRecv, CommAmt);

			// C[i] = w * P[i]
			paillier_mul(mPubKey, cC+i, cP+i, &cW);

			#ifdef _DEBUG_SMIN
			DebugDec("w = SM(q, C[i]) \t\t", &cW, idx);
			DebugDec("C[i] = w * P[i] \t\t", cC+i, idx);
			#endif
		}
	}

	// ------------------------------------------------------------------------------------
	if (DATA_SQUARE_LENGTH%2 == 0) {
		ulSync.lock();
		(*c2)++;
		if (*c2 < iNumActThds)
			cv->wait(ulSync, [&] {return *c2>=iNumActThds;});
		else {
			*c1=0;
			cv->notify_all();
		}
		ulSync.unlock();
	}

	return;
}

// SMin^PSV protocol for single-thread
void PaillierCrypto::SmaxSmin_PSV(cTxt_t cV[][DATA_SQUARE_LENGTH],
		const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bW[ENC_SIZE]={0,};

	int iClusterSize = tTbl->n;
	int iSize = iClusterSize-1;
	unsigned int uiR[iSize], uiT, tmp;
	cTxt_t cS, cW, cZ, cC[iSize], cU[iSize];
	pTxt_t pR;
	int perPos;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cS.c, cW.c, cZ.c, pR.m, NULL);
	for (int i=0 ; i<iSize ; i++)
		mpz_inits(cC[i].c, cU[i].c, NULL);
	#else
	mpz_init2(cS.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(pR.m, 2*GMP_N_SIZE*2);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iSize ; i++) {
		mpz_init2(cC[i].c, 2*GMP_N_SIZE*2);
		mpz_init2(cU[i].c, 2*GMP_N_SIZE*2);
	}
	#endif

	srand((unsigned int)time(0));

	// ----------------------------------------------------------------------------------------
	#ifdef _DEBUG_Detail
	printf("\n************   (single-thread) SMAX/SMIN^(P,S,V) - %d-th row   ************\n\n", row);
	#endif

	// part candidate mode: C[i] = Gamma[i]
	for (int j=0 ; j<iSize ; j++) {
		tmp = tDat->DTposST3[iClusterSize][row][j];
		mpz_set(cC[j].c, tTbl->cGam1[tmp].c);
	}

	#ifdef _DEBUG_SMIN_PSV
	for (int j=0 ; j<iSize ; j++) {
		printf("cDT position (%d, %d) ", row, tDat->DTposST3[iClusterSize][row][j]);
		DebugDec("cC[j]: \t\t", cC+j, idx);
	}
	#endif
	// ----------------------------------------------------------------------------------------

	// SMIN protocol
	SmaxSmin(cC, iSize, tDat->DTposST3[iClusterSize][row], row, tTbl, tDat, tVal, idx, tRecv, CommAmt);

	// --------------------------------------------------------------------
	// E(u[i]) = sum_j=1~i-1 E(C[j])
	mpz_set_ui(cU[0].c, 1);
	for (int i=1 ; i<iSize ; i++)
		paillier_mul(mPubKey, cU+i, cU+(i-1), cC+(i-1));

	#ifdef _DEBUG_SMIN_PSV
	for (int i=0 ; i<iSize ; i++)
		DebugDec("cU: \t\t", cU+i, idx);
	#endif

	// phi permutation idx / depermutation idx
	for (int i=0 ; i<iSize ; i++)
		uiR[i] = i;
	for (int i=iSize-1 ; i>0 ; i--) {
		tmp = rand() % (i+1);
		uiT = uiR[i];
		uiR[i] = uiR[tmp];
		uiR[tmp] = uiT;
	}

	#ifdef _DEBUG_SMIN_PSV
	printf("[DH-%03d] phi permutation : ", idx);
	for (int i=0 ; i<iSize ; i++)
		printf(" %d ", uiR[i]);
	std::cout << std::endl;
	#endif

	// --------------------------------------------------------------------
	for (int i=0 ; i<iSize ; i++) {
		perPos = uiR[i];
		#ifdef _DEBUG_SMIN_PSV
		DebugDec("(permuted) cC: \t\t", cC+perPos, idx);
		DebugDec("(permuted) cU: \t\t", cU+perPos, idx);
		#endif

		paillier_mul(mPubKey, &cS, cU+perPos, cU+perPos);
		paillier_mul(mPubKey, &cS, cC+perPos, &cS);
		paillier_mul(mPubKey, &cS, &(tVal->cN1), &cS);

		#ifdef _DEBUG_SMIN_PSV
			DebugDec("(permuted) cS: \t\t", &cS, idx);
		#endif

		// E(w[i]) = E(s[i])^r_i
		mpz_urandomm(pR.m, state, mPubKey->n);
		paillier_exp(mPubKey, &cW, &cS, &pR);

		#ifdef _DEBUG_SMIN_PSV
			DebugDec("cW: \t\t", &cW, idx);
		#endif

		//------------------------------------------------------------------------------
		// send w[i]'
		bW[0]=0;
		paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cW);
		SetSendMsg(ucSendPtr, bW, idx, COM_SING, ENC_SIZE);
		#ifdef _DEBUG_SMIN_PSV
	    DebugCom("(buf) E(x_i) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
		#endif

	    // DH --> CSP
		try {
			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
			CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
		}
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		// DH <-- CSP
		ulRecv.lock();
		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
		ulRecv.unlock();

		// receive E(y[i])
		short sLen ;
		sLen = Byte2Short(ucRecvPtr+HED_LEN);
		CommAmt[1] += sLen;
		#ifdef _DEBUG_SMIN_PSV
	    DebugCom("(buf) E(y[i]) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
		#endif

		#ifdef _DEBUG_Assert
	    sLen = Byte2Short(ucRecvPtr+HED_LEN);
	    assert((*(ucRecvPtr+2)==COM_SING)&&(sLen==ENC_SIZE));
		#endif

	    paillier_ciphertext_from_bytes(&cZ, ucRecvPtr+HED_SIZE, ENC_SIZE);
		#ifdef _DEBUG_SMIN_PSV
		DebugOut("E(y[i]) (Copied)\t\t\t", cZ.c, idx);
		DebugDec("E(y[i]) \t\t\t", &cZ, idx);
		#endif

		// memory initialization
		memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
		ucRecvPtr[1] = 0xff;
		//------------------------------------------------------------------------------

	    int DTcol = tDat->DTposST3[iClusterSize][row][perPos];
		#ifdef _DEBUG_SMIN_PSV
		DebugDec("cY: \t\t", &cZ, idx);
		printf("permuted position: \t\t(%02d, %02d) \n", row, DTcol);
		printf("resultant cV position: \t(%02d) \n", perPos);
	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			DebugDec("(permuted) cDT[i,j,k]: \t\t", tTbl->cDT4[row][DTcol]+k, idx);
	    }
		#endif

	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
	    	SecMul(cV[perPos]+k, &cZ, tTbl->cDT4[row][DTcol]+k, idx, tRecv, CommAmt);
	    }

		#ifdef _DEBUG_SMIN_PSV
	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			DebugDec("(permuted) cV[row][col]: \t\t", cV[perPos]+k, idx);
	    }
		#endif
	}

	//------------------------------------------------------------------------------
	// print the results of SMIN
	#ifdef _DEBUG_SMIN_PSV
	for (int j=0 ; j<iSize ; j++) {
		for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			printf("(%02d, %02d, %02d) ", row,j,k);
			DebugDec("cV[j,k]: \t\t", cV[j]+k, idx);
		}
	}
	#endif

	return;
}

// SMIN protocol for SMin^PSV (single-thread)
void PaillierCrypto::SmaxSmin(cTxt_t* cC, const int iSize, const int *DTpos,
		const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	cTxt_t *cD, cS, cQ, cW, cP[iSize] ;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cS.c, cQ.c, cW.c, NULL);
	for (int i=0 ; i<iSize ; i++)
		mpz_inits(cP[i].c, NULL);
	#else
	mpz_init2(cS.c, 2*GMP_N_SIZE*2);
	mpz_init2(cQ.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iSize ; i++)
		mpz_init2(cP[i].c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	#ifdef _DEBUG_Detail
	printf("************  (single-thread) SMAX/SMIN  ************\n");
	#endif

	for (int j=DATA_SQUARE_LENGTH-1 ; j>=0 ; j--) {

		#ifdef _DEBUG_Detail
		printf("\n---------  Round %02d  ---------\n", j);
		#endif

		// --------------------------------------------------------------------
		// Step 1: compute P[i] = SM(C[i], d[i,j])
		mpz_set_ui(cS.c, 1);		// initialize local S

		for (int i=0 ; i<iSize ; i++) {
			cD = tTbl->cDT4[row][DTpos[i]]+j ;

			#ifdef _DEBUG_SMIN
			DebugDec("d[i,j] \t\t\t", cD, idx);
			DebugDec("C[i] \t\t\t", cC+i, idx);
			#endif

			// P[i] = SM(C[i], d[i,j])
			SecMul(cP+i, cC+i, cD, idx, tRecv, CommAmt);

			// step 2: local S = Sum_(i=1~n) P[i]
			paillier_mul(mPubKey, &cS, &cS, cP+i);

			#ifdef _DEBUG_SMIN
			DebugDec("P[i] = SM(C[i], d[i,j]) \t", cP+i, idx);
			DebugDec("local S = Sum_(i=1~n) P[i] \t", &cS, idx);
			#endif
		}
		// --------------------------------------------------------------------
		// Step 3: compare q = SEQ(S, 0)
		SEQ(&cQ, &cS, &(tVal->c0), tVal, idx, tRecv, CommAmt);

		#ifdef _DEBUG_SMIN
		DebugDec("Q = \t\t\t\t", &cQ, idx);
		#endif

		// --------------------------------------------------------------------
		// Step 4: compute C[i] = SM(q, C[i]) * P[i]

		for (int i=0 ; i<iSize ; i++) {
			// w = SM(q, C[i])
			SecMul(&cW, &cQ, cC+i, idx, tRecv, CommAmt);

			// C[i] = w * P[i]
			paillier_mul(mPubKey, cC+i, cP+i, &cW);

			#ifdef _DEBUG_SMIN
			DebugDec("w = SM(q, C[i]) \t\t", &cW, idx);
			DebugDec("C[i] = w * P[i] \t\t", cC+i, idx);
			#endif
		}
	}

	// ----------------------------------------------------------------------------------------
	// print results of SMIN
	#ifdef _DEBUG_Detail
	for (int i=0 ; i<iSize ; i++) {
		printf("(%02d, %02d) ", row, DTpos[i]);
		DebugDec("C[i]: \t\t", cC+i, idx);
	}
	#endif

	return;
}

// SMin^PSV protocol for multi-thread
void PaillierCrypto::SmaxSmin_PSV(cTxt_t cV[][DATA_SQUARE_LENGTH], cTxt_t cDmin_i[][DATA_SQUARE_LENGTH],
		const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bW[ENC_SIZE]={0,};
	unsigned int *d1, *d2;
	std::condition_variable *dv;
	int iNumActThds = tDat->iNumActThds;

	int iClusterSize = tTbl->n;
	int iSizeTotal = iClusterSize-1;
	int RemTask = row % iNumActThds;
	unsigned int thSize = tDat->thNumS3RRow[iClusterSize][RemTask];
	int tmpIdx = (idx-RemTask) / (iClusterSize%THREAD1_NUM);
	int iSize = tDat->thPosS3RRow[iSizeTotal][thSize][tmpIdx+1] - tDat->thPosS3RRow[iSizeTotal][thSize][tmpIdx];
	unsigned int uiR[iSize], uiT, tmp;
	cTxt_t cS, cW, cZ, *cC, *cU;
	cTxt_t *cDmin_iRow;
	pTxt_t pR;
	int Start;
	int perPos;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cS.c, cW.c, pR.m, cZ.c, NULL);
	#else
	mpz_init2(cS.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	mpz_init2(pR.m, 2*GMP_N_SIZE*2);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	#endif

	srand((unsigned int)time(0));

	tSync->cd=0;

	// ----------------------------------------------------------------------------------------
	d1 = &(tSync->d1[RemTask]);
	d2 = &(tSync->d2[RemTask]);
	dv = &(tSync->dv[RemTask]);

	cC = tDat->cC[RemTask] ;
	cU = tDat->cU[RemTask] ;
	// ----------------------------------------------------------------------------------------

	Start = tDat->thPosS3RRow[iSizeTotal][thSize][tmpIdx] ;

	// part candidate mode: C[i] = Gamma[i]
	for (int j=0 ; j<iSize ; j++) {
		tmp = tDat->DTposST3[iClusterSize][row][Start+j];
		mpz_set(cC[Start+j].c, tTbl->cGam1[tmp].c);
	}

	#ifdef _DEBUG_SMIN_PSV
	for (int j=0 ; j<iSize ; j++) {
		printf("cDT position (%d, %d) ", row, tDat->DTposST3[iClusterSize][row][Start+j]);
		DebugDec("cC[j]: \t\t", cC+(Start+j), idx);
	}
	#endif

	SmaxSmin(cC+Start, iSize, tDat->DTposST3[iClusterSize][row], Start, tDat->cS[RemTask], &(tDat->cSum[RemTask]),
			RemTask, tmpIdx, thSize, row, tTbl, tDat, tVal, idx, tRecv, tSync, CommAmt);

	// ------------------------------------------------------------------------------------
	ulSync.lock();
	(*d1)++;
	if (*d1 < thSize)
		dv->wait(ulSync, [&] {return *d1>=thSize;});
	else {
		// --------------------------------------------------------------------
    	// print results of SMIN
		#ifdef _DEBUG_Detail
		for (int i=0 ; i<iSizeTotal ; i++) {
			printf("(%02d) ", i);
			DebugDec("C[i]: \t\t", cC+i, idx);
		}
		#endif
		// --------------------------------------------------------------------
		// E(u[i]) = sum_j=1~i-1 E(C[j])
		cU = tDat->cU[RemTask] ;
		mpz_set_ui(cU[0].c, 1);
		for (int i=1 ; i<iSizeTotal ; i++)
			paillier_mul(mPubKey, cU+i, cU+(i-1), cC+(i-1));

		#ifdef _DEBUG_SMIN_PSV
		for (int i=0 ; i<iSizeTotal ; i++)
			DebugDec("cU: \t\t", cU+i, idx);
		#endif

		// --------------------------------------------------------------------
		*d2=0;		dv->notify_all();
	}
	ulSync.unlock();
	// ------------------------------------------------------------------------------------

	// phi permutation idx / depermutation idx
	for (int i=0 ; i<iSize ; i++)
		uiR[i] = i;
	for (int i=iSize-1 ; i>0 ; i--) {
		tmp = rand() % (i+1);
		uiT = uiR[i];
		uiR[i] = uiR[tmp];
		uiR[tmp] = uiT;
	}

	#ifdef _DEBUG_SMIN_PSV
	printf("[DH-%03d] phi permutation : ", idx);
	for (int i=0 ; i<iSize ; i++)
		printf(" %d ", uiR[i]);
	std::cout << std::endl;
	#endif

	for (int i=0 ; i<iSize ; i++) {
		perPos = uiR[i];
		#ifdef _DEBUG_SMIN_PSV
		DebugDec("(permuted) cC: \t\t", cC+Start+perPos, idx);
		DebugDec("(permuted) cU: \t\t", cU+Start+perPos, idx);
		#endif

		paillier_mul(mPubKey, &cS, cU+Start+perPos, cU+Start+perPos);
		paillier_mul(mPubKey, &cS, cC+Start+perPos, &cS);
		paillier_mul(mPubKey, &cS, &(tVal->cN1), &cS);

		#ifdef _DEBUG_SMIN_PSV
			DebugDec("(permuted) cS: \t\t", &cS, idx);
		#endif

		// E(w[i]) = E(s[i])^r_i
		mpz_urandomm(pR.m, state, mPubKey->n);
		paillier_exp(mPubKey, &cW, &cS, &pR);

		#ifdef _DEBUG_SMIN_PSV
			DebugDec("cW: \t\t", &cW, idx);
		#endif

		//------------------------------------------------------------------------------
		// send w[i]'
		bW[0]=0;
		paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cW);
		SetSendMsg(ucSendPtr, bW, idx, COM_SING, ENC_SIZE);
		#ifdef _DEBUG_SMIN_PSV
	    DebugCom("(buf) E(x_i) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
		#endif

	    // DH --> CSP
		try {
			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
			CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
		}
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		// DH <-- CSP
		ulRecv.lock();
		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
		ulRecv.unlock();

		// receive E(y[i])
		short sLen ;
		sLen = Byte2Short(ucRecvPtr+HED_LEN);
		CommAmt[1] += sLen;
		#ifdef _DEBUG_SMIN_PSV
		//sLen = Byte2Short(ucRecvPtr+HED_LEN);
	    DebugCom("(buf) E(y[i]) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
		#endif

		#ifdef _DEBUG_Assert
	    sLen = Byte2Short(ucRecvPtr+HED_LEN);
	    assert((*(ucRecvPtr+2)==COM_SING)&&(sLen==ENC_SIZE));
		#endif

	    paillier_ciphertext_from_bytes(&cZ, ucRecvPtr+HED_SIZE, ENC_SIZE);
		#ifdef _DEBUG_SMIN_PSV
		DebugOut("E(y[i]) (Copied)\t\t\t", cZ.c, idx);
		DebugDec("E(y[i]) \t\t\t", &cZ, idx);
		#endif
		//------------------------------------------------------------------------------

		// memory initialization
		memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
		ucRecvPtr[1] = 0xff;

	    int DTcol = tDat->DTposST3[iClusterSize][row][Start+perPos];
	    int Vcol = Start+perPos;
		#ifdef _DEBUG_SMIN_PSV
		DebugDec("cY: \t\t", &cZ, idx);
		printf("permuted position: \t\t(%02d, %02d) \n", row, DTcol);
		printf("resultant cV position: \t(%02d) \n", Vcol);
	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			DebugDec("(permuted) cDT[i,j,k]: \t\t", tTbl->cDT4[row][DTcol]+k, idx);
	    }
		#endif

	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
	    	SecMul(cV[Vcol]+k, &cZ, tTbl->cDT4[row][DTcol]+k, idx, tRecv, CommAmt);
	    }

		#ifdef _DEBUG_SMIN_PSV
	    for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			DebugDec("(permuted) cV[row][col]: \t\t", cV[Vcol]+k, idx);
	    }
		#endif
	}

	// ------------------------------------------------------------------------------------
	ulSync.lock();
	(tSync->cd)++;
	if (tSync->cd < iNumActThds)
		tSync->cv.wait(ulSync, [&] {return tSync->cd>=iNumActThds;});
	else {
		// --------------------------------------------------------------------
    	// print results of SMIN_PSV
		#ifdef _DEBUG_Detail
    	for (int j=0 ; j<iSizeTotal ; j++) {
    		for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
    			printf("(%02d, %02d, %02d) ", row,j,k);
    			DebugDec("cV[j,k]: \t\t", cV[j]+k, idx);
    		}
    	}
		#endif
    	// --------------------------------------------------------------------
    	cDmin_iRow = tDat->cDmin_i[row];
		for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
			for (int j=0 ; j<iClusterSize-1 ; j++) {
				paillier_mul(mPubKey, cDmin_iRow+k, cDmin_iRow+k, cV[j]+k);
			}
			#ifdef _DEBUG_Step3
			DebugDec("cMin[i,k]: \t\t", cDmin_iRow+k, idx);
			#endif
		}
    	// --------------------------------------------------------------------
		tSync->cv.notify_all();
	}
	ulSync.unlock();

	return;
}

// SMIN protocol for SMin^PSV (multi-thread)
void PaillierCrypto::SmaxSmin(cTxt_t *cC, const int iSize, const int *DTpos, int Start, cTxt_t *cLocalS_th, cTxt_t *cGlobalS, const int remTask, const int tmpIdx, const unsigned int thNum,
		const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	unsigned int *d1, *d2;
	std::condition_variable *dv;
	unsigned int *pSyncCnt2, *pSyncCnt1;

	cTxt_t *cD, *cLocalS, *cQ, cW, cP[iSize] ;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cW.c, NULL);
	for (int i=0 ; i<iSize ; i++)
		mpz_inits(cP[i].c, NULL);
	#else
	mpz_init2(cW.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iSize ; i++)
		mpz_init2(cP[i].c, 2*GMP_N_SIZE*2);
	#endif

	// ----------------------------------------------------------------------------------------
	cLocalS = cLocalS_th+tmpIdx;
	cQ = &(tDat->cQ[remTask]) ;
	d1 = &(tSync->d1[remTask]);
	d2 = &(tSync->d2[remTask]);
	dv = &(tSync->dv[remTask]);
	// ----------------------------------------------------------------------------------------

	for (int j=DATA_SQUARE_LENGTH-1 ; j>=0 ; j--) {

		#ifdef _DEBUG_Detail
		printf("\n---------  Round %02d  ---------\n", j);
		#endif

		// --------------------------------------------------------------------
		// Step 1: compute P[i] = SM(C[i], d[i,j])
		mpz_set_ui(cLocalS->c, 1);		// initialize local S

		for (int i=0 ; i<iSize ; i++) {
			cD = tTbl->cDT4[row][DTpos[Start+i]]+j ;

			#ifdef _DEBUG_SMIN
			DebugDec("d[i,j] \t\t\t", cD, idx);
			DebugDec("C[i] \t\t\t", cC+i, idx);
			#endif

			// P[i] = SM(C[i], d[i,j])
			SecMul(cP+i, cC+i, cD, idx, tRecv, CommAmt);

			// step 2: local S = Sum_(i=1~n) P[i]
			paillier_mul(mPubKey, cLocalS, cLocalS, cP+i);

			#ifdef _DEBUG_SMIN
			DebugDec("P[i] = SM(C[i], d[i,j]) \t", cP+i, idx);
			DebugDec("local S = Sum_(i=1~n) P[i] \t", cLocalS, idx);
			#endif
		}
		// --------------------------------------------------------------------
		ulSync.lock();
		#if (DATA_SQUARE_LENGTH-1) % 2 == 0
		if (j%2 == 0) 		{	pSyncCnt1 = d1;	pSyncCnt2 = d2;	}
		else 					{	pSyncCnt1 = d2;	pSyncCnt2 = d1;	}
		#else
		if (j%2 == 0) 		{	pSyncCnt1 = d2;	pSyncCnt2 = d1;	}
		else 					{	pSyncCnt1 = d1;	pSyncCnt2 = d2;	}
		#endif
		(*pSyncCnt2)++;
		if (*pSyncCnt2 < thNum)
			dv->wait(ulSync, [&] {return *pSyncCnt2>=thNum;});
		else {
			// -----------------------------------------------------------------
			// Step 2: compute global S
			mpz_set_ui(cGlobalS->c, 1);		// initialize global S
			for (unsigned int i=0 ; i<thNum ; i++) {
				#ifdef _DEBUG_SMIN
				DebugDec("local S = Sum_(i=1~n) P[i] \t", cLocalS_th+i, idx);
				#endif
				paillier_mul(mPubKey, cGlobalS, cGlobalS, cLocalS_th+i);
			}

			#ifdef _DEBUG_SMIN
			DebugDec("compute global S \t", cGlobalS, idx);
			#endif
			// -----------------------------------------------------------------
			// Step 3: compare q = SEQ(S, 0)
			SEQ(cQ, cGlobalS, &(tVal->c0), tVal, idx, tRecv, CommAmt);

			#ifdef _DEBUG_SMIN
			DebugDec("Q = \t\t\t\t", cQ, idx);
			#endif
			// -----------------------------------------------------------------
			*pSyncCnt1=0;		dv->notify_all();
		}
		ulSync.unlock();
		// --------------------------------------------------------------------
		// Step 4: compute C[i] = SM(q, C[i]) * P[i]

		for (int i=0 ; i<iSize ; i++) {
			// w = SM(q, C[i])
			SecMul(&cW, cQ, cC+i, idx, tRecv, CommAmt);

			// C[i] = w * P[i]
			paillier_mul(mPubKey, cC+i, cP+i, &cW);

			#ifdef _DEBUG_SMIN
			DebugDec("w = SM(q, C[i]) \t\t", &cW, idx);
			DebugDec("C[i] = w * P[i] \t\t", cC+i, idx);
			#endif
		}
	}

	// ------------------------------------------------------------------------------------
	if (DATA_SQUARE_LENGTH%2 == 0) {
		ulSync.lock();
		(*d2)++;
		if (*d2 < thNum)
			dv->wait(ulSync, [&] {return *d2>=thNum;});
		else {
			*d1=0;
			dv->notify_all();
		}
		ulSync.unlock();
	}

	return;
}

// secure equality protocol
int PaillierCrypto::SEQ(cTxt_t* cQ, cTxt_t* cA, cTxt_t* cB, shdVal_t* tVal, unsigned short idx, recv_t* tRecv, unsigned long* CommAmt)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE*DATA_NUMBER_LENGTH]={0,};
	unsigned char bX[ENC_SIZE]={0,};

	int iDelta;
	cTxt_t cX, cR, cT, cD, cY[DATA_NUMBER_LENGTH], cC[DATA_NUMBER_LENGTH], cU[DATA_NUMBER_LENGTH];
	pTxt_t pR;
	mpz_t	mN2;
	int iRand[DATA_NUMBER_LENGTH]={0,};
	unsigned int uiR[DATA_NUMBER_LENGTH], uiT, k;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cX.c, cR.c, cT.c, cD.c, pR.m, mN2, NULL);
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
		mpz_inits(cY[i].c, cC[i].c, cU[i].c, NULL);
	#else
	mpz_init2(cX.c, 2*GMP_N_SIZE*2);
	mpz_init2(cR.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cD.c, 2*GMP_N_SIZE*2);
	mpz_init2(pR.m, 2*GMP_N_SIZE*2);
	mpz_init2(mN2, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		mpz_init2(cY[i].c, 2*GMP_N_SIZE);
		mpz_init2(cC[i].c, 2*GMP_N_SIZE);
		mpz_init2(cU[i].c, 2*GMP_N_SIZE);
	}
	#endif

	srand((unsigned int)time(0));

	// E(r), 0 <= r < M  (Pick random number)
	mpz_urandomm(pR.m, state, mPubKey->n);
	paillier_enc(&cR, mPubKey, &pR, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SEQ
	DebugOut("r \t\t\t\t", pR.m, idx);
	DebugOut("E(r)\t\t\t\t", cR.c, idx);
	#endif

	// E(x) = E(a)*E(b)^(M-1)*E(r)
	//paillier_exp(mPubKey, &c_B, cB, &(tVal->pN1));
	//paillier_mul(mPubKey, &cX, cA, &c_B);
	//paillier_mul(mPubKey, &cX, &cX, &cR);
	paillier_mul(mPubKey, &cX, cA, &cR);
	#ifdef _DEBUG_SEQ
	DebugOut("E(x) = E(a)*E(b)^(M-1)*E(r) (ciphertext)\t\t", cX.c, idx);
	DebugDec("E(x) = E(a)*E(b)^(M-1)*E(r) \t\t\t", &cX, idx);
	#endif

	//------------------------------------------------------------------------------
	// send E(x)
	paillier_ciphertext_to_bytes(bX, ENC_SIZE, &cX);
	SetSendMsg(ucSendPtr, bX, idx, COM_SEQ1, ENC_SIZE);
	#ifdef _DEBUG_SEQ
    DebugCom("(buf) E(x) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

    // (1st) DH --> CSP
	try {
		mSSocket->send(ucSendPtr, HED_SIZE+ENC_SIZE);
		CommAmt[0] += ENC_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	//-----------------------------------------------------
	// bit-decomposition of random value r
	#ifdef _DEBUG_SEQ
	DebugOut("rand = \t\t\t\t", pR.m, idx);
	#endif

	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		iRand[i] = mpz_fdiv_ui(pR.m, 2);
		mpz_fdiv_q_2exp(pR.m, pR.m, 1);
	}

	#ifdef _DEBUG_SEQ
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		printf("%d ", iRand[i]) ;
	}
	printf("\n");
	#endif
	//-----------------------------------------------------

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive E(x[n])
	short sLen ;
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_SEQ
    DebugCom("(buf) E(x[n]) (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
    sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_SEQ1)&&(sLen==ENC_SIZE*DATA_NUMBER_LENGTH));
	#endif

    for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		paillier_ciphertext_from_bytes(&cY[i], ucRecvPtr+HED_SIZE+(ENC_SIZE*i), ENC_SIZE);
		#ifdef _DEBUG_SEQ
		DebugOut("E(x[i]) (Copied)\t\t\t", cY[i].c, idx);
		DebugDec("E(x[i]) \t\t\t", cY+i, idx);
		#endif
    }

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;
	//------------------------------------------------------------------------------

	// E(y[i]) = E(x[i] xor r[i])
    for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
    	if (iRand[i] == 1) {
        	// E(y[i]) = E(x[i] xor 1)
    		paillier_exp(mPubKey, cY+i, cY+i, &(tVal->pN1));
    		paillier_mul(mPubKey, cY+i, &(tVal->c1), cY+i);
    	}
    		// E(y[i]) = E(x[i])

		#ifdef _DEBUG_SEQ
		DebugDec("E(y[i]) = E(x[i] xor r[i]) \t\t\t", cY+i, idx);
		#endif
    }

	// delta = {0, 1}
	iDelta = rand() % 2;
	#ifdef _DEBUG_SEQ
	DebugOut("Delta = \t\t\t\t\t\t\t", iDelta, idx);
	#endif

	if (iDelta==0) {
		// E(c[0]) = E( Sum_j=0~l-1 E(Y[j]) )
		mpz_set_ui(cT.c, 1);		// initialize T
		for (int j=0 ; j<DATA_NUMBER_LENGTH ; j++)
			paillier_mul(mPubKey, &cT, &cT, cY+j);
		#ifdef _DEBUG_SEQ
		DebugDec("(Alpha=0) E(c[0]) = E( Sum_j=0~l-1 E(Y[j]) ) \t", &cT, idx);
		#endif

		// E(C[0]) = E(c[0])^r
		mpz_urandomm(pR.m, state, mPubKey->n);
		paillier_exp(mPubKey, &cC[0], &cT, &pR);
		#ifdef _DEBUG_SEQ
		DebugDec("E(C[0]) = E(c[0])^r \t\t\t\t", &cC[0], idx);
		#endif

		// mN2 = N^2
		mpz_mul(mN2, mPubKey->n, mPubKey->n);

		// E(c[i]) = E(r) for j=1~l-1
		for (int i=1 ; i<DATA_NUMBER_LENGTH ; i++) {
			mpz_urandomm(cC[i].c, state, mN2);
			#ifdef _DEBUG_SEQ
			DebugDec("E(c[i]) = E(r) for j=1~l-1 \t\t\t\t", cC+i, idx);
			#endif
		}
	}
	else {
		// E(u[i]) = sum_j=i+1~l-1 E(y[j])
		mpz_set_ui(cU[DATA_NUMBER_LENGTH-1].c, 1);
		for (int j=DATA_NUMBER_LENGTH-2 ; j>=0 ; j--)
			paillier_mul(mPubKey, cU+j, cU+(j+1), cY+(j+1));

		for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
			// E(T) = E(-1) * E(Y[i]) * (Sum_j=i+1~l-1 E(Y[j]) )^2
			paillier_mul(mPubKey, &cT, cU+i, cU+i);
			paillier_mul(mPubKey, &cT, cY+i, &cT);
			paillier_mul(mPubKey, &cT, &(tVal->cN1), &cT);
			#ifdef _DEBUG_SEQ
			printf("[DH-%03d] *** (Alpha=1) [i=%d] ***\n", (int)idx, i);
			DebugDec("(Alpha=1) E(c[i]) = E(-1) * E(Y[i]) * (Sum_j=i+1~l-1 E(Y[j]))^2 \t", &cT, idx);
			#endif

			// E(c[i]) = E(c)^r
			mpz_urandomm(pR.m, state, mPubKey->n);
			paillier_exp(mPubKey, cC+i, &cT, &pR);
			#ifdef _DEBUG_SEQ
			DebugDec("E(c[i]) = E(c[i])^r \t\t\t\t", cC+i, idx);
			#endif
		}
	}

	// compute permutation phi
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
		uiR[i] = i;
	for (int i=DATA_NUMBER_LENGTH-1 ; i>0 ; i--) {
		k = rand() % (i+1);
		uiT = uiR[i];
		uiR[i] = uiR[k];
		uiR[k] = uiT;
	}

	#ifdef _DEBUG_SEQ
	printf("[DH-%03d] phi permutation : ", idx);
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
		printf(" %d ", uiR[i]);
	std::cout << std::endl;
	#endif

    // send phi( D[n] )
	#ifdef _DEBUG_SEQ
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		DebugOut("cD[i] (ciphertext)\t", cC[i].c, idx);
		DebugDec("cD[i] \t\t", cC+i, idx);
	}
	#endif

	//------------------------------------------------------------------------------
	SetSendMsg(ucSendPtr, bX, idx, COM_SEQ2, ENC_SIZE);
	for (int j=0 ; j<DATA_NUMBER_LENGTH ; j++) {
		bX[0]=0;
		paillier_ciphertext_to_bytes(bX, ENC_SIZE, cC+uiR[j]);
		memcpy(ucSendPtr+HED_SIZE+(j*ENC_SIZE), bX, ENC_SIZE);
	}
	ucSendPtr[3] = (unsigned char)((DATA_NUMBER_LENGTH*ENC_SIZE) & 0x000000ff);
	ucSendPtr[4] = (unsigned char)(((DATA_NUMBER_LENGTH*ENC_SIZE) & 0x0000ff00) >> 8);
	#ifdef _DEBUG_SEQ
	DebugCom("(buf) D (Hex):\t\t\t", ucSendPtr, HED_SIZE+(DATA_NUMBER_LENGTH*ENC_SIZE), idx);
	#endif

	// (2nd) DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
		CommAmt[0] += sizeof(ucSendPtr) - HED_SIZE;
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive E(delta_b)
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	CommAmt[1] += sLen;
	#ifdef _DEBUG_SEQ
    DebugCom("(buf) E(delta_b) (Hex)\t\t\t\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
    sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_SEQ2)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(&cD, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SEQ
	//DebugOut("delta_b' (Copied)", cD.c, idx);
	DebugDec("delta_b' \t\t", &cD, idx);
	#endif

	// memory initialization
	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;
	//------------------------------------------------------------------------------

	if (iDelta==0) {
		// E(gamma) = E(delta_b)
		mpz_set(cQ->c, cD.c);
		#ifdef _DEBUG_SEQ
		DebugDec("E(gamma) = E(delta_b) \t\t\t\t", cQ, idx);
		#endif
	}
	else {
		// E(gamma) = E(1) * E(delta_b)^(N-1)
		paillier_exp(mPubKey, &cT, &cD, &(tVal->pN1));
		paillier_mul(mPubKey, cQ, &(tVal->c1), &cT);
		#ifdef _DEBUG_SEQ
		//DebugOut("E(gamma) = E(1) * E(delta_b)^(N-1)", cD->c, idx);
		DebugDec("E(gamma) = E(1) * E(delta_b)^(N-1) \t", cQ, idx);
		#endif
	}

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
//	assert(   p0.m->_mp_alloc ==   1);
//	assert(   pR.m->_mp_alloc ==   GMP_N_SIZE);
//	assert(   cP.c->_mp_alloc == 2*GMP_N_SIZE);
//	assert(   cT.c->_mp_alloc == 2*GMP_N_SIZE+1);
//	assert(   cV.c->_mp_alloc == 2*GMP_N_SIZE*2);
//	assert(cU[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	#endif

	return 0;
}

// tasks for program termination
short PaillierCrypto::TerminatePgm(unsigned short idx, recv_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+1]={0,};
	unsigned char bA[1]={0,};
	short sRes;

	// send a'
	bA[0] = 0xFF;
	SetSendMsg(ucSendPtr, bA, idx, COM_TERM, 1);
	#ifdef _DEBUG_TERMINATE_PGM
    DebugCom("(buf) Terminate Program (sending)\t", ucSendPtr, 1+HED_SIZE, idx);
	#endif

    // DH --> CSP
	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// DH <-- CSP
	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive message
	short sLen ;
	#ifdef _DEBUG_TERMINATE_PGM
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
    DebugCom("(buf) Terminate Program (received)\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
    assert((*(ucRecvPtr+2)==COM_TERM)&&(sLen==3));
	#endif

    sRes = ucRecvPtr[HED_SIZE+2];

//	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
//	ucRecvPtr[1] = 0xff;

	return sRes;
}

// print the running time and communication amount of ppAHC
void PaillierCrypto::PrintResult(short sTrd, shdDat_t *tDat)
{
	struct tm *date;
	const time_t t = time(NULL);
	double TotalTime[4]={0,}, dTmp, MB[2];
	unsigned long TotalAmount[4][2]={0,}, uiTmp[2];
	int Min, Sec;


	printf("\n\n************************   Running Time (seconds)  ************************\n");
	printf("initialization\t\t   : %d \n", (int)tDat->dTimeEtc[TIME_INIT][0]);

	for (int i=0 ; i<DATA_NUM-PARAM_K ; i++) {
		dTmp = 0;
		for (int j=0 ; j<4 ; j++) {
			dTmp += tDat->dTimeRnd[i][j];
			TotalTime[j] += tDat->dTimeRnd[i][j];
		}

		printf("%02d round (cluster size: %02d): %d  \t ( %d, \t %d, \t %d, \t %d ) \n",
				i+1, DATA_NUM-i-1, (int)dTmp, (int)tDat->dTimeRnd[i][0], (int)tDat->dTimeRnd[i][1], (int)tDat->dTimeRnd[i][2], (int)tDat->dTimeRnd[i][3]);
	}

	printf("\n\n************************   Communication Amount (bytes)  ************************\n");
	printf("initialization\t\t   : %lu / %lu \n", tDat->uiCommAmtIni[0], tDat->uiCommAmtIni[1]);

	for (int i=0 ; i<DATA_NUM-PARAM_K ; i++) {
		uiTmp[0] = uiTmp[1] = 0;
		for (int j=0 ; j<4 ; j++) {
			uiTmp[0] += tDat->uiCommAmt[i][j][0];
			uiTmp[1] += tDat->uiCommAmt[i][j][1];
			TotalAmount[j][0] += tDat->uiCommAmt[i][j][0];
			TotalAmount[j][1] += tDat->uiCommAmt[i][j][1];
		}

		printf("%02d round (cluster size: %02d): %lu / %lu  \t ( %lu / %lu, \t %lu / %lu, \t %lu / %lu, \t %lu / %lu ) \n",
				i+1, DATA_NUM-i-1, uiTmp[0], uiTmp[1],
				tDat->uiCommAmt[i][0][0], tDat->uiCommAmt[i][0][1], tDat->uiCommAmt[i][1][0], tDat->uiCommAmt[i][1][1],
				tDat->uiCommAmt[i][2][0], tDat->uiCommAmt[i][2][1], tDat->uiCommAmt[i][3][0], tDat->uiCommAmt[i][3][1]);
	}

	date = localtime(&t);
	dTmp = TotalTime[0] + TotalTime[1] + TotalTime[2] + TotalTime[3];
	Min = (int)(dTmp/60);
	Sec = (int)dTmp%60;
	uiTmp[0] = TotalAmount[0][0] + TotalAmount[1][0] + TotalAmount[2][0] + TotalAmount[3][0];
	uiTmp[1] = TotalAmount[0][1] + TotalAmount[1][1] + TotalAmount[2][1] + TotalAmount[3][1];
	MB[0] = (double)uiTmp[0]/1048576;
	MB[1] = (double)uiTmp[1]/1048576;

	printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
	printf( "[%d/%d/%d %d:%d:%d] %d min %d sec (%d sec) \t ( %d, \t %d, \t %d, \t %d ) \n" ,
			date->tm_year + 1900 , date->tm_mon + 1 , date->tm_mday , date->tm_hour , date->tm_min , date->tm_sec,
			Min, Sec, (int)dTmp, (int)TotalTime[0], (int)TotalTime[1], (int)TotalTime[2], (int)TotalTime[3]);
	printf( "[%d/%d/%d %d:%d:%d] %.2f / %.2f MB (%lu / %lu B) \t ( %lu / %lu, \t %lu / %lu, \t %lu / %lu, \t %lu / %lu )\n" ,
			date->tm_year + 1900 , date->tm_mon + 1 , date->tm_mday , date->tm_hour , date->tm_min , date->tm_sec,
			MB[0], MB[1], uiTmp[0], uiTmp[1], TotalAmount[0][0], TotalAmount[0][1], TotalAmount[1][0], TotalAmount[1][1],
			TotalAmount[2][0], TotalAmount[2][1], TotalAmount[3][0], TotalAmount[3][1]);
	printf( "\t\t    (key size, rounds, N, threads) = ( %d , %d , %d , %d / %d )\n", MOD_SIZE, DATA_NUM-PARAM_K, DATA_NUM, THREAD1_NUM, sTrd);

	return ;
}

// set sending message
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

// set sending message : two data & two length
inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data1, unsigned char* Data2,
		const unsigned short Idx, const unsigned char Tag, const unsigned short Len1, const unsigned short Len2)
{
	ucSendPtr[0] = (unsigned char)(Idx & 0x000000ff);
	ucSendPtr[1] = (unsigned char)((Idx & 0x0000ff00) >> 8);
	ucSendPtr[2] = Tag;
	ucSendPtr[3] = (unsigned char)((Len1+Len2) & 0x000000ff);
	ucSendPtr[4] = (unsigned char)(((Len1+Len2) & 0x0000ff00) >> 8);
	memcpy(ucSendPtr+HED_SIZE, Data1, Len1);
	memcpy(ucSendPtr+HED_SIZE+Len1, Data2, Len2);

	return;
}

// set sending message
inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, const unsigned char* ucRecvPtr, unsigned char bData)
{
    memcpy(ucSendPtr, ucRecvPtr, HED_LEN);
    ucSendPtr[HED_LEN]  = 1;
    ucSendPtr[HED_LEN+1]= 0;
    ucSendPtr[HED_SIZE] = bData;

	return;
}

inline unsigned short PaillierCrypto::Byte2Short(unsigned char* In)
{
	return (In[1]<<8) + In[0];
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const mpz_t pData, const int idx)
{
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, pData);
	return;
}

void PaillierCrypto::DebugOut(const char* pMsg, short pData, const int idx)
{
	printf("[DH-%03d] %s : %d\n", idx, pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugDec(const char* pMsg, cTxt_t* cTxt, const int idx)
{
	pTxt_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif

	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 2*GMP_N_SIZE+1);
	#endif

    return;
}

inline void PaillierCrypto::DebugDec(cTxt_t *cTxt, const int idx)
{
	pTxt_t m;
	int iBit, iTotal=0;

	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif

	for (int i=0 ; i<DATA_SQUARE_LENGTH ; i++) {
		paillier_dec(&m, mPubKey, mPrvKey, cTxt+i);
		iBit = 1 - mpz_get_ui(m.m);
		iTotal += iBit<<i;
	}
	printf("[DH-%03d] %X -\t", idx, iTotal);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <=2*GMP_N_SIZE+1);
	#endif

    return;
}

inline void PaillierCrypto::DebugDecBit(cTxt_t* cTxt)
{
	pTxt_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif

	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("%ZX", &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 2*GMP_N_SIZE+1);
	#endif

    return;
}

inline void PaillierCrypto::DebugCom(const char* pMsg, unsigned char *data, int len, const int idx)
{
	constexpr char hexmap[] = {'0', '1', '2', '3', '4', '5', '6', '7',
							   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	char s[len*2+1]={0,};
	for (int i = 0; i < len; ++i) {
		s[2*i]	 = hexmap[(data[i] & 0xF0) >> 4];
		s[2*i+1] = hexmap[data[i] & 0x0F];
	}
	printf("[DH-%03d] %s : %s\n", idx, pMsg, s);

	return;
}

// print DT
void PaillierCrypto::DebugDT(const char* pMsg, shdTbl_t *tTbl, int pLast, int iStep, const int idx)
{
	int iClusterSize = tTbl->n + pLast;
	int iDT[iClusterSize][iClusterSize]={0,}, iBit, DTElement;
	pTxt_t pT;

	#ifdef _DEBUG_INIT_1
	mpz_init(pT.m);
	#else
	mpz_init2(pT.m, 3*GMP_N_SIZE);
	#endif

	for (int i=0 ; i<iClusterSize ; i++) {
		for (int j=0 ; j<iClusterSize ; j++) {
			DTElement=0;
			if (i < j) {
				for (int k=0 ; k<DATA_SQUARE_LENGTH ; k++) {
					if (iStep == 3)
						paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cDT3[i][j]+k);
					else if (iStep == 4)
						paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cDT4[i][j]+k);
					else
						assert(0);
					iBit = 1 - mpz_get_ui(pT.m);
					DTElement += iBit<<k;
				}
			}
			iDT[i][j] = DTElement;
		}
	}

	printf("[DH-%03d] (DT) %s \n", idx, pMsg);
	for (int i=0 ; i<iClusterSize ; i++) {
		printf("%d -\t", i+1);
		for (int j=0 ; j<iClusterSize ; j++) {
			printf("%X \t", iDT[i][j]);
		}
		printf("\n");
	}

	#ifdef _DEBUG_Assert
	assert(pT.m->_mp_alloc <= 3*GMP_N_SIZE);
	#endif

    return;
}

// print MT
void PaillierCrypto::DebugMT(const char* pMsg, shdTbl_t *tTbl, const int idx)
{
	int iClusterSize = tTbl->n;
	int iBeta[iClusterSize][iClusterSize]={0,}, iGamma2[iClusterSize][2]={0,}, iGamma1[iClusterSize]={0,};
	pTxt_t pT;

	#ifdef _DEBUG_INIT_1
	mpz_init(pT.m);
	#else
	mpz_init2(pT.m, 3*GMP_N_SIZE);
	#endif

	for (int i=0 ; i<iClusterSize ; i++) {
		for (int j=0 ; j<iClusterSize ; j++) {
			paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cBeta[i]+j);
			iBeta[i][j] = mpz_get_ui(pT.m);
		}

		paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cGam2[i]);
		iGamma2[i][0] = mpz_get_ui(pT.m);
		paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cGam2[i]+1);
		iGamma2[i][1] = mpz_get_ui(pT.m);

		paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cGam1+i);
		iGamma1[i] = mpz_get_ui(pT.m);
	}

	printf("[DH-%03d] (MT) %s \n", idx, pMsg);

	printf("***  beta-table  *** \n");
	for (int i=0 ; i<iClusterSize ; i++) {
		printf("%d -\t", i+1);
		for (int j=0 ; j<iClusterSize ; j++)
			printf("%X \t", iBeta[i][j]);
		printf("\n");
	}

	printf("***  gamma-table (%d) *** \n", iClusterSize);
	for (int j=0 ; j<2 ; j++) {
		for (int i=0 ; i<iClusterSize ; i++)
			printf("\t%X", iGamma2[i][j]);
		printf("\n");
	}

	printf("***  gamma-array (%d) *** \n", iClusterSize);
	for (int i=0 ; i<iClusterSize ; i++)
		printf("\t%X", iGamma1[i]);
	printf("\n");

	#ifdef _DEBUG_Assert
	assert(pT.m->_mp_alloc <= 3*GMP_N_SIZE);
	#endif

    return;
}

// print CT
void PaillierCrypto::DebugCT(const char* pMsg, shdTbl_t *tTbl, int pLast, const int idx)
{
	int iClusterSize = tTbl->n + pLast;
	int iCT[iClusterSize][DATA_NUM]={0,};
	pTxt_t pT;

	#ifdef _DEBUG_INIT_1
	mpz_init(pT.m);
	#else
	mpz_init2(pT.m, 3*GMP_N_SIZE);
	#endif

	for (int i=0 ; i<iClusterSize ; i++) {
		for (int j=0 ; j<DATA_NUM ; j++) {
			paillier_dec(&pT, mPubKey, mPrvKey, tTbl->cCT[i]+j);
			iCT[i][j] = mpz_get_ui(pT.m);
		}
	}

	printf("[DH-%03d] (CT) %s \n", idx, pMsg);
	for (int i=0 ; i<iClusterSize ; i++) {
		printf("%d -\t", i+1);
		for (int j=0 ; j<DATA_NUM ; j++) {
			printf("%X \t", iCT[i][j]);
		}
		printf("\n");
	}

	#ifdef _DEBUG_Assert
	assert(pT.m->_mp_alloc <= 3*GMP_N_SIZE);
	#endif

    return;
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

int PaillierCrypto::Receiver(recv_t* tRecv, shdDat_t *tDat)
{
	unsigned char* ucRecvBuf[THREAD1_NUM+1];
	unsigned short sLen ;
	short sIdx, sbx=0;
	int iReceivedLen;

	for (int i=0 ; i<THREAD1_NUM+1 ; i++) {
		ucRecvBuf[i] = new unsigned char[MAXRECV];
		memset(ucRecvBuf[i], 0, MAXRECV);
		ucRecvBuf[i][1] = 0xff;
	}

	while (1) {
		while (ucRecvBuf[sbx][1]!=0xff) {
			sbx++;
			if (sbx>=tDat->iNumActThds+1)
				sbx=0;
		}
		iReceivedLen = 0;

	    // DH --> CSP
		try {
			iReceivedLen = mSSocket->recv(ucRecvBuf[sbx], HED_SIZE);
			sLen = Byte2Short(ucRecvBuf[sbx]+HED_LEN);
			while (iReceivedLen != HED_SIZE+sLen)
				iReceivedLen += mSSocket->recv(ucRecvBuf[sbx]+iReceivedLen, HED_SIZE+sLen-iReceivedLen);
	    }
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		sIdx = Byte2Short(ucRecvBuf[sbx]);

		#ifdef _DEBUG_Communication
		DebugCom("(Receiver Thread) RecvBuf (Hex)", ucRecvBuf[sbx], HED_SIZE+sLen, sIdx);
		#endif

		//-----------------------------------------------------
		if (*(ucRecvBuf[sbx]+2) == COM_TERM) {
			assert(Byte2Short(ucRecvBuf[sbx]+HED_SIZE) == MOD_SIZE);
			memcpy(ucKillMsg, ucRecvBuf[sbx], HED_SIZE+3);
		    tRecv->m.lock();
		    tRecv->pa[sIdx] = ucKillMsg;
		    tRecv->cv.notify_all();
		    tRecv->m.unlock();

			#ifdef _DEBUG_Communication
			DebugCom("(Receiver Thread) RecvBuf (Hex)", ucKillMsg, HED_SIZE+3, sIdx);
			#endif

			for (int i=0 ; i<THREAD1_NUM+1 ; i++)
				delete[] ucRecvBuf[i];
			return 0;
		}
		//-----------------------------------------------------

		tRecv->m.lock();
	    tRecv->pa[sIdx] = ucRecvBuf[sbx];
	    tRecv->cv.notify_all();
	    tRecv->m.unlock();
	}

	return 0;
}
