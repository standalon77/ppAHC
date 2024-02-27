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
//#define _DEBUG_Setting
//#define _DEBUG_Detail
//#define _DEBUG_Initialization
//#define _DEBUG_SM
//#define _DEBUG_SBD
//#define _DEBUG_LSB
//#define _DEBUG_SVR
//#define _DEBUG_Step1
//#define _DEBUG_Step2
//#define _DEBUG_Step3
//#define _DEBUG_Step4
//#define _DEBUG_SMIN_ASB
//#define _DEBUG_SMIN_PSV
//#define _DEBUG_SMIN
//#define _DEBUG_SEQ
//#define _DEBUG_TERMINATE_PGM
//#define _DEBUG_Communication

#include <gmp.h>
extern "C"{
	#include "paillier.h"
}
#include "ServerSocket.h"
#include <assert.h>
#include <queue>
#include <mutex>
#include <condition_variable>

#include <thread>
#include <vector>

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

// structure for sync thread
typedef struct {
	std::mutex m;
	std::condition_variable cv;
	unsigned int c1;
	unsigned int c2;
	std::condition_variable dv[THREAD1_NUM];
	unsigned int d1[THREAD1_NUM];
	unsigned int d2[THREAD1_NUM];
	unsigned int cd;
} sync_t;

// structure for recv thread
typedef struct {
	std::mutex m;
	std::condition_variable cv;
	unsigned char *pa[THREAD1_NUM];
} recv_t;

// structure for tables (CT, DT, MT)
typedef struct {
	int 	n;		// cluster size
	cTxt_t cCT[DATA_NUM][DATA_NUM];
	cTxt_t cDT3[DATA_NUM][DATA_NUM][DATA_SQUARE_LENGTH];
	cTxt_t cDT4[DATA_NUM][DATA_NUM][DATA_SQUARE_LENGTH];
	cTxt_t cBeta[DATA_NUM][DATA_NUM];
	cTxt_t cGam2[DATA_NUM][2];
	cTxt_t cGam1[DATA_NUM];
} shdTbl_t;

// Running Time Tag
const unsigned char TIME_INIT 	= 0 ;

// shared Data Tag
const int NUM_THD_TASK = 4;						// Thread Task의 개수
const unsigned char POS_TRI 	= 0 ;
const unsigned char POS_SQU 	= 1 ;
const unsigned char POS_LIN 	= 2 ;
const unsigned char POS_ST3 	= 3 ;

typedef struct {
	int 	thSrtPos[NUM_THD_TASK][DATA_NUM+1][THREAD1_NUM+1][2];
	// starting position of each thread's task
	// [i] data type    (4: 0~3)
	// 		[POS_TRI]				(triang) n*(n-1)/2 	(row, column)
	// 		[POS_SQU]				(square) n*n 			(row, column)
	// 		[POS_LIN][j][k][0]	(linear) n*(n-1)/2 	(index)
	// 		[POS_LIN][j][k][1]	(linear) n 			(index)
	//		[POS_ST3][j][k][0]	(Step 3) n row 		(index)
	//		[POS_ST3][j][k][1]	(Step 3) n row 		(row number assigned to each thread for remaining rows of step 3)
	// [j] number of clusters
	// [k] thread index

	int 	thNumS3RRow[DATA_NUM+1][THREAD1_NUM-1];
	// (SMIN^PSV execution for remaining row of step 3) number of threads assigned to each row
	// [i] number of clusters
	// [j] remaining row number

	int 	thPosS3RRow[DATA_NUM][THREAD1_NUM+1][THREAD1_NUM+1];
	// (SMIN^PSV execution for remaining row of step 3) starting position of thread's task (data)
	// [i] number of data for each row
	// [j] number of computing threads in each row
	// [k] thread index

	int 	DTposST1[DATA_NUM+1][DTRI_SIZ][2];
	// (SMIN^ASB) data index of SMIN -> coordinate of DT
	// [i] data size (1~N), [j] data size of SMIN (1~N*(N-1)/2), [k] (x, y) coordinate

	int		DTposST3[DATA_NUM+1][DATA_NUM][DATA_NUM-1];
	// (SMIN^PSV) data coordinate of SMIN -> coordinate of DT
	// [i] data size (1~N), [j] row number (0~(N-1)), [k] data position of SMIN (1~(N-1))

	// initialization
	cTxt_t cD [DATA_NUM][DATA_DIM];					// input dataset
	cTxt_t _cD[DATA_NUM][DATA_DIM];					// -cD[][]

	// Step 2
	cTxt_t cGamS[DATA_NUM-1];						// summation array of gamma_(i,2)

	// Step 4
	cTxt_t cBs[DATA_NUM-1][DATA_NUM-1];			// b^s-table
	cTxt_t cGs[DATA_NUM-2];							// g^s-array

	// SMIN
	cTxt_t cC[THREAD1_NUM][DTRI_SIZ];				// SMIN result
	cTxt_t cQ[THREAD1_NUM];							// SEQ result
	cTxt_t cS[THREAD1_NUM][THREAD1_NUM];			// local sum of P
	cTxt_t cSum[THREAD1_NUM];						// sum of local s

	// SMIN^mode
	cTxt_t cU[THREAD1_NUM][DTRI_SIZ];					// summation of C[i]
	cTxt_t cMDT[DATA_NUM][DATA_NUM];					// result of SMIN^ASB
	cTxt_t cDmin_ij[THREAD1_NUM][DATA_NUM-1][DATA_SQUARE_LENGTH]; 	// result of SMIN^PSV
	cTxt_t cDmin_i				[DATA_NUM]	 [DATA_SQUARE_LENGTH];	// distance of new cluster

	// computation of running time
	time_t start;										// staring point
	time_t end;										// ending point
	double dTimeRnd[DATA_NUM-PARAM_K][4];			// running time for each round
	double dTimeEtc[DATA_NUM-PARAM_K][2];			// running time for initialization

	// computation of communication amount
	unsigned long uiCommAmt[DATA_NUM-PARAM_K][4][2]={0,};	// communication amount for each round (DH/CSP)
	unsigned long uiCommAmtIni[2]={0,};						// communication amount for initialization (DH/CSP)

	// number of active threads
	int iNumActThds;
} shdDat_t;

// frequently used values
typedef struct {
	pTxt_t pN1;
	cTxt_t	c0;
	cTxt_t	c1;
	cTxt_t	cN1;
} shdVal_t;

class PaillierCrypto {
private:
	gmp_randstate_t 		state;
	ServerSocket 			*mSSocket;
	paillier_prvkey_t		*mPrvKey;
	paillier_pubkey_t		*mPubKey;
	unsigned char 		ucKillMsg[HED_SIZE+3];

	inline void SetSendMsg(unsigned char *ucSendPtr, unsigned char *Data, 						 const unsigned short Idx, const unsigned char Tag, const unsigned short Len);
	inline void SetSendMsg(unsigned char *ucSendPtr, unsigned char *Data1, unsigned char *Data2, const unsigned short Idx, const unsigned char Tag, const unsigned short Len1, const unsigned short Len2);
	inline void SetSendMsg(unsigned char *ucSendPtr, const unsigned char *ucRecvPtr, unsigned char bData);
	inline unsigned short Byte2Short(unsigned char *In);
	inline void DebugDecBit					(cTxt_t *cTxt);
	inline void DebugDec(const char *pMsg, cTxt_t *cTxt,		  						const int idx);
	inline void DebugDec						(cTxt_t* cTxt, 							const int idx);
	inline void DebugOut(const char *pMsg, const mpz_t pData, 						const int idx);
	inline void DebugOut(const char *pMsg, short pData, 		  						const int idx);
	inline void DebugCom(const char *pMsg, unsigned char *data, int len, 			const int idx);
			void DebugDT (const char *pMsg, shdTbl_t *tTbl, int pLast, int step,	const int idx);
			void DebugMT (const char *pMsg, shdTbl_t *tTbl, 						 	const int idx);
			void DebugCT (const char *pMsg, shdTbl_t *tTbl, int pLast, 			 	const int idx);
	inline void paillier_ciphertext_from_bytes(cTxt_t *ct, void *c, int len);
	inline void paillier_ciphertext_to_bytes(unsigned char *buf, int len, cTxt_t *ct );

	void SecMul		(cTxt_t *cRes, cTxt_t *cEa, cTxt_t *cEb, 							  unsigned short idx, recv_t *tRecv, unsigned long* CommAmt);
	void SecMul		(cTxt_t *cRes, cTxt_t *cEa, 			  							  unsigned short idx, recv_t *tRecv, unsigned long* CommAmt);
	void SBD			(cTxt_t *cXb,  cTxt_t *cXc, cTxt_t *cX, int iLen,shdVal_t *tVal, unsigned short idx, recv_t *tRecv, unsigned long* CommAmt);
	void EncryptedLSB	(cTxt_t *cRes, cTxt_t *cT, 							shdVal_t *tVal, unsigned short idx, recv_t *tRecv, unsigned long* CommAmt);
	int  SVR			(cTxt_t *cEX,  cTxt_t *cEXi, int len, 				shdVal_t *tVal, unsigned short idx, recv_t *tRecv, unsigned long* CommAmt);
	inline void ComputeTaskAmount(int *SizeEachTask, 			 int SizeTotalTask);
	inline void ComputeTaskAmount(int *SizeEachTask, int thSize, int SizeTotalTask);
	void SmaxSmin_ASB(cTxt_t cMDT[][DATA_NUM], 									 		 shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt);
	void SmaxSmin_PSV(cTxt_t cV[][DATA_SQUARE_LENGTH], 					const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv				  , unsigned long* CommAmt);
	void SmaxSmin_PSV(cTxt_t cV[][DATA_SQUARE_LENGTH], cTxt_t cDmin_i[][DATA_SQUARE_LENGTH],
																				const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt);
	void SmaxSmin 	(cTxt_t *cC, const int iSize, const int DTpos[][2], 	int Start, cTxt_t *cLocalS_th, cTxt_t *cGlobalS,
																								 shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt);
	void SmaxSmin 	(cTxt_t *cC, const int iSize, const int *DTpos,	const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv				  , unsigned long* CommAmt);
	void SmaxSmin		(cTxt_t *cC, const int iSize, const int *DTpos,		int Start, cTxt_t *cLocalS_th, cTxt_t *cGlobalS, const int remTask, const int tmpIdx, const unsigned int thNum,
																				const int row, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync, unsigned long* CommAmt);
	int SEQ			(cTxt_t *cQ, cTxt_t *cA, cTxt_t *cB, 																	 shdVal_t *tVal, unsigned short idx, recv_t *tRecv				  , unsigned long* CommAmt);

public:
	PaillierCrypto();
	virtual ~PaillierCrypto();

	int  Receiver		(recv_t *tRecv, shdDat_t *tDat);
	void assignTask	(shdDat_t *tDat);
	void computeDTpos	(shdDat_t *tDat);
	void SetPubKey();
	void SetPrvKey();			// for debugging
	void InitializeShdDat(sync_t *tSync, recv_t *tRecv, shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal);
	bool readDataset		(cTxt_t cD[DATA_DIM][DATA_DIM]);

	void initialization	(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync);
	void step_1			(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync);
	void step_2			(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync);
	void step_3			(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync);
	void step_4			(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t *tVal, unsigned short idx, recv_t *tRecv, sync_t *tSync);
	short TerminatePgm														   (unsigned short idx, recv_t *tRecv);
	void PrintResult		(short sTrd, shdDat_t *tDat);
};

#endif /* PAILLIERCRYPTO_H_ */
