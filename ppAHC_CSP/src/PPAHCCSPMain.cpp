//============================================================================
// Name        : ppAHC_CSP_Main.cpp
// Author      : PJS
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include "PaillierCrypto.h"
#include <vector>
#include <thread>
using namespace std;

void Receiver_Trd(PaillierCrypto *oCrypto, recv_t* tRecv)
{
	oCrypto->Receiver(tRecv);
	return;
}

// CSP thread for ppAHC protocol
void ppAHC_Trd(PaillierCrypto *oCrypto, unsigned short idx, recv_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr;

	#ifdef _DEBUG_THREAD
	printf("<<<  %02d - th Thread start  >>>\n", idx);
	#endif

	while(1) {
		ulRecv.lock();
		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
		ulRecv.unlock();

		// according to the instruction of received message, select the corresponding function
		switch(ucRecvPtr[HED_INS]) {
		case COM_MUL1  :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Sqare Command >>>\n\n", idx);
			#endif
			oCrypto->SecMul1	 (idx, ucRecvPtr);			break;
		case COM_MUL2  :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Multiplication Command >>>\n\n", idx);
			#endif
			oCrypto->SecMul2	 (idx, ucRecvPtr);			break;
		case COM_LSB   :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Encrypted LSB Command >>>\n\n", idx);
			#endif
			oCrypto->EncryptedLSB(idx, ucRecvPtr);			break;
		case COM_SVR   :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - SVR Command >>>\n\n", idx);
			#endif
			oCrypto->SVR		 (idx, ucRecvPtr);			break;
		case COM_SING :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Single Output mode Command >>>\n\n", idx);
			#endif
			oCrypto->SingleOutput(idx, ucRecvPtr);			break;
		case COM_SEQ1 :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - SEQ1 Command >>>\n\n", idx);
			#endif
			oCrypto->SEQ1		 (idx, ucRecvPtr);			break;
		case COM_SEQ2 :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - SEQ2 Command >>>\n\n", idx);
			#endif
			oCrypto->SEQ2		 (idx, ucRecvPtr);			break;
		case COM_STP2 :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Step 2 Command >>>\n\n", idx);
			#endif
			oCrypto->Step_2	 (idx, ucRecvPtr, tRecv);	break;
		case COM_ST2D1 :
		case COM_ST2D2 :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Step 2 (delta) Command >>>\n\n", idx);
			#endif
			oCrypto->Step_2_Delta(idx, ucRecvPtr);			break;
		case COM_ST2V :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - Step 2 (v) Command >>>\n\n", idx);
			#endif
			oCrypto->Step_2_v	 (idx, ucRecvPtr);			break;
		case COM_TERM :
			#ifdef _DEBUG_THREAD
			printf("\n<<<  %03d - TERMINATE Command >>>\n\n", idx);
			#endif
			oCrypto->TerminatePgm(idx, ucRecvPtr);
			printf("\n<<<  %02d - th Thread END !!!  >>>\n", idx);
			return;
		default:
			short sIdx = (ucRecvPtr[1]<<8) + ucRecvPtr[0];
			oCrypto->DebugComMain("RecvBuf Error Msg (Hex)\t", ucRecvPtr, HED_SIZE+2*ENC_SIZE, idx, sIdx);
			assert(0);
		}
	}

	return;
}

int main() {
	vector<thread> ppAHC_Threads, Receiver_Threads;
	recv_t tRecv;
	PaillierCrypto oCrypto(MOD_SIZE);

	oCrypto.distributePubKey();			// distribute public key
	oCrypto.distributePrvKey();			// distribute private key
	for (int i=0 ; i<THREAD_NUM ; i++)
		tRecv.pa[i] = NULL;

	// generate ppAHC threads and Communication threads
	for (int i=0 ; i<THREAD_NUM ; i++)
		ppAHC_Threads.push_back(thread(ppAHC_Trd, &oCrypto, (unsigned short)i, &tRecv));
	Receiver_Threads.push_back(thread(Receiver_Trd, &oCrypto, &tRecv));

	// terminate ppAHC threads and Communication threads
	for (int i=0 ; i<THREAD_NUM ; i++)
		ppAHC_Threads[i].join();
	Receiver_Threads[0].join();
	std::cout << std::endl << "<<<  (ppAHC_CSP) Receiving Thread TERMINATE !!!  >>>"<< std::endl << std::endl;

	return 0;
}

