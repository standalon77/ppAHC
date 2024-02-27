//============================================================================
// Name        : ppAHC_DH_Main.cpp
// Author      : PJS
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include "PaillierCrypto.h"
#include <time.h>
#include <vector>
#include <thread>
using namespace std;

void Receiver_Trd(PaillierCrypto *oCrypto, recv_t* tRecv, shdDat_t *tDat)
{
	oCrypto->Receiver(tRecv, tDat);
	sleep(1);

	return;
}

// DH thread for ppAHC protocol
void ppAHC_Trd(shdTbl_t *tTbl, shdDat_t *tDat, shdVal_t* tVal, unsigned short idx, sync_t* tSync, recv_t* tRecv, PaillierCrypto *oCrypto)
{
	oCrypto->initialization(tTbl, tDat, tVal, idx, tRecv, tSync);

	do {
		printf("<<<  [%02d-thread] Round Start !!! (cluster size: %d)  >>>\n", idx, tTbl->n);

		oCrypto->step_1(tTbl, tDat, tVal, idx, tRecv, tSync);
		oCrypto->step_2(tTbl, tDat, tVal, idx, tRecv, tSync);
		oCrypto->step_3(tTbl, tDat, tVal, idx, tRecv, tSync);
		oCrypto->step_4(tTbl, tDat, tVal, idx, tRecv, tSync);

		// if number of threads > number of data, terminate the last thread
		if (idx >= tTbl->n) {
			printf("\n<<<  %02d - th Thread END !!!  >>>\n", idx);
			return ;
		}
	} while (tTbl->n > PARAM_K);

	// print the running time and the communication amount
	std::unique_lock<std::mutex> ulSync(tSync->m, std::defer_lock);
	ulSync.lock();
	(tSync->c1)++;
	if (tSync->c1 >= tDat->iNumActThds) {
    	short sTrd = oCrypto->TerminatePgm(idx, tRecv);
    	oCrypto->PrintResult(sTrd, tDat);
    	tSync->cv.notify_all();
    }
	ulSync.unlock();

	printf("\n<<<  %02d - th Thread END !!!  >>>\n", idx);

	return;
}

int main() {
	vector<thread> ppAHC_Threads, Receiver_Threads;
	sync_t 	*tSync = new sync_t;
	recv_t  	*tRecv = new recv_t;
	PaillierCrypto *oCrypto = new PaillierCrypto;
	shdTbl_t	*tTbl = new shdTbl_t;
	shdDat_t	*tDat = new shdDat_t;
	shdVal_t	*tVal = new shdVal_t;

	oCrypto->assignTask(tDat);				// assign tasks to each thread
	oCrypto->computeDTpos(tDat);			// change the data index of SMIN to the coordinate of DT
	oCrypto->SetPubKey();					// get public key from CSP
	oCrypto->SetPrvKey();					// get private key from CSP for debugging
	oCrypto->InitializeShdDat(tSync, tRecv, tTbl, tDat, tVal);	// initialize shared data
	oCrypto->readDataset(tDat->cD);			// input dataset

	// generate ppAHC threads and Communication threads
	for (int i=0 ; i<THREAD1_NUM ; i++)
		ppAHC_Threads.push_back(thread(ppAHC_Trd, tTbl, tDat, tVal, (unsigned short)i, tSync, tRecv, oCrypto));
	Receiver_Threads.push_back(thread(Receiver_Trd, oCrypto, tRecv, tDat));

	// terminate ppAHC threads and Communication threads
	for (int i=0 ; i<THREAD1_NUM ; i++)
		ppAHC_Threads[i].join();
	Receiver_Threads[0].join();
	std::cout << std::endl << "<<<  (ppAHC_DH) Receiving Thread TERMINATE !!!  >>>"<< std::endl << std::endl;

	delete tSync;
	delete tRecv;
	delete oCrypto;
	delete tTbl;
	delete tDat;
	delete tVal;

	return 0;
}
