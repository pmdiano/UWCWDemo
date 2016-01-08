#include "uwcw.h"
#include "DllUWCW.h"
#include <Windows.h>
#include <time.h>

extern svmData svm;
extern UWCW_FLOAT FFT_DR;
extern bool theta_first_frame;
extern UWCW_FLOAT theta_previous;
extern UWCW_FLOAT rhoScore, rhoScorePrev;
extern UWCW_INT rhoNum, rhoNumPrev;
extern cableLine *rhoCableLine;
extern cableLine *rhoCablePrev;
extern bool useTrackingAlgorithm;



UWCW_FLOAT testDll()
{
	return (UWCW_FLOAT)svm.M;
}

void InitUWCW()
{
	cleanUpData();
}

void SetAlgorithm(bool useOriginal)
{
	useTrackingAlgorithm = !useOriginal;
	cleanUpData();
}

void getScore(UWCW_FLOAT *currentScore, UWCW_FLOAT *prevScore,
		UWCW_FLOAT *minRange, UWCW_FLOAT *maxRange)
{
	UWCW_INT i, midPoint;
	UWCW_FLOAT range;
	*currentScore = rhoScore;
	*prevScore = rhoScorePrev;

	*minRange = 0;
	*maxRange = 0;

	if (rhoNum > 0)
	{
		midPoint = rhoCableLine[0].nPix / 2;
		if (midPoint > 0)
		{
			range = (UWCW_FLOAT)(rhoCableLine[0].range[midPoint]) * FFT_DR;
			*minRange = range;
			*maxRange = range;
		}

		for (i=1; i<rhoNum; i++)
		{
			midPoint = rhoCableLine[i].nPix / 2;
			if (midPoint > 0)
			{
				range = (UWCW_FLOAT)(rhoCableLine[i].range[midPoint]) * FFT_DR;
				if (range < *minRange) *minRange = range;
				if (range > *maxRange) *maxRange = range;
			}
		}
	}
}


UWCW_INT detectUWCW(UWCW_FLOAT *RP_dBm, UWCW_INT M, UWCW_INT N,
		UWCW_FLOAT *Ang_plot, UWCW_FLOAT Ang_endpt, UWCW_FLOAT fft_dR,
		UWCW_FLOAT MagMax, UWCW_FLOAT MagMin, bool isFirstFrame,
		UWCW_INT *rangeIndices, UWCW_INT *angleIndices)
{
	if (isFirstFrame)
	{
		srand(time(NULL));
	}

	theta_first_frame = isFirstFrame;
	singleFrameDetect(RP_dBm, M, N, Ang_plot, Ang_endpt, MagMax, MagMin, &svm, 0, M, fft_dR);
	if (rhoScore >= 0.5)
		return catCableIndices(rhoCableLine, rhoNum, 1, rangeIndices, angleIndices);
	else
		return 0;
}

UWCW_INT catCableIndices(cableLine *cables, UWCW_INT numCable, UWCW_FLOAT score,
	UWCW_INT *rangeIndices, UWCW_INT *angleIndices)
{
	UWCW_INT totalNum = 0;
	UWCW_INT i, j, k;
	UWCW_INT *rangeDst = rangeIndices;
	UWCW_INT *angleDst = angleIndices;

	if (score < 0.5)
		return 0;

	for (i=0; i<numCable; i++)
	{
		if (cables[i].isCable)
		{
			totalNum += cables[i].nPix;

			memcpy(rangeDst, cables[i].range, cables[i].nPix * sizeof(UWCW_INT));
			memcpy(angleDst, cables[i].angle, cables[i].nPix * sizeof(UWCW_INT));

			rangeDst += cables[i].nPix;
			angleDst += cables[i].nPix;
		}
	}

	return totalNum;
}

// ********************************************************************************** //
// entry point of dll, not part of core implementation
// ********************************************************************************** //
BOOL WINAPI DllMain(
    HINSTANCE hinstDLL,  // handle to DLL module
    DWORD fdwReason,     // reason for calling function
    LPVOID lpReserved )  // reserved
{
    // Perform actions based on the reason for calling.
    switch( fdwReason ) 
    { 
        case DLL_PROCESS_ATTACH:
         // Initialize once for each new process.
         // Return FALSE to fail DLL load.
			svmRead("uwcw_svm.dat", &svm);
            break;

        case DLL_THREAD_ATTACH:
         // Do thread-specific initialization.
            break;

        case DLL_THREAD_DETACH:
         // Do thread-specific cleanup.
            break;

        case DLL_PROCESS_DETACH:
         // Perform any necessary cleanup.
            break;
    }
    return TRUE;  // Successful DLL_PROCESS_ATTACH.
}