// ********************************************************************************** //
// external interface functions, to be called in C# GUI
// ********************************************************************************** //
extern "C"
{
	__declspec(dllexport) UWCW_FLOAT testDll();


	__declspec(dllexport) UWCW_INT detectUWCW(UWCW_FLOAT *RP_dBm, UWCW_INT M, UWCW_INT N,
		UWCW_FLOAT *Ang_plot, UWCW_FLOAT Ang_endpt, UWCW_FLOAT fft_dR,
		UWCW_FLOAT MagMax, UWCW_FLOAT MagMin, bool isFirstFrame,
		UWCW_INT *rangeIndices, UWCW_INT *angleIndices);

	__declspec(dllexport) void InitUWCW();

	__declspec(dllexport) void getScore(UWCW_FLOAT *currentScore, UWCW_FLOAT *prevScore,
		UWCW_FLOAT *minRange, UWCW_FLOAT *maxRange);

	__declspec(dllexport) void SetAlgorithm(bool useOriginal);
}

UWCW_INT catCableIndices(cableLine *cables, UWCW_INT numCable, UWCW_FLOAT score,
	UWCW_INT *rangeIndices, UWCW_INT *angleIndices);