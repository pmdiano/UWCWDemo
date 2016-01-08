// ********************************************************************************** //
// UWCW core implementation
// Qirong Ma
// 11/1/11
// ********************************************************************************** //


#include "DefUWCW.h"


// ********************************************************************************** //
// data structure definition
// ********************************************************************************** //
// svm data structure
// ********************************************************************************** //
typedef struct _svmData
{
	unsigned short M;		// # of support vectors
	unsigned short N;		// feature dimension
	float **sv;
	float *alpha;
	float bias;
	float *shift;
	float *scale;
	float *svss;	// sum of squares of sv
	float *b;		// numel: 2N+1
} svmData;

// ********************************************************************************** //
// cable line data structure
// ********************************************************************************** //
typedef struct _cableLine
{
	bool isCable;	// svm result
	UWCW_FLOAT score;	// linear combination score
	UWCW_FLOAT rho;		// rho and theta of BWcart
	UWCW_FLOAT theta;
	UWCW_INT nPix;		// number of pixels on this line
	UWCW_FLOAT *data;	// data on this line
	UWCW_INT *range;
	UWCW_INT *angle;
} cableLine;

// ********************************************************************************** //
// theta particle filter data structure
// ********************************************************************************** //
typedef struct _tparticle
{
	UWCW_FLOAT theta;
	UWCW_FLOAT thetap;
	UWCW_FLOAT thetapp;
	UWCW_FLOAT w;
} tparticle;

// ********************************************************************************** //
// rho particle filter and tracker data structures
// ********************************************************************************** //
enum tracker_status{
	head,		// head for the linked list
	candidate,
	healthy,
	dormant,
	dead
};

typedef struct _rparticle
{
	UWCW_FLOAT rho;
	UWCW_FLOAT rhop;
	UWCW_FLOAT rhopp;
	UWCW_FLOAT theta;
	UWCW_FLOAT thetap;
	UWCW_FLOAT thetapp;
	UWCW_FLOAT w;
	UWCW_FLOAT hVal;
	UWCW_FLOAT svmVal;
	UWCW_FLOAT linearVal;
} rparticle;

typedef struct _rtracker
{
	rparticle *rpts_;	// particles in this tracker
	UWCW_INT rpts_num;	// num of particles in this tracker

	tracker_status rstatus;
	UWCW_INT life_time;
	UWCW_INT status_time;

	// tracked rho of current time
	UWCW_FLOAT rho;		
	UWCW_FLOAT rhop;
	UWCW_FLOAT rhopp;
	UWCW_FLOAT theta;
	UWCW_FLOAT thetap;
	UWCW_FLOAT thetapp;
	
	// previous tracked data
	UWCW_FLOAT *data_original_;
	UWCW_INT *data_threshold_;
	UWCW_INT data_num;

	// previous tracked Hough data
	UWCW_INT *hough_;
	UWCW_INT hough_num;	// number of hough points saved; be an odd number

	// to make the list of rtrackers
	struct _rtracker *next_;
} rtracker;



// ********************************************************************************** //
// internal implementation functions
// ********************************************************************************** //
// high-level processing functions
// ********************************************************************************** //
UWCW_FLOAT singleFrameDetect(UWCW_FLOAT *RP_dBm, UWCW_INT M, UWCW_INT N,
	UWCW_FLOAT *Ang_plot, UWCW_FLOAT Ang_endpt, UWCW_FLOAT MagMax, UWCW_FLOAT MagMin,
	svmData *svm, UWCW_INT y1, UWCW_INT y2, UWCW_FLOAT fft_dR);

void completeCableData(cableLine *cables, UWCW_INT num,
	UWCW_FLOAT *RP_dBm, UWCW_FLOAT *Ang_plot, UWCW_INT M, UWCW_INT N,
	UWCW_INT *BW, UWCW_INT x0);

UWCW_FLOAT doubleFrameScore(cableLine *current, UWCW_INT currentNum,
	cableLine *prev, UWCW_INT prevNum, UWCW_FLOAT prevScore, bool firstFrame);

void cleanUpData();
void adjustParameters(UWCW_INT M, UWCW_INT N);
void freeLineArray(cableLine *cLs, UWCW_INT num);

// ********************************************************************************** //
// double frame score functions
// ********************************************************************************** //
UWCW_INT numCable(cableLine *cl, UWCW_INT num);

UWCW_INT checkParallel(cableLine *cl, UWCW_INT num);

UWCW_FLOAT lineCorrelate(cableLine l1, cableLine l2,
	UWCW_FLOAT Rsigma, UWCW_FLOAT Tsigma, UWCW_FLOAT C);

UWCW_INT temporalCorrelate(cableLine *current, UWCW_INT currentNum,
	cableLine *prev, UWCW_INT prevNum);


// ********************************************************************************** //
// SVM related functions
// ********************************************************************************** //
void svmRead(const char* filename, svmData *svm);
UWCW_FLOAT svmCalc(svmData *svm, UWCW_FLOAT *f);
UWCW_FLOAT linearCalc(svmData *svm, UWCW_FLOAT *f);

// ********************************************************************************** //
// feature computation and getting data functions
// ********************************************************************************** //
UWCW_INT getLineData(UWCW_FLOAT *Bscope, UWCW_FLOAT *Ang_plot, UWCW_INT M, UWCW_INT N,
	UWCW_INT *BWBscope,	UWCW_INT x0, UWCW_FLOAT rho, UWCW_FLOAT theta,
	UWCW_FLOAT *data, UWCW_INT *range, UWCW_INT *angle,	bool makeLong);

void powerSpectrum(UWCW_FLOAT *input, UWCW_INT size, UWCW_FLOAT *spec);

void autoCorrelation(UWCW_FLOAT *input, UWCW_INT size, UWCW_FLOAT *autoCo);

UWCW_FLOAT sumArray(UWCW_FLOAT *data, UWCW_INT n1, UWCW_INT n2);

UWCW_INT compareThreshold(UWCW_FLOAT *data, UWCW_INT dataLen, UWCW_FLOAT T,
	UWCW_INT *num, UWCW_INT *loc);

UWCW_FLOAT uwStd(UWCW_FLOAT *data, UWCW_INT dataLen);

void thresholdFeature(UWCW_FLOAT *data, UWCW_INT size, UWCW_FLOAT T,
	UWCW_FLOAT *f1, UWCW_FLOAT *f2, UWCW_FLOAT *f3);

void computeFeature(UWCW_FLOAT *data, UWCW_INT dataSize, UWCW_INT specSize,
	UWCW_INT autoSize, UWCW_FLOAT *feature, UWCW_INT featSize);

// ********************************************************************************** //
// thresholding functions
// ********************************************************************************** //
void blockThreshold(UWCW_FLOAT *Bscope,	UWCW_INT M,	UWCW_INT N, UWCW_INT *BWBscope,
	UWCW_INT sliceHeight, UWCW_INT sliceWidth, UWCW_FLOAT percentage,
	UWCW_FLOAT threshGate, UWCW_INT y1, UWCW_INT y2);

UWCW_FLOAT blockPercent(UWCW_FLOAT *Bscope,	UWCW_INT M,  UWCW_INT N,
	UWCW_INT h1, UWCW_INT h2, UWCW_INT w1, UWCW_INT w2, UWCW_FLOAT percentage,
	UWCW_FLOAT *T2,	UWCW_FLOAT *p2);

int dcomp(const void *a, const void *b);

// ********************************************************************************** //
// coordinate transform functions
// ********************************************************************************** //
void coorTransform(UWCW_INT *BWBscope, UWCW_INT M, UWCW_INT N, UWCW_INT *BWcart,
	UWCW_INT *pNC, UWCW_INT *px0, UWCW_FLOAT Ang_endpt, UWCW_FLOAT *Ang_plot,
	UWCW_INT interFac, UWCW_INT y1, UWCW_INT y2);

void coorTransAssign(UWCW_INT *BWcart, UWCW_INT MC, UWCW_INT NC, UWCW_FLOAT rho,
	UWCW_FLOAT theta, UWCW_INT x0, UWCW_INT val);

// ********************************************************************************** //
// Hough Transform functions
// ********************************************************************************** //
void houghTransform(UWCW_INT *BWcart, UWCW_INT MC, UWCW_INT NC, UWCW_INT y1, UWCW_INT y2,
	UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_FLOAT *theta, UWCW_INT thetaLen,	UWCW_INT *H);

// ********************************************************************************** //
// theta particle filtering functions
// ********************************************************************************** //
UWCW_FLOAT entropyMax3(UWCW_INT *h, UWCW_INT rhoLen, UWCW_INT ss);

tparticle* theta_initialize_particles(UWCW_INT num, UWCW_FLOAT theta,
	UWCW_INT* H, UWCW_INT rhoLen, UWCW_INT thetaLen, UWCW_INT idx);

void theta_free_particles(tparticle* particles, UWCW_INT num);

UWCW_FLOAT theta_average(tparticle* theta_particles, UWCW_INT num);

void theta_transition_particles(tparticle* particles, UWCW_INT num, UWCW_FLOAT theta_std);

void gauss(UWCW_FLOAT std1, UWCW_FLOAT std2, UWCW_FLOAT *y1, UWCW_FLOAT *y2);

void theta_compute_likelihood(tparticle* theta_particles, UWCW_INT num,
	UWCW_INT* H, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT(*similarity)(UWCW_INT *h1, UWCW_INT *h2, UWCW_INT n));

UWCW_INT thetaToIdx(UWCW_FLOAT theta, UWCW_INT thetaEnd, UWCW_INT thetaLen);

tparticle* tresample(tparticle* ptcs, UWCW_INT n);

int tparticleCmp(const void *p1, const void *p2);

UWCW_FLOAT similarity_1(UWCW_INT *h1, UWCW_INT *h2, UWCW_INT n);

UWCW_FLOAT trackTheta( UWCW_INT *BWcart, UWCW_INT MC, UWCW_INT NC,
	UWCW_FLOAT *theta, UWCW_INT thetaLen,
	UWCW_INT **Hout, UWCW_FLOAT **rhoOut, UWCW_INT *rhoLenOut);

// ********************************************************************************** //
// rho particle filtering functions
// ********************************************************************************** //
void getMaxVal(	UWCW_INT *H, UWCW_INT M, UWCW_INT N,
	UWCW_FLOAT *maxVal, UWCW_INT *maxRow, UWCW_INT *maxCol);

UWCW_FLOAT get_rho_thresh(UWCW_INT *h, UWCW_INT rhoLen);
UWCW_FLOAT getRhoThreshRaw(UWCW_FLOAT *RP_dBm, 
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, 
	UWCW_FLOAT MagMax, UWCW_FLOAT MagMin);

void rho_transition_particles(rparticle *rho_particles, UWCW_INT num,
	UWCW_FLOAT rho_std, UWCW_FLOAT theta_new, UWCW_INT ss, rtracker *r_track);

UWCW_FLOAT rho_compute_likelihood(rparticle *rho_particles, UWCW_INT num, UWCW_FLOAT rhop,
	UWCW_FLOAT theta, UWCW_INT *h, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0, svmData *svm, rtracker *r_track);

UWCW_INT rho_to_idx(UWCW_FLOAT rho, UWCW_FLOAT *rhos, UWCW_INT rhoLen);

void rho_normalize_weights(rparticle *rpts, UWCW_INT num);

int rparticleCmp(const void *p1, const void *p2);

rparticle* rresample(rparticle* ptcs, UWCW_INT n);

void rtracker_free(rtracker *rtracker_);

// ********************************************************************************** //
// rho tracker list functions
// ********************************************************************************** //
void process_rtrackers(rtracker *rtracker_,
	UWCW_FLOAT theta, UWCW_INT *h, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0, svmData *svm, UWCW_INT margin);

rtracker *rtracker_init_zero();

rparticle* rho_initialize_particles(UWCW_INT num, UWCW_FLOAT rho, UWCW_FLOAT theta,
	UWCW_FLOAT hVal, UWCW_FLOAT svmVal, UWCW_FLOAT linearVal);

rtracker *rtracker_init(UWCW_FLOAT rho, UWCW_FLOAT theta,
	UWCW_FLOAT *lineData, UWCW_INT *lineRange, UWCW_INT *lineAngle, UWCW_INT dataNum,
	UWCW_INT *BW, UWCW_INT M, UWCW_INT N, UWCW_FLOAT hVal, UWCW_FLOAT svmVal,
	UWCW_FLOAT linearVal, UWCW_INT *hough_data, UWCW_INT hnum);

UWCW_INT rtracker_free_list(rtracker *head);

UWCW_INT get_rtracker_num(rtracker *rtracker_);

UWCW_INT add_rtrackers(rtracker **rtrackers__, UWCW_INT *h, UWCW_INT rhoLen, UWCW_FLOAT theta,
	UWCW_FLOAT *rho, UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT x0, svmData *svm, UWCW_INT margin, UWCW_FLOAT thresh);

UWCW_INT rhoTrackerToCable(rtracker *rtrackers_, UWCW_INT totalNum, cableLine **cls);

UWCW_INT get_rtracker_rho(rtracker *rtracker_, UWCW_FLOAT *rho, UWCW_INT *is_cable);

UWCW_INT rhoTrack(UWCW_FLOAT theta, UWCW_INT *H, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot, UWCW_FLOAT MagMax,
	UWCW_FLOAT MagMin, UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0,
	svmData *svm, UWCW_INT margin);

// ********************************************************************************** //
// training functions
// ********************************************************************************** //
void saveCableData(cableLine *cls, UWCW_INT num);
void saveLineData(UWCW_FLOAT *lineData, UWCW_INT dataNum);


// ********************************************************************************** //
// previous algorithm functions
// ********************************************************************************** //
UWCW_INT houghPeaks(UWCW_INT *H,			// Hough data
				UWCW_INT rhoLen,			// # of rows in H
				UWCW_INT thetaLen,			// # of columns in H
				UWCW_INT numPeaks,			// # of peaks to find
				UWCW_FLOAT hGate,			// gate value for a peak
				UWCW_INT nHoodRowHalf,		// neighborhood size half row
				UWCW_INT nHoodColHalf,		// neighborhood size half column
				UWCW_INT *peaks);

UWCW_INT detectLines( UWCW_INT *BWcart,		// input b-n-w cartesian image
				 UWCW_INT MC,				// size of BWcart,
				 UWCW_INT NC,
				 UWCW_INT ystart,			// the region to detect
				 UWCW_INT yend,
				 UWCW_FLOAT *theta,			// HT parameter
				 UWCW_INT thetaLen,
				 UWCW_INT sliceHeight,		// slice-processing height
				 UWCW_INT sliceOverlap,		// overlapping between slices
				 UWCW_INT numPeaks,			// # of peaks in each slice
				 UWCW_FLOAT *lines,
				 UWCW_FLOAT ratio);

UWCW_INT removeRepLines(UWCW_FLOAT *cLs, UWCW_INT num);

void detectCablesOld(UWCW_FLOAT *RP_dBm, UWCW_INT M, UWCW_INT N,
	UWCW_FLOAT *Ang_plot, svmData *svm, UWCW_INT y1, UWCW_INT y2,
	UWCW_INT *BW, UWCW_INT *BWcart, UWCW_INT NC, UWCW_INT x0);

void voteTheta(cableLine *cables, UWCW_INT cableNum);

UWCW_INT thetaToIdx2(UWCW_FLOAT theta);