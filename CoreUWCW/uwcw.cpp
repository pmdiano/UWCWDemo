#include "uwcw.h"
#include "ParametersUWCW.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fftw3.h"
#include <time.h>

// ********************************************************************************** //
// global variables
// ********************************************************************************** //
// SVM classifier coefficients, directly loaded from "svm_data.dat"
// ********************************************************************************** //
svmData svm;
UWCW_FLOAT FFT_DR;

// ********************************************************************************** //
// theta tracking variables
// ********************************************************************************** //
bool theta_first_frame = false;
tparticle *theta_particles = NULL;
UWCW_FLOAT theta_previous = -100;

// ********************************************************************************** //
// rho tracking variables
// ********************************************************************************** //
rtracker *g_rtrackers_ = NULL;
// these are supposed to be rho tracking variables, but now be used as general variables
UWCW_FLOAT rhoScore(0), rhoScorePrev(0);
UWCW_INT rhoNum(0), rhoNumPrev(0);
cableLine *rhoCableLine = NULL;
cableLine *rhoCablePrev = NULL;
UWCW_FLOAT PREV_RHO_THRESH = 0;

bool useTrackingAlgorithm = true;

UWCW_FLOAT singleFrameDetect(UWCW_FLOAT *RP_dBm,
					  UWCW_INT M,			// # of rows of RP_dBm
					  UWCW_INT N,			// # of columns of RP_dBm
					  UWCW_FLOAT *Ang_plot,	// angle values
					  UWCW_FLOAT Ang_endpt,	// angle value end
					  UWCW_FLOAT MagMax, UWCW_FLOAT MagMin,
					  svmData *svm,
					  UWCW_INT y1,
					  UWCW_INT y2,
					  UWCW_FLOAT fft_dR)
{
	UWCW_INT interFactor;
	UWCW_INT NC, x0, rhoLen, margin, *BW(NULL), *BWcart(NULL), *H(NULL);
	UWCW_FLOAT *rho(NULL);
	// 0. copy previous frame's result
	FFT_DR = fft_dR;
	if (!theta_first_frame)
	{		
		freeLineArray(rhoCablePrev, rhoNumPrev);
		rhoCablePrev = rhoCableLine;
		rhoCableLine = NULL;
		rhoScorePrev = rhoScore;
		rhoNumPrev = rhoNum;
	}
	else
	{
		freeLineArray(rhoCableLine, rhoNum);
		freeLineArray(rhoCablePrev, rhoNumPrev);
		rhoCableLine = NULL;
		rhoCablePrev = NULL;
	}

	// 1. thresholding
	BW = (UWCW_INT*)calloc((int)M*(int)N, sizeof(UWCW_INT));
	adjustParameters(M, N);
	blockThreshold(RP_dBm, M, N, BW, SLICE_HEIGHT, SLICE_WIDTH, TH_PCT, MagMin, y1, y2);
	// 2. coordinate transform
	interFactor = (UWCW_INT)((UWCW_FLOAT)M/512.0+0.5);	
	UWCW_FLOAT x0D = M*sin(Ang_endpt*UW_PI/180);
	x0 = (UWCW_INT)(x0D)+1;
	NC = 2*x0+1;
	BWcart = (UWCW_INT*)calloc((int)M*(int)NC, sizeof(UWCW_INT));

	coorTransform(BW, M, N, BWcart, &NC, &x0, Ang_endpt, Ang_plot, interFactor, y1, y2);

	// 3. theta tracking
	margin = 8;
	if (M>2000) margin = 15;
	if (M>4000) margin = 30;
	margin = (UWCW_INT)((UWCW_FLOAT)margin * 0.286 / fft_dR);

	if (useTrackingAlgorithm)
	{
		UWCW_FLOAT maxTheta = trackTheta(BWcart, M, NC, PTHETA, PTHETA_LEN, &H, &rho, &rhoLen);

		// 4. rho tracking, also detect the cables
		rhoTrack(maxTheta, H, rho, rhoLen, PTHETA_LEN, RP_dBm, BW, BWcart,
			Ang_plot, MagMax, MagMin, M, N, NC, x0, svm, margin);		
	}
	else
	{
		detectCablesOld(RP_dBm, M, N, Ang_plot, svm, y1, y2, BW, BWcart, NC, x0);
	}

	// 5. double frame score & update result
	rhoScore = doubleFrameScore(rhoCableLine, rhoNum, rhoCablePrev, rhoNumPrev, rhoScorePrev, theta_first_frame);
	
	// training: save cable data
	//saveCableData(rhoCableLine, rhoNum);

	// clean up
	if (BW) free(BW);
	if (BWcart) free(BWcart);
	if (H) free(H);

	return 0;
}


// calculate the frame score from detection results of two frames
UWCW_FLOAT doubleFrameScore(cableLine *current, UWCW_INT currentNum,
						cableLine *prev, UWCW_INT prevNum, UWCW_FLOAT prevScore,
						bool firstFrame)	// to indicate whether this is the first frame
{
	static UWCW_FLOAT AN = 0.0;
	static UWCW_FLOAT A = 0.5;
	static UWCW_FLOAT BOOST = 1.5;
	static UWCW_FLOAT DECLINE = 1.0/BOOST;
	static UWCW_FLOAT THETA_BOOST[5] = {1.0, 1.2, 1.3, 1.4, 1.5};
	static UWCW_FLOAT CORRE_BOOST[5] = {1.2, 1.3, 1.4, 1.5, 1.5};
	UWCW_FLOAT s;
	UWCW_INT nParallel;
	UWCW_INT nCorrelated;

	// intital score based on # of cables
	s = AN;
	if (numCable(current, currentNum))
	{
		s = A;
		// parallel check
		nParallel = checkParallel(current, currentNum);
		if (nParallel >= 5)
			s = s*BOOST;
		else
			s = s*THETA_BOOST[nParallel];
	}

	// factor previous score
	if (!firstFrame)
	{
		s = (s+prevScore)/2.0;
	}

	// target correlation
	if (!firstFrame)
	{
		nCorrelated = temporalCorrelate(current, currentNum, prev, prevNum);
		if (nCorrelated == 0)
			s = s*DECLINE;
		else if (nCorrelated >= 5)
			s = s*BOOST;
		else
			s = s*CORRE_BOOST[nCorrelated-1];
	}

	// clip
	if (s>1) s=1.0;

	return s;
}

void cleanUpData()
{
	theta_previous = -100;
	rhoScore = 0;
	rhoScorePrev = 0;
	
	freeLineArray(rhoCableLine, rhoNum);
	freeLineArray(rhoCablePrev, rhoNumPrev);
	rhoCableLine = NULL;
	rhoCablePrev = NULL;
	rhoNum = 0;
	rhoNumPrev = 0;

	theta_free_particles(theta_particles, THETA_NUM_PARTICLES);
	rtracker_free_list(g_rtrackers_);
	theta_particles = NULL;
	g_rtrackers_ = NULL;
}

void freeLineArray(cableLine *cLs, UWCW_INT num)
{
	UWCW_INT i;
	if (cLs)
	{	
		for (i=0; i<num; i++)
		{
			if (cLs[i].data) free(cLs[i].data);
			if (cLs[i].angle) free(cLs[i].angle);
			if (cLs[i].range) free(cLs[i].range);
		}
	}

	if (cLs) free(cLs);
}

void svmRead(const char* filename, svmData *svm)
{
	FILE *fp;
	short i,j, M, N;
	float *buf = NULL;

	fp = fopen(filename, "rb");
	if (!fp)
	{
		fprintf(stderr, "SVM read error\n");
		exit(1);
	}		

	fread(&(svm->M), sizeof(unsigned short), 1, fp);
	fread(&(svm->N), sizeof(unsigned short), 1, fp);
	M = svm->M;
	N = svm->N;

	// allocate space
	svm->sv = (float**)calloc(sizeof(float*), M);
	for (i=0; i<M; i++)
		(svm->sv)[i] = (float*)calloc(sizeof(float), N);
	svm->alpha = (float*)calloc(sizeof(float), M);
	svm->shift = (float*)calloc(sizeof(float), N);
	svm->scale = (float*)calloc(sizeof(float), N);
	svm->b =(float*)calloc(sizeof(float), 2*N+1);

	// read sv
	buf = (float*)calloc(sizeof(float), M*N);
	fread(buf, sizeof(float), M*N, fp);
	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			svm->sv[i][j] = buf[j*M+i];
		}
	}

	// alpha, bias, shift, scale, and b
	fread(svm->alpha, sizeof(float), M, fp);
	fread(&(svm->bias), sizeof(float), 1, fp);
	fread(svm->shift, sizeof(float), N, fp);
	fread(svm->scale, sizeof(float), N, fp);
	fread(svm->b, sizeof(float), 2*N+1, fp);

	// get sum of squares of sv
	svm->svss = (float*)calloc(sizeof(float), M);
	for (i=0; i<M; i++)
	{
		svm->svss[i] = 0.0;
		for (j=0; j<N; j++)
		{
			svm->svss[i] += (svm->sv[i][j])*(svm->sv[i][j]);
		}
	}

	if (buf) free(buf);
}

UWCW_FLOAT svmCalc(svmData *svm, UWCW_FLOAT *f)
{
	UWCW_FLOAT sigma, a, val;
	UWCW_FLOAT *s;
	UWCW_INT M = svm->M;
	UWCW_INT N = svm->N;
	sigma = 1.0;
	
	// change the feature
	UWCW_INT i, j;
	for (i=0; i<svm->N; i++)
	{
		f[i] += svm->shift[i];
		f[i] *= svm->scale[i];
	}

	// 
	s = (UWCW_FLOAT*)calloc(sizeof(UWCW_FLOAT), svm->M);
	for (i=0; i<M; i++)
	{
		s[i] = 0;
		for (j=0; j<N; j++)
		{
			s[i] += svm->sv[i][j]*f[j];
		}
		s[i] *= -2.0;
	}
	a = 0;
	for (j=0; j<N; j++)
	{
		a += f[j]*f[j];
	}
	for (i=0; i<M; i++)
	{
		s[i] += svm->svss[i];
		s[i] += a;
		s[i] = exp(-1/(2*sigma*sigma)*s[i]);
	}

	//
	val = 0;
	for (i=0; i<M; i++)
	{
		val += svm->alpha[i] * s[i];
	}

	val += svm->bias;

	if (s) free(s);

	return val;
}

UWCW_FLOAT linearCalc(svmData *svm, UWCW_FLOAT *f)
{
	UWCW_FLOAT val;
	UWCW_INT i;
	val = svm->b[0];

	for (i=0; i<svm->N; i++)
	{
		val += f[i]*svm->b[i+1];
	}

	for (i=0; i<svm->N; i++)
	{
		val += f[i]*f[i]*svm->b[i+1+svm->N];
	}

	return val;
}

UWCW_INT getLineData(UWCW_FLOAT *Bscope, UWCW_FLOAT *Ang_plot, UWCW_INT M, UWCW_INT N,
	UWCW_INT *BWBscope,	UWCW_INT x0, UWCW_FLOAT rho, UWCW_FLOAT theta,
	UWCW_FLOAT *data, UWCW_INT *range, UWCW_INT *angle,	bool makeLong)
{
	UWCW_INT num;	// number of points on that line
	UWCW_INT T1, T2;	// first and last non-zero point theta index
	UWCW_INT Ridx, Tidx, i;	
	UWCW_FLOAT R, x0cost, *sint;

	if (theta>=90) theta-= 180;

	x0cost = (UWCW_FLOAT)x0 * cos(theta*UW_PI/180);
	sint = (UWCW_FLOAT*)calloc(N, sizeof(UWCW_FLOAT));
	for (i=0; i<N; i++)
	{
		sint[i] = sin((Ang_plot[i]+theta)*UW_PI/180);
	}

	if (BWBscope)
	{
		// first non-zero point
		for (T1=1; T1<N-1; T1++)
		{
			R = (rho-x0cost)/sint[T1];
			Ridx = (UWCW_INT)(R+0.5);
			if (Ridx<1 || Ridx>M) continue;
			if (BWBscope[(Ridx-1)*N+T1]) break;
		}
		
		// last non-zero point
		for (T2=N-2; T2>=1; T2--)
		{
			R = (rho-x0cost)/sint[T2];
			Ridx = (UWCW_INT)(R+0.5);
			if (Ridx<1 || Ridx>M) continue;
			if (BWBscope[(Ridx-1)*N+T2]) break;
		}
		
		// swap if necessary
		if (T1>T2)
		{
			Tidx=T1;
			T1=T2;
			T2=Tidx;
		}

		// make long
		if (makeLong)
		{
			if (T1 > CABLE_TIP_END)
				T1 = CABLE_TIP_END;
			if (N > 300)
			{
				if (T2 < 350-CABLE_TIP_END-1)
					T2 = 350-CABLE_TIP_END-1;
			}
			else
			{
				if (T2 < 176-CABLE_TIP_END-1)
					T2 = 176-CABLE_TIP_END-1;
			}
		}
	}
	else
	{
		T1 = CABLE_TIP_END;
		if (N > 300)
			T2 = 350-CABLE_TIP_END-1;
		else
			T2 = 176-CABLE_TIP_END-1;
	}

	// retrieve data
	num = 0;
	for (Tidx=T1; Tidx<=T2; Tidx++)
	{
		R = (rho-x0cost)/sint[Tidx];
		Ridx = (UWCW_INT)(R-0.5);
		if (Ridx<0 || Ridx>=M) continue;

		range[num] = Ridx;
		angle[num] = Tidx;
		data[num]  = Bscope[Ridx*N+Tidx];
		num++;
	}

	if (sint) free(sint);
	return num;
}

void powerSpectrum(UWCW_FLOAT *input,		// input data
				   UWCW_INT size,			// data size (which is also DFT size)
				   UWCW_FLOAT *spec)
{
	fftw_complex *data, *fftRes;
	fftw_plan planFwd;
	UWCW_INT i;

	data	= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
	fftRes	= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);

	planFwd = fftw_plan_dft_1d(size, data, fftRes, FFTW_FORWARD, FFTW_ESTIMATE);

	// populate input data
	for (i=0; i<size; i++)
	{
		data[i][0] = input[i];
		data[i][1] = 0.0;
	}

	fftw_execute(planFwd);

	// get spectrum
	for (i=0; i<size; i++)
	{
		spec[i] = fftRes[i][0]*fftRes[i][0] + fftRes[i][1]*fftRes[i][1];
	}

	// free memory
	fftw_destroy_plan(planFwd);

	fftw_free(data);
	fftw_free(fftRes);
}

void autoCorrelation(UWCW_FLOAT *input,
					 UWCW_INT size,
					 UWCW_FLOAT *autoCo)
{
	fftw_complex *data, *fftRes, *ifftRes;
	fftw_plan planFwd, planBwd;
	UWCW_INT i;
	UWCW_FLOAT a1;

	data	= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
	fftRes	= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
	ifftRes	= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
	
	planFwd = fftw_plan_dft_1d(size, data, fftRes,  FFTW_FORWARD,  FFTW_ESTIMATE);
	planBwd = fftw_plan_dft_1d(size, data, ifftRes, FFTW_BACKWARD, FFTW_ESTIMATE);

	// input data
	for (i=0; i<size; i++)
	{
		data[i][0] = input[i];
		data[i][1] = 0.0;
	}

	fftw_execute(planFwd);

	// compute abs square
	for (i=0; i<size; i++)
	{
		data[i][0] = fftRes[i][0]*fftRes[i][0] + fftRes[i][1]*fftRes[i][1];
		data[i][1] = 0.0;
	}

	fftw_execute(planBwd);

	// auto correlation
	for (i=0; i<size; i++)
	{
		autoCo[i] = ifftRes[i][0];
	}

	a1 = autoCo[0];
	if (a1!=0.0)
	{
		for (i=0; i<size; i++)
		{
			autoCo[i] = autoCo[i] / a1;
		}
	}
	
	fftw_destroy_plan(planFwd);
	fftw_destroy_plan(planBwd);

	fftw_free(data);
	fftw_free(fftRes);
	fftw_free(ifftRes);
}

UWCW_FLOAT sumArray(UWCW_FLOAT *data, UWCW_INT n1, UWCW_INT n2)
{
	UWCW_FLOAT s = 0.0;
	UWCW_INT i;
	for (i=n1; i<n2; i++)
	{
		s += data[i];
	}

	return s;
}

UWCW_INT compareThreshold(UWCW_FLOAT *data,
					 UWCW_INT dataLen,
					 UWCW_FLOAT T,
					 UWCW_INT *num,
					 UWCW_INT *loc)
{
	bool newSeg = 0;
	UWCW_INT locLen = 0;
	UWCW_INT i;
	loc[locLen++] = 0;
	num[0] = 0;

	for (i=0; i<dataLen; i++)
	{
		if (data[i]>=T)
		{
			if (!newSeg)
			{
				newSeg=1;
				num[0]++;
				loc[locLen++]=i+1;
			}
		}
		else
		{
			if (newSeg)
			{
				newSeg=0;
				loc[locLen++]=i+1;
			}
		}
	}

	return locLen;
}

UWCW_FLOAT uwStd(UWCW_FLOAT *data, UWCW_INT dataLen)
{
	UWCW_INT i;
	UWCW_FLOAT m, std;
	std = 0;
	if (dataLen <=1 )
	{
		return 0;
	}
	else
	{
		m = sumArray(data, 0, dataLen) / (UWCW_FLOAT)dataLen;
		for (i=0; i<dataLen; i++)
		{
			std += (data[i]-m)*(data[i]-m);
		}
		std /= (UWCW_FLOAT)(dataLen-1);
		std = sqrt(std);
		return std;
	}
}

void thresholdFeature(UWCW_FLOAT *data,		// data and size
					  UWCW_INT size,
					  UWCW_FLOAT T,			// threshold
					  UWCW_FLOAT *f1,		// three features
					  UWCW_FLOAT *f2,
					  UWCW_FLOAT *f3)
{
	UWCW_INT numSeg, locLen, i, distLen;
	UWCW_INT *loc;
	UWCW_FLOAT *dist, *center;

	loc = (UWCW_INT*)calloc(size+1, sizeof(UWCW_INT));
	locLen = compareThreshold(data, size, T, &numSeg, loc);

	for (i=0; i<locLen-1; i++)
	{
		loc[i]=loc[i+1];
	}
	locLen--;

	if (locLen%2)
	{
		loc[locLen]=size;
		locLen++;
	}

	center = (UWCW_FLOAT*)calloc(locLen/2, sizeof(UWCW_FLOAT));
	for (i=0; i<locLen/2; i++)
	{
		center[i] = loc[2*i+1]+loc[2*i];
		center[i] = center[i]/2.0;
	}

	if (locLen/2 <= 1)
	{
		dist = (UWCW_FLOAT*)calloc(1, sizeof(UWCW_FLOAT));
		dist[0] = 0;
		distLen = 1;
	}
	else
	{
		dist = (UWCW_FLOAT*)calloc(locLen/2-1, sizeof(UWCW_FLOAT));
		distLen = locLen/2-1;
		for (i=0; i<distLen; i++)
		{
			dist[i] = center[i+1]-center[i];
		}
	}


	*f1 = (UWCW_FLOAT)numSeg/(UWCW_FLOAT)size*100.0;
	*f2 = sumArray(dist, 0, distLen) / (UWCW_FLOAT)distLen;
	*f3 = uwStd(dist, distLen);

	if (loc) free(loc);
	if (center) free(center);
	if (dist) free(dist);
}

void computeFeature(UWCW_FLOAT *data,		// input data
					UWCW_INT dataSize,		// data size
					UWCW_INT specSize,		// DFT size for computing power spectrum
					UWCW_INT autoSize,		// DFT size for computing autocorrelation
					UWCW_FLOAT *feature,	// output feature
					UWCW_INT featSize)
{
	UWCW_INT i;

	UWCW_FLOAT *specData, *autoData, *sortData;
	UWCW_FLOAT *spec, *autoCo;
	UWCW_FLOAT iD;

	specData= (UWCW_FLOAT*)calloc(specSize, sizeof(UWCW_FLOAT));
	autoData= (UWCW_FLOAT*)calloc(autoSize, sizeof(UWCW_FLOAT));
	spec	= (UWCW_FLOAT*)calloc(specSize, sizeof(UWCW_FLOAT));
	autoCo	= (UWCW_FLOAT*)calloc(autoSize, sizeof(UWCW_FLOAT));

	memcpy(specData, data, __min(dataSize, specSize)*sizeof(UWCW_FLOAT));
	memcpy(autoData, data, __min(dataSize, autoSize)*sizeof(UWCW_FLOAT));

	// feature 1&2: power spectrum distribution
	powerSpectrum(specData, specSize, spec);
	if (spec[0]!=0)
	{
		feature[0] = sumArray(spec, 1, 5)/spec[0]*100.0;
		feature[1] = sumArray(spec, 5, specSize)/spec[0]*100.0;
	}

	// feature 3~6: statistics
	sortData = (UWCW_FLOAT*)calloc(dataSize, sizeof(UWCW_FLOAT));
	memcpy(sortData, data, dataSize*sizeof(UWCW_FLOAT));

	qsort(sortData, dataSize, sizeof(UWCW_FLOAT), dcomp);
	
	feature[2] = sumArray(data, 0, dataSize)/(UWCW_FLOAT)dataSize;
	feature[3] = sortData[dataSize-1];
	iD = (UWCW_FLOAT)dataSize*0.95;
	i = (UWCW_INT)iD;
	if ((UWCW_FLOAT)i < iD) i++;
	feature[4] = sortData[i-1];
	iD = (UWCW_FLOAT)dataSize*0.68;
	i = (UWCW_INT)iD;
	if ((UWCW_FLOAT)i < iD) i++;
	feature[5] = sortData[i-1];

	// feature 7~9: threshold with 95%
	thresholdFeature(data, dataSize, feature[4], &feature[6], &feature[7], &feature[8]);
	
	// feature 10~12: threshold with 68%
	thresholdFeature(data, dataSize, feature[5], &feature[9], &feature[10], &feature[11]);

	// feature 13~14: autocorrelation
	for (i=0; i<__min(dataSize,autoSize); i++)
	{
		autoData[i] += 72.0;
	}
	autoCorrelation(autoData, autoSize, autoCo);
	feature[13] = sumArray(autoCo, 0, autoSize)/(UWCW_FLOAT)autoSize;	// mean
	qsort(autoCo, autoSize, sizeof(UWCW_FLOAT), dcomp);	
	
	// median
	if (autoSize % 2)		// exact median
		feature[12] = autoCo[(UWCW_INT)(autoSize/2)];
	else
		feature[12] = (autoCo[autoSize/2] + autoCo[autoSize/2-1])/2.0;


	// clean up
	if (specData) free(specData);
	if (autoData) free(autoData);
	if (spec) free(spec);
	if (autoCo) free (autoCo);
	if (sortData) free(sortData);
}

void blockThreshold(UWCW_FLOAT *Bscope,
					UWCW_INT M,
					UWCW_INT N,
					UWCW_INT *BWBscope,
					UWCW_INT sliceHeight,
					UWCW_INT sliceWidth,
					UWCW_FLOAT percentage,
					UWCW_FLOAT threshGate,
					UWCW_INT y1,
					UWCW_INT y2)
{
	UWCW_INT nVer, nHor;
	UWCW_FLOAT nVerD, nHorD, T, T2;
	UWCW_INT i, j, ii, jj, h1, h2, w1, w2;
	UWCW_FLOAT p2 = LOWBOUND_PCT;

	nVerD = (UWCW_FLOAT)M / (UWCW_FLOAT)sliceHeight;
	nVer  = (UWCW_INT)nVerD;
	if ((UWCW_FLOAT)nVer < nVerD) nVer++;

	nHorD = (UWCW_FLOAT)N / (UWCW_FLOAT)sliceWidth;
	nHor  = (UWCW_INT)nHorD;
	if ((UWCW_FLOAT)nHor < nHorD) nHor++;

	for (i=0; i<nVer; i++)
	{
		for (j=0; j<nHor; j++)
		{
			h1 = i*sliceHeight;
			h2 = (i+1)*sliceHeight;
			w1 = j*sliceWidth;
			w2 = (j+1)*sliceWidth;
			if (h2>M) h2 = M;
			if (w2>N) w2 = N;
			if (y1>=h2 || y2<=h1) continue;

			T = blockPercent(Bscope, M, N, h1, h2, w1, w2, percentage, &T2, &p2);
			if (T<threshGate) T=threshGate;
			if (LOWBOUND_PCT<100 && T2<threshGate) T2 = threshGate;

			for (ii=h1; ii<h2; ii++)
			{
				for (jj=w1; jj<w2; jj++)
				{
					if (LOWBOUND_PCT>=100)
						BWBscope[ii*N+jj] = (Bscope[ii*N+jj]>=T);
					else
						BWBscope[ii*N+jj] = (UWCW_INT)(Bscope[ii*N+jj]-T2+0.5)*(Bscope[ii*N+jj]>=T);
				}
			}
		}
	}
}

UWCW_FLOAT blockPercent(UWCW_FLOAT *Bscope,
					UWCW_INT M,  UWCW_INT N,
					UWCW_INT h1, UWCW_INT h2,
					UWCW_INT w1, UWCW_INT w2,
					UWCW_FLOAT percentage,
					UWCW_FLOAT *T2,
					UWCW_FLOAT *p2)
{
	UWCW_FLOAT T, *blockData, *pt;
	UWCW_INT i, j, M2, N2;	

	M2 = (UWCW_INT)(UWCW_FLOAT(h2-h1)/2+0.5);
	N2 = (UWCW_INT)(UWCW_FLOAT(w2-w1)/2+0.5);
	blockData = (UWCW_FLOAT*)calloc(sizeof(UWCW_FLOAT), M2*N2);

	pt = blockData;
	for (i=h1; i<h2; i+=2)
	{
		for (j=w1; j<w2; j+=2)
		{
			*pt = Bscope[i*N+j];
			pt++;
		}
	}

	qsort(blockData, M2*N2, sizeof(UWCW_FLOAT), dcomp);
	i = (UWCW_INT)((UWCW_FLOAT)M2*(UWCW_FLOAT)N2*percentage/100.0);

	T = blockData[i];

	if (p2 && T2 && *p2 < 100)
	{
		i = (UWCW_INT)((UWCW_FLOAT)M2*(UWCW_FLOAT)N2*(*p2)/100.0);
		*T2 = blockData[i];
	}

	free(blockData);
	return T;
}

int dcomp(const void *a, const void *b)
{
	if (*(UWCW_FLOAT*)a < *(UWCW_FLOAT*)b)
		return -1;
	else if (*(UWCW_FLOAT*)a > *(UWCW_FLOAT*)b)
		return 1;
	else
		return 0;
}

void coorTransform(UWCW_INT *BWBscope,		// input b-n-w Bscope image
				   UWCW_INT M,				// # of rows of BWBscope
				   UWCW_INT N,				// # of columns of BWBscope
				   UWCW_INT *BWcart,			// output b-n-w cartesian image, shall be allocated outside
				   UWCW_INT *pNC,			// # of columns of BWcart
				   UWCW_INT *px0,			// middle point in BWcart
				   UWCW_FLOAT Ang_endpt,	// parameter, angle endpoint, in degrees
				   UWCW_FLOAT *Ang_plot,	// 1 by N vector angle value, in degrees
				   UWCW_INT interFac,		// interpolation factor
				   UWCW_INT y1,				// starting and ending row
				   UWCW_INT y2)
{
	UWCW_FLOAT x0D, t2, t;
	UWCW_INT i, x0, NC, m, n;

	// cartesian image size
	x0D = M*sin(Ang_endpt*UW_PI/180);
	x0 = (UWCW_INT)(x0D)+1;
	NC = 2*x0+1;

	// loop through input Bscope
	for (m=y1; m<y2-1; m++)
	{
		for (n=0; n<N; n++)
		{
			if (BWBscope[m*N+n])
			{
				UWCW_INT val = BWBscope[m*N+n];
				coorTransAssign(BWcart, M, NC, m, Ang_plot[n]*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);

				// angle interpolation
				if (interFac >= 1)
				{
					// interpolate before
					if (n>0)
					{
						t2 = (Ang_plot[n-1]+Ang_plot[n])/2;	// midpoint before
						for (i=0; i<interFac; i++)
						{
							t = ((interFac-i)*t2+i*Ang_plot[n])/(UWCW_FLOAT)(interFac);
							coorTransAssign(BWcart, M, NC, m, t*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);
						}
					}

					// interpolate after
					if (n<N-1)
					{
						t2 = (Ang_plot[n]+Ang_plot[n+1])/2;	// midpoint after
						for (i=1; i<=interFac; i++)
						{
							t = ((interFac-i)*Ang_plot[n]+i*t2)/(UWCW_FLOAT)(interFac);
							coorTransAssign(BWcart, M, NC, m, t*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);
						}
					}
				}
			}
		}
	}
}

void coorTransAssign(UWCW_INT *BWcart,		// b-n-w cartesian image
					 UWCW_INT MC,			// # of rows of BWcart
					 UWCW_INT NC,			// # of columns of BWcart
					 UWCW_FLOAT rho,		// rho in Bscope
					 UWCW_FLOAT theta,		// theta in Bscope
					 UWCW_INT x0,
					 UWCW_INT val)
{
	UWCW_FLOAT xD, yD;
	UWCW_INT x, y;

	xD = (rho+1)*sin(theta)+x0;
	yD = (rho+1)*cos(theta);

	x = (UWCW_INT)(xD-0.5);
	y = (UWCW_INT)(yD-0.5);

	if (y>=0 && y<MC)
	{
		if (x>=0 && x<NC)
		{
			if (val)
				BWcart[y*NC+x] = __max(val, BWcart[y*NC+x]);
			else
				BWcart[y*NC+x] = 1;
		}
	}
}

void houghTransform(UWCW_INT *BWcart,		// input b-n-w cartesian image
					UWCW_INT MC,				// # of rows of BWcart
					UWCW_INT NC,				// # of columns of BWcart
					UWCW_INT y1,				// starting and ending processing region
					UWCW_INT y2,
					UWCW_FLOAT *rho,		// rho
					UWCW_INT rhoLen,
					UWCW_FLOAT *theta,		// theta
					UWCW_INT thetaLen,
					UWCW_INT *H)			// output image, sizeof (rhoLen rows, thetaLen columns)
{
	UWCW_INT i, thetaIdx, rhoIdx, x, y;	
	UWCW_FLOAT *cost = NULL, *sint = NULL, firstRho, slope, oneRho;	

	// cos/sin lookup table
	cost = (UWCW_FLOAT *)calloc(thetaLen, sizeof(UWCW_FLOAT));
	sint = (UWCW_FLOAT *)calloc(thetaLen, sizeof(UWCW_FLOAT));
	for (i=0; i<thetaLen; i++)
	{
		cost[i] = cos(theta[i]*UW_PI/180);
		sint[i] = sin(theta[i]*UW_PI/180);
	}

	// rho conversion from value to index
	firstRho = rho[0];
	slope = ((UWCW_FLOAT)rhoLen-1) / (rho[rhoLen-1]-firstRho);


	// Hough Transform
	for (y=y1; y<y2; y++)
	{
		for (x=0; x<NC; x++)
		{
			if (BWcart[y*NC+x])
			{
				for (thetaIdx=0; thetaIdx<thetaLen; thetaIdx++)
				{
					oneRho = (UWCW_FLOAT)x*cost[thetaIdx] + (UWCW_FLOAT)(y-y1)*sint[thetaIdx];
					rhoIdx = (UWCW_INT)(slope*(oneRho-firstRho)+0.5);
					H[rhoIdx*thetaLen+thetaIdx] += BWcart[y*NC+x];
				}
			}
		}
	}

	if (cost) free(cost);
	if (sint) free(sint);
}

UWCW_FLOAT entropyMax3(UWCW_INT *h, UWCW_INT rhoLen, UWCW_INT ss)
{
	UWCW_FLOAT w = 0;
	UWCW_INT j, k, l, side(50);
	UWCW_INT localMaxSum, localMax, maxIdx[3];
	bool ignoreThis(false);
	
	// inverse entropy
	for (j=0; j<rhoLen; j++)
	{
		if (h[j*ss])
			w += h[j*ss] * log((UWCW_FLOAT)h[j*ss]);
	}

	// sum of three local maxima
	localMaxSum = 0;
	for (k=0; k<3; k++)
	{
		localMax = -1;
		maxIdx[k] = -1;
		for (j=0; j<rhoLen; j++)
		{
			ignoreThis = false;
			for (l=0; l<k; l++)
			{
				if (maxIdx[l]-j >= -side && maxIdx[l]-j <= side)
				{
					ignoreThis = true;
				}
			}

			if (!ignoreThis && h[j*ss]>localMax)
			{
				localMax = h[j*ss];
				maxIdx[k] = j;
			}
		}
		localMaxSum += localMax;
	}	

	w *= (UWCW_FLOAT)localMaxSum;	
	return w;
}

tparticle* theta_initialize_particles(UWCW_INT num, UWCW_FLOAT theta, UWCW_INT* H, UWCW_INT rhoLen, UWCW_INT thetaLen, UWCW_INT idx)
{
	tparticle* theta_particles = (tparticle*)calloc(num, sizeof(tparticle));
	
	if (theta<0) theta = theta+180;
	for (UWCW_INT i=0; i<num; i++)
	{
		theta_particles[i].theta = theta;
		theta_particles[i].thetap = 0;
		theta_particles[i].thetapp = 0;
		theta_particles[i].w = 1;
	}

	return theta_particles;
}

void theta_free_particles(tparticle* particles, UWCW_INT num)
{
	free(particles);
}

UWCW_FLOAT theta_average(tparticle* theta_particles, UWCW_INT num)
{
	UWCW_FLOAT theta_return = 0;
	UWCW_INT i;
	for (i=0; i<num; i++)
	{
		theta_return += theta_particles[i].w * theta_particles[i].theta;
	}

	return theta_return;
}

void theta_transition_particles(tparticle* particles, UWCW_INT num, UWCW_FLOAT theta_std)
{
	UWCW_FLOAT tn, t, tp, tpp, dum1, dum2;
	UWCW_INT i;

	for (i=0; i<num; i++)
	{
		t	= particles[i].theta;
		tp	= particles[i].thetap;
		tpp	= particles[i].thetapp;
		if (tp>0 && tpp>0)
			tn = 3*t - 3*tp + tpp;
		else if (tp>0 && tpp<=0)
			tn = 2*t - tp;
		else
			tn = t;
		gauss(theta_std, theta_std, &dum1, &dum2);
		tn += dum1;
		//tn += (i - num/2);

		if (tn<PTHETA_END) tn = (UWCW_FLOAT)PTHETA_END;
		if (tn>180-PTHETA_END) tn = (UWCW_FLOAT)(180-PTHETA_END);

		particles[i].thetapp = tp;
		particles[i].thetap  = t;
		particles[i].theta   = tn;
	}
}

void gauss(UWCW_FLOAT std1, UWCW_FLOAT std2, UWCW_FLOAT *y1, UWCW_FLOAT *y2)
{
	UWCW_FLOAT x1, x2, w;

	do
	{
		x1 = 2.0f * rand() / (UWCW_FLOAT)(RAND_MAX) - 1.0;
		x2 = 2.0f * rand() / (UWCW_FLOAT)(RAND_MAX) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);

	w = sqrt((-2.0f * log(w)) / w);
	*y1 = x1 * w * std1;
	*y2 = x2 * w * std2;
}

void theta_compute_likelihood(tparticle* theta_particles, UWCW_INT num, UWCW_INT* H,
	UWCW_INT rhoLen, UWCW_INT thetaLen, UWCW_FLOAT(*similarity)(UWCW_INT *h1, UWCW_INT *h2, UWCW_INT n))
{
	UWCW_INT i, j, *h, idx;
	UWCW_FLOAT w;

	h = (UWCW_INT*)calloc(rhoLen, sizeof(UWCW_INT));
	for (i=0; i<num; i++)
	{
		idx = thetaToIdx(theta_particles[i].theta, PTHETA_END, PTHETA_LEN);

		for (j=0; j<rhoLen; j++)
		{
			h[j] = H[j*thetaLen+idx];
		}

		w = 1;//(*similarity)(theta_particles[i].h, h, rhoLen);
		w *= entropyMax3(h, rhoLen, 1);
		w *= exp(-2*(theta_particles[i].theta-theta_previous)*(theta_particles[i].theta-theta_previous)/(100*THETA_VARIANCE*THETA_VARIANCE));
		theta_particles[i].w = w;
	}

	if (h) free(h);
}

UWCW_INT thetaToIdx(UWCW_FLOAT theta, UWCW_INT thetaEnd, UWCW_INT thetaLen)
{
	UWCW_INT thetaI, idx;

	thetaI = (UWCW_INT)(theta+0.5);
	if (thetaI < thetaEnd) thetaI = thetaEnd;
	if (thetaI > 180-thetaEnd) thetaI = 180 - thetaEnd;
	
	if (thetaI >= 90) idx = thetaI - 90;
	else idx = thetaI - (90 - thetaLen);

	return idx;
}

void theta_normalize_weights(tparticle *particles, UWCW_INT num)
{
	UWCW_FLOAT wsum = 0;
	UWCW_INT i;

	for (i=0; i<num; i++)
	{
		wsum += particles[i].w;
	}

	for (i=0; i<num; i++)
	{
		particles[i].w /= wsum;
	}
}

tparticle* tresample(tparticle* ptcs, UWCW_INT n)
{
	tparticle *ptcsn;
	UWCW_INT i, j, k=0, np;

	qsort(ptcs, n, sizeof(tparticle), &tparticleCmp);
	ptcsn = (tparticle*)calloc(n, sizeof(tparticle));
	for (i=0; i<n; i++)
	{
		np = (UWCW_INT)(ptcs[i].w * n + 0.5);
		for (j=0; j<np; j++)
		{
			ptcsn[k] = ptcs[i];
			k++;
			if (k==n)
				return ptcsn;
		}
	}

	while (k<n)
	{
		ptcsn[k] = ptcs[0];
		k++;
	}
	return ptcsn;
}

int tparticleCmp(const void *p1, const void *p2)
{
	tparticle *_p1 = (tparticle*)p1;
	tparticle *_p2 = (tparticle*)p2;

	if (_p1->w > _p2->w)
		return -1;
	if (_p1->w < _p2->w)
		return 1;
	return 0;
}

void adjustParameters(UWCW_INT M, UWCW_INT N)
{
	if (M>4000)
		hFactor = 4;
	else if (M>2000)
		hFactor = 2;
	else
		hFactor = 1;
	
	if (N>300)
		wFactor = 1;
	else
		wFactor = 2;

	SLICE_HEIGHT = M / 32;
	SLICE_WIDTH = 64/wFactor;
	if (M>4000) HGATE_FACTOR = 24;
}

UWCW_INT numCable(cableLine *cl, UWCW_INT num)
{
	int n(0), i;
	for (i=0; i<num; i++)
	{
		if (cl[i].isCable && cl[i].score>=0.5)
		{
			n++;
		}
	}

	return n;
}

UWCW_INT checkParallel(cableLine *cl, UWCW_INT num)
{
	UWCW_INT i, j;
	UWCW_INT *nArray = (UWCW_INT*)calloc(num, sizeof(UWCW_INT));
	for (i=0; i<num; i++)
	{
		if (cl[i].isCable && cl[i].score>=.5)
		for (j=0; j<num; j++)
		{
			if (!cl[j].isCable) continue;
			if (cl[j].score<.5) continue;
			if (fabs(cl[j].theta-cl[i].theta)<=.5)
				nArray[i]++;
		}
	}

	j = 0;
	for (i=0; i<num; i++)
	{
		if (j<nArray[i])
			j = nArray[i];
	}

	free(nArray);
	return (j-1);
}

UWCW_FLOAT lineCorrelate(cableLine l1, cableLine l2, UWCW_FLOAT Rsigma, UWCW_FLOAT Tsigma, UWCW_FLOAT C)
{
	UWCW_FLOAT cr, ca, score;
	cr = exp(-(fabs(l1.rho)-fabs(l2.rho))*(fabs(l1.rho)-fabs(l2.rho))/(2*Rsigma*Rsigma));
	ca = exp(-(fabs(l1.theta)-fabs(l2.theta))*(fabs(l1.theta)-fabs(l2.theta))/(2*Tsigma*Tsigma));
	score = cr*ca*C;

	return score;
}

UWCW_INT temporalCorrelate(cableLine *current, UWCW_INT currentNum,
	cableLine *prev, UWCW_INT prevNum)
{
	static UWCW_FLOAT RSIGMA = 150.0;
	static UWCW_FLOAT TSIGMA = 3.0;
	static UWCW_FLOAT C = exp(1.0);
	UWCW_INT n(0), i, j;
	UWCW_FLOAT *corre = NULL;

	if (!prev || !prevNum)
	{
		return 0;
	}

	corre = (UWCW_FLOAT*)calloc(prevNum, sizeof(UWCW_FLOAT));

	for (i=0; i<currentNum; i++)
	{
		if (current[i].isCable && current[i].score>=.5)
		{
			memset(corre, 0, sizeof(UWCW_FLOAT)*prevNum);
			for (j=0; j<prevNum; j++)
			{
				corre[j] = lineCorrelate(current[i], prev[j], RSIGMA, TSIGMA, C);
			}

			for (j=1; j<prevNum; j++)
			{
				if (corre[j]>corre[0])
					corre[0]=corre[j];
			}
			if (corre[0]>=1) n++;
			current[i].score = current[i].score*corre[0];
			current[i].isCable = current[i].isCable && (current[i].score>=.5);
		}
	}

	free(corre);
	return n;
}

UWCW_FLOAT trackTheta( UWCW_INT *BWcart,		// input b-n-w cartesian image
					UWCW_INT MC,				// size of BWcart,
					UWCW_INT NC,
				 	UWCW_FLOAT *theta,			// HT parameter
					UWCW_INT thetaLen,
					UWCW_INT **Hout, UWCW_FLOAT **rhoOut, UWCW_INT *rhoLenOut)
{
	UWCW_INT q, rhoLen, M, N, i, maxIdx, *H;
	UWCW_FLOAT d, *rho, maxSum, *thetaSum, thetaToReturn(0);

	M = MC;
	N = NC;
	d = sqrt((UWCW_FLOAT)(M-1)*(M-1)+(UWCW_FLOAT)(N-1)*(N-1));
	q = (UWCW_INT)(d)+1;
	rhoLen = 2*q+1;
	rho = (UWCW_FLOAT*)calloc(rhoLen, sizeof(UWCW_FLOAT));
	H = (UWCW_INT*)calloc(rhoLen*thetaLen, sizeof(UWCW_INT));
	thetaSum = (UWCW_FLOAT*)calloc(thetaLen, sizeof(UWCW_FLOAT));

	memset(H, 0, sizeof(UWCW_INT)*rhoLen*thetaLen);
	memset(thetaSum, 0, sizeof(UWCW_INT)*thetaLen);
	for (i=-q; i<=q; i++)
	{
		rho[i+q] = (UWCW_FLOAT)i;
	}

	houghTransform(BWcart, MC, NC, 0, MC, rho, rhoLen, theta, thetaLen, H);	
	
	if (theta_first_frame)
	{	
		// get theta sum
		maxSum = -1;
		maxIdx = -1;
		for (i=0; i<thetaLen; i++)
		{	
			thetaSum[i] = entropyMax3(H+i, rhoLen, thetaLen);
			if (maxSum < thetaSum[i])
			{
				maxSum = thetaSum[i];
				maxIdx = i;
			}
		}
		thetaToReturn = theta[maxIdx];
		if (thetaToReturn<0) thetaToReturn += 180;
		theta_particles = theta_initialize_particles(THETA_NUM_PARTICLES, thetaToReturn, H, rhoLen, thetaLen, maxIdx);
	}
	else
	{
		theta_transition_particles(theta_particles, THETA_NUM_PARTICLES, THETA_VARIANCE);
		theta_compute_likelihood(theta_particles, THETA_NUM_PARTICLES, H, rhoLen, thetaLen, &similarity_1);
		theta_normalize_weights(theta_particles, THETA_NUM_PARTICLES);
		tparticle *theta_particles_n = tresample(theta_particles, THETA_NUM_PARTICLES);
		
		theta_free_particles(theta_particles, THETA_NUM_PARTICLES);
		theta_particles = theta_particles_n;
		thetaToReturn = theta_particles[0].theta;
		thetaToReturn = (UWCW_INT)(thetaToReturn + 0.5);
	}	
	
	theta_previous = thetaToReturn;

	//if (rho) free(rho);
	//if (H) free(H);
	*rhoOut = rho;
	*rhoLenOut = rhoLen;
	*Hout = H;
	if (thetaSum) free(thetaSum);

	return thetaToReturn;
}

UWCW_FLOAT similarity_1(UWCW_INT *h1, UWCW_INT *h2, UWCW_INT n)
{
	return (UWCW_FLOAT)1;
}

UWCW_FLOAT get_rho_thresh(UWCW_INT *h, UWCW_INT rhoLen)
{
	UWCW_FLOAT vMax, thresh;
	UWCW_INT iMax, i, iDum;

	vMax = -1;
	iMax = -1;
	thresh = 0;
	for (i=0; i<rhoLen; i++)
	{
		thresh += h[i];
		if (h[i]>vMax)
		{
			vMax = h[i];
			iMax = i;
		}
	}


	UWCW_INT *h2 = (UWCW_INT*)calloc(rhoLen, sizeof(UWCW_INT));
	memcpy(h2, h, rhoLen*sizeof(UWCW_INT));

	for (UWCW_INT k=0; k<0; k++)
	{
		for (i=iMax-15; i<=iMax+15; i++)
		{
			if (i>=0 && i<rhoLen)
				h2[i] = 0;
		}
		
		getMaxVal(h2, rhoLen, 1, &vMax, &iMax, &iDum);
	}
	
	thresh /= rhoLen;
	thresh += vMax / 3;

	free(h2);
	return thresh;
}

void getMaxVal(	UWCW_INT *H, UWCW_INT M, UWCW_INT N,
	UWCW_FLOAT *maxVal, UWCW_INT *maxRow, UWCW_INT *maxCol)
{
	UWCW_INT i, j;
	*maxVal = -1;
	*maxRow = -1;
	*maxCol = -1;

	UWCW_INT val;
	
	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			val = H[i*N+j];
			if (*maxVal < val)
			{
				*maxVal = val;
				*maxRow = i;
				*maxCol = j;
			}
		}
	}
}

void rho_transition_particles(rparticle *rho_particles, UWCW_INT num, UWCW_FLOAT rho_std, UWCW_FLOAT theta_new, UWCW_INT ss, rtracker *r_track)
{
	UWCW_FLOAT rn, r, rp, rpp, dum1, dum2;
	UWCW_INT i;
	bool turn_over = false;

	// condition of turning over 90 degree line
	if (r_track->theta >= 90 && theta_new < 90
		|| r_track->theta < 90 && theta_new >= 90)
		turn_over = true;

	if (r_track)
	{
		r	= r_track->rho;
		rp	= r_track->rhop;
		rpp = r_track->rhopp;

		if (turn_over)
			rn = -r;
		else
			rn = r;
		
		for (i=0; i<num; i++)
		{
			rho_particles[i].rhopp = rp;
			rho_particles[i].rhop  = r;
			rho_particles[i].rho   = rn + (i-num/2)*ss;

			rho_particles[i].thetapp= rho_particles[i].thetap;
			rho_particles[i].thetap	= rho_particles[i].theta;
			rho_particles[i].theta	= theta_new;
		}
	}
	else
	{
		for (i=0; i<num; i++)
		{
			r	= rho_particles[i].rho;
			rp	= rho_particles[i].rhop;
			rpp	= rho_particles[i].rhopp;

			if (turn_over)
				rn = -r;
			else
				rn = r;

			gauss(rho_std, rho_std, &dum1, &dum2);
			rn += dum1;

			rho_particles[i].rhopp	= rp;
			rho_particles[i].rhop	= r;
			rho_particles[i].rho	= (int)(rn+0.5);

			rho_particles[i].thetapp= rho_particles[i].thetap;
			rho_particles[i].thetap	= rho_particles[i].theta;
			rho_particles[i].theta	= theta_new;
		}
	}
}

UWCW_FLOAT rho_compute_likelihood(rparticle *rho_particles, UWCW_INT num, UWCW_FLOAT rhop,
	UWCW_FLOAT theta, UWCW_INT *h, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0, svmData *svm, rtracker *r_track)
{
	UWCW_INT i, ridx, ridxp, j, dataNum, lineRange[512], lineAngle[512];
	UWCW_FLOAT linearVal, svmVal, w, lineData[512], f[FEAT_SIZE], base = 0;
	UWCW_FLOAT data1[512], data2[512];
	UWCW_FLOAT asso_s = -1, rhod;

	for (j=0; j<2*HOUGH_SIDE+1; j++)
	{
		base += r_track->hough_[j] * r_track->hough_[j];
	}

	for (i=0; i<num; i++)
	{
		ridx = rho_to_idx(rho_particles[i].rho, rho, rhoLen);
		w = 1;
		
		// current SVM
		dataNum = getLineData(RP_dBm, Ang_plot, M, N,
			BW, x0, rho_particles[i].rho, rho_particles[i].theta<90 ? rho_particles[i].theta : rho_particles[i].theta-180,
			lineData, lineRange, lineAngle, false);
		computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
		linearVal = linearCalc(svm, f);
		svmVal = svmCalc(svm, f);
		if (dataNum<15 || linearVal<0) linearVal = 0;
		if (svmVal > SVM_BASE) svmVal = SVM_BASE;
		
		w = 0;
		for (j=-HOUGH_SIDE; j<=HOUGH_SIDE; j++)
		{
			if (ridx+j>=0 && ridx+j<rhoLen)
			{
				w += r_track->hough_[j+HOUGH_SIDE] * h[ridx+j];
			}
		}

		w /= base;
		rhod = fabs(rho_particles[i].rho) - fabs(rhop);
		w *= exp(-2*rhod*rhod/(400*RHO_VARIANCE*RHO_VARIANCE));
		w *= linearVal * (-svmVal + SVM_BASE);

		if (w > asso_s)
			asso_s = w;

		// update
		rho_particles[i].w = w;
		rho_particles[i].hVal = h[ridx];		
	}

	return asso_s;
}

UWCW_INT rho_to_idx(UWCW_FLOAT rho, UWCW_FLOAT *rhos, UWCW_INT rhoLen)
{
	UWCW_FLOAT idxD;
	UWCW_INT idx;
	idxD = (rho-rhos[0])/(rhos[1]-rhos[0]);
	idx = (UWCW_INT)(idxD+0.5);
	return idx;
}

void rho_normalize_weights(rparticle *rpts, UWCW_INT num)
{
	UWCW_FLOAT sum = 0;
	UWCW_INT i;
	for (i=0; i<num; i++)
		sum += rpts[i].w;
	for (i=0; i<num; i++)
		rpts[i].w /= sum;
}

int rparticleCmp(const void *p1, const void *p2)
{
	rparticle *_p1 = (rparticle*)p1;
	rparticle *_p2 = (rparticle*)p2;

	if (_p1->w > _p2->w)
		return -1;
	if (_p1->w < _p2->w)
		return 1;
	return 0;
}

rparticle* rresample(rparticle* ptcs, UWCW_INT n)
{
	rparticle *ptcsn;
	UWCW_INT i, j, k=0, np;

	qsort(ptcs, n, sizeof(rparticle), &rparticleCmp);
	ptcsn = (rparticle*)calloc(n, sizeof(rparticle));
	for (i=0; i<n; i++)
	{
		np = (UWCW_INT)(ptcs[i].w * n + 0.5);
		for (j=0; j<np; j++)
		{
			ptcsn[k] = ptcs[i];
			k++;
			if (k==n)
				return ptcsn;
		}
	}

	while (k<n)
	{
		ptcsn[k] = ptcs[0];
		k++;
	}
	return ptcsn;
}

void rtracker_free(rtracker *rtracker_)
{
	if (rtracker_->data_original_) free(rtracker_->data_original_);
	if (rtracker_->data_threshold_) free(rtracker_->data_threshold_);
	if (rtracker_->rpts_) free(rtracker_->rpts_);
	if (rtracker_->hough_) free(rtracker_->hough_);
	
	rtracker_->data_original_ = NULL;
	rtracker_->data_threshold_ = NULL;
	rtracker_->rpts_ = NULL;
	rtracker_->hough_ = NULL;
	rtracker_->rpts_num = 0;
	rtracker_->hough_num = 0;
	rtracker_->rho = 0;
	rtracker_->rhop = 0;
	rtracker_->rhopp = 0;
	rtracker_->theta = 0;
	rtracker_->thetap = 0;
	rtracker_->thetapp = 0;
}

void process_rtrackers(rtracker *rtracker_,
	UWCW_FLOAT theta, UWCW_INT *h, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0, svmData *svm, UWCW_INT margin)
{
	rtracker *current = rtracker_;
	rtracker *previous;
	rparticle *rptsn;
	UWCW_FLOAT asso_s;	// associatio score
	UWCW_INT dataNum, i, j, ridx;
	UWCW_FLOAT lineData[512];
	UWCW_INT lineRange[512];
	UWCW_INT lineAngle[512];
	UWCW_INT hidx;
	UWCW_FLOAT f[FEAT_SIZE];
	bool del = false;

	UWCW_INT ss = 1;
	ss = rhoLen>6000 ? 2 : 1;
	if (rhoLen>6000) RHO_VARIANCE = 4;

	while (current)
	{
		del = false;
		if (current->rstatus == candidate || current->rstatus == healthy)
		{
			// run the particle filtering algorithm
			rho_transition_particles(current->rpts_, current->rpts_num, RHO_VARIANCE, theta, ss, current);
			asso_s = rho_compute_likelihood(current->rpts_, current->rpts_num, current->rho,
				theta, h, rho, rhoLen, thetaLen,
				RP_dBm, BW, BWcart, Ang_plot,
				M, N, NC, x0, svm, current);
			rho_normalize_weights(current->rpts_, current->rpts_num);
			rptsn = rresample(current->rpts_, current->rpts_num);

			free(current->rpts_);
			current->rpts_ = rptsn;

			// update tracker
			current->rhopp = current->rhop;
			current->rhop  = current->rho;
			current->rho   = current->rpts_[0].rho;
			
			current->thetapp = current->thetap;
			current->thetap  = current->theta;
			current->theta   = theta;

			// update hough datas
			hidx = rho_to_idx(current->rho, rho, rhoLen);
			memcpy(current->hough_, h+hidx-HOUGH_SIDE, (2*HOUGH_SIDE+1)*sizeof(UWCW_INT));
			
			dataNum = getLineData(RP_dBm, Ang_plot, M, N,
				0, x0, current->rho, current->theta<90 ? current->theta : current->theta-180,
				lineData, lineRange, lineAngle, false);
			//computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
			current->data_num = dataNum;

			free(current->data_original_);
			free(current->data_threshold_);
			current->data_original_ = (UWCW_FLOAT*)calloc(dataNum, sizeof(UWCW_FLOAT));
			current->data_threshold_ = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_FLOAT));
			memcpy(current->data_original_, lineData, sizeof(UWCW_FLOAT)*dataNum);
			for (i=0; i<dataNum; i++)
			{
				current->data_threshold_[i] = BW[ lineRange[i]*N+lineAngle[i] ];
			}
			
			// output results
			if (asso_s >= ASSOCIATION_THRESHOLD)
			{
				current->rstatus = healthy;
				current->life_time++;
				current->status_time++;
				// suppress h
				ridx = rho_to_idx(current->rho, rho, rhoLen);
				for (j=ridx-margin; j<=ridx+margin; j++)
				{
					if (j>=0 && j<rhoLen)
						h[j] = 0;
				}
			}
			else
			{
				
				if (current->rstatus == healthy)
				{
					current->rstatus = dormant;
					current->life_time++;
					current->status_time = 1;
				}
				else
				{
					// candidate with no association: kill it
					previous->next_ = current->next_;
					rtracker_free(current);
					del = true;
				}
			}
		}
		else if (current->rstatus == dormant)
		{
			if (current->status_time >= DORMANCY_THRESHOLD)
			{
				previous->next_ = current->next_;
				rtracker_free(current);
				del = true;
			}
		}

		if (!del) previous = current;
		current = current->next_;
	}
}

rtracker *rtracker_init_zero()
{
	rtracker *rtracker_ = (rtracker*)calloc(1, sizeof(rtracker));
	rtracker_->data_num = 0;
	rtracker_->data_original_ = NULL;
	rtracker_->data_threshold_ = NULL;
	rtracker_->life_time = 0;
	rtracker_->rho = 0;
	rtracker_->rhop = 0;
	rtracker_->rhopp = 0;
	rtracker_->rpts_num = 0;
	rtracker_->rpts_ = NULL;
	rtracker_->rstatus = head;
	rtracker_->status_time = 0;
	rtracker_->next_ = NULL;


	return rtracker_;
}

rparticle* rho_initialize_particles(UWCW_INT num, UWCW_FLOAT rho, UWCW_FLOAT theta, UWCW_FLOAT hVal, UWCW_FLOAT svmVal, UWCW_FLOAT linearVal)
{
	rparticle *rpts;
	UWCW_INT i;

	rpts = (rparticle*)calloc(num, sizeof(rparticle));
	for (i=0; i<num; i++)
	{
		rpts[i].rho		= rho;//>0 ? rho : -rho;
		rpts[i].rhop	= 0;
		rpts[i].rhopp	= 0;
		rpts[i].theta	= theta;
		rpts[i].thetap	= 0;
		rpts[i].thetapp	= 0;
		rpts[i].w		= 1;

		rpts[i].hVal = hVal;
		rpts[i].svmVal = svmVal;
		rpts[i].linearVal = linearVal;
	}

	return rpts;
}

rtracker *rtracker_init(UWCW_FLOAT rho, UWCW_FLOAT theta,
	UWCW_FLOAT *lineData, UWCW_INT *lineRange, UWCW_INT *lineAngle, UWCW_INT dataNum,
	UWCW_INT *BW, UWCW_INT M, UWCW_INT N, UWCW_FLOAT hVal, UWCW_FLOAT svmVal, UWCW_FLOAT linearVal,
	UWCW_INT *hough_data, UWCW_INT hnum)
{
	UWCW_INT i;
	rtracker *rtracker_ = (rtracker*)calloc(1, sizeof(rtracker));
	
	// initialize particles code
	rtracker_->rpts_ = rho_initialize_particles(RHO_NUM_PARTICLES, rho, theta, hVal, svmVal, linearVal);
	rtracker_->rpts_num = RHO_NUM_PARTICLES;
	
	// other stuff
	rtracker_->rstatus = candidate;
	rtracker_->life_time = 1;
	rtracker_->status_time = 1;

	rtracker_->rho = rho;//>0 ? rho : -rho;
	rtracker_->rhop = 0;
	rtracker_->rhopp = 0;
	rtracker_->theta = theta;
	rtracker_->thetap = 0;
	rtracker_->thetapp = 0;

	rtracker_->data_num = dataNum;
	rtracker_->data_original_ = (UWCW_FLOAT*)calloc(dataNum, sizeof(UWCW_FLOAT));
	rtracker_->data_threshold_ = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_INT));
	memcpy(rtracker_->data_original_, lineData, dataNum*sizeof(UWCW_FLOAT));
	for (i=0; i<dataNum; i++)
	{
		rtracker_->data_threshold_[i] = BW[ lineRange[i]*N+lineAngle[i] ];
	}

	// hough data storage
	rtracker_->hough_ = (UWCW_INT*)calloc(hnum, sizeof(UWCW_INT));
	memcpy(rtracker_->hough_, hough_data, hnum*sizeof(UWCW_INT));
	rtracker_->hough_num = hnum;

	rtracker_->next_ = NULL;

	return rtracker_;
}

UWCW_INT rtracker_free_list(rtracker *head)
{
	UWCW_INT num = 0;
	rtracker *current, *tmp;
	current = head;
	while (current)
	{
		rtracker_free(current);
		tmp = current->next_;
		free(current);
		current = tmp;
		num++;
	}
	return num;
}

UWCW_INT get_rtracker_num(rtracker *rtracker_)
{
	UWCW_INT num = 0;
	while (rtracker_)
	{
		if (fabs(rtracker_->rho)>0.01)
			num++;
		rtracker_ = rtracker_->next_;
	}
	return num;
}

UWCW_INT add_rtrackers(rtracker **rtrackers__, UWCW_INT *h, UWCW_INT rhoLen, UWCW_FLOAT theta,
	UWCW_FLOAT *rho, UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_FLOAT *Ang_plot,
	UWCW_INT M, UWCW_INT N, UWCW_INT x0, svmData *svm, UWCW_INT margin, UWCW_FLOAT thresh)
{
	rtracker *current;
	UWCW_FLOAT vMax, thetaN, linearVal, svmVal;
	UWCW_INT iMax, i, j, dataNum, num_added, numNow;
	UWCW_INT lineRange[512];
	UWCW_INT lineAngle[512];
	UWCW_FLOAT lineData[512];
	UWCW_FLOAT f[FEAT_SIZE];

	if (!(*rtrackers__))
	{
		*rtrackers__ = rtracker_init_zero();
	}
	
	// get to the place
	current = *rtrackers__;
	while (current->next_) current = current->next_;
	num_added=0;

	vMax = -1;
	iMax = -1;
	for (i=0; i<rhoLen; i++)
	{
		if (vMax < h[i])
		{
			vMax = h[i];
			iMax = i;
		}
	}

	thetaN = theta;
	if (thetaN>=90) thetaN -= 180;
	// get rho candidates:1. local maxima; 2. over threshold; 3. pass svm and linear
	numNow = get_rtracker_num(*rtrackers__);
	
	while (vMax>=thresh)
	{
		// constrain the total number of trackers
		numNow++;
		if (numNow > MAX_TRACKERS) break;

		dataNum = getLineData(RP_dBm, Ang_plot, M, N,
			BW, x0, rho[iMax], thetaN,
			lineData, lineRange, lineAngle, false);
		computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
		linearVal = linearCalc(svm, f);
		svmVal = svmCalc(svm, f);

		// save line Data here
		//saveLineData(lineData, dataNum);

		if (dataNum>=15 && (svmVal<0 && linearVal>=0.5))
		{
			num_added++;
			dataNum = getLineData(RP_dBm, Ang_plot, M, N,
				NULL, x0, rho[iMax], thetaN,
				lineData, lineRange, lineAngle, false);
			current->next_ = rtracker_init(rho[iMax], theta,
				lineData, lineRange, lineAngle, dataNum,
				BW, M, N, vMax, svmVal, linearVal, h+iMax-HOUGH_SIDE, 2*HOUGH_SIDE+1);
			current = current->next_;		
		}

		// get next one
		for (j=iMax-margin; j<=iMax+margin; j++)
		{
			if (j>=0 && j<rhoLen)
				h[j] = 0;
		}
		vMax = -1;
		iMax = -1;
		for (j=0; j<rhoLen; j++)
		{
			if (vMax<h[j])
			{
				vMax = h[j];
				iMax = j;
			}
		}
	}

	return num_added;
}

UWCW_INT rhoTrackerToCable(rtracker *rtrackers_, UWCW_INT totalNum, cableLine **cls)
{
	UWCW_INT num(0), i(0), k(0), dataNum;
	rtracker *current(rtrackers_);

	// get the number of cables
	while (current && !current->rho) current = current->next_;
	for (i=0; i<totalNum; i++)
	{
		if (current->rstatus == dormant)
			;
		else
			num++;
		
		current = current->next_;
	}

	// transfer data, only status and score are important
	*cls = (cableLine*)calloc(num, sizeof(cableLine));
	current = rtrackers_;
	while (current && !current->rho) current = current->next_;
	for (i=0; i<totalNum; i++)
	{
		if (current->rstatus == dormant)
			;
		else
		{
			(*cls)[k].isCable = true;
			(*cls)[k].score = current->rpts_[0].linearVal;
			(*cls)[k].rho = current->rho;
			(*cls)[k].theta = current->theta;
						
			k++;
		}

		current = current->next_;
	}

	return num;
}

UWCW_FLOAT getRhoThreshRaw(UWCW_FLOAT *RP_dBm, 
	UWCW_INT M, UWCW_INT N, UWCW_INT NC, 
	UWCW_FLOAT MagMax, UWCW_FLOAT MagMin)
{
	UWCW_INT interFactor = (UWCW_INT)((UWCW_FLOAT)M/512.0+0.5);
	int i;
	UWCW_FLOAT m;

	m = RP_dBm[0];
	for (i=0; i<(int)M*(int)N; i++)
	{
		if (m < RP_dBm[i])
			m = RP_dBm[i];
	}

	if (m > MagMax) m = MagMax;
	m -= MagMin;

	return m * (UWCW_FLOAT)interFactor * 3.;
}


UWCW_INT rhoTrack(UWCW_FLOAT theta, UWCW_INT *H, UWCW_FLOAT *rho, UWCW_INT rhoLen, UWCW_INT thetaLen,
	UWCW_FLOAT *RP_dBm, UWCW_INT *BW, UWCW_INT *BWcart, UWCW_FLOAT *Ang_plot, UWCW_FLOAT MagMax,
	UWCW_FLOAT MagMin, UWCW_INT M, UWCW_INT N, UWCW_INT NC, UWCW_INT x0,
	svmData *svm, UWCW_INT margin)
{
	UWCW_INT thetaI, idx, i, *h, num_trackers, num_trackers2, iMax0;
	UWCW_FLOAT theta_show = theta;
	UWCW_FLOAT thresh;

	// get theta index in H
	if (theta<0) theta+=180;
	thetaI = (UWCW_INT)(theta+0.5);
	idx;
	if (thetaI<PTHETA_END) thetaI = PTHETA_END;
	if (thetaI>180-PTHETA_END) thetaI = 180-PTHETA_END;
	if (thetaI>=90) idx = thetaI-90;
	else idx = thetaI-(90-PTHETA_LEN);

	// copy H(theta) into a separate array
	h = (UWCW_INT*)calloc(rhoLen, sizeof(UWCW_INT));
	for (i=0; i<rhoLen; i++)
	{
			h[i] = H[i*thetaLen+idx];
	}
	thresh = get_rho_thresh(h, rhoLen);
	// code for new threshold
	thresh = getRhoThreshRaw(RP_dBm, M, N, NC, MagMax, MagMin);
	PREV_RHO_THRESH = thresh;
	if (!theta_first_frame)		
	{
		thresh = (thresh + PREV_RHO_THRESH) / 2.0;
	}
	

	// process tracker list here
	process_rtrackers(g_rtrackers_, theta, h, rho, rhoLen, thetaLen, RP_dBm, BW, BWcart, Ang_plot, M, N, NC, x0, svm, margin);

	// add new trackers
	num_trackers = add_rtrackers(&g_rtrackers_, h, rhoLen, theta_show, rho, RP_dBm, BW, Ang_plot, M, N,x0, svm, margin, thresh);
	num_trackers2 = get_rtracker_num(g_rtrackers_);
	
	// show tracker image, output rho
	UWCW_FLOAT *rho_show = (UWCW_FLOAT*)calloc(num_trackers2, sizeof(UWCW_FLOAT));
	UWCW_INT *is_show = (UWCW_INT*)calloc(num_trackers2, sizeof(UWCW_INT));
	get_rtracker_rho(g_rtrackers_, rho_show, is_show);

	rhoNum = rhoTrackerToCable(g_rtrackers_, num_trackers2, &rhoCableLine);
	completeCableData(rhoCableLine, rhoNum, RP_dBm, Ang_plot, M, N, BW, x0);
	
	free(h);
	return num_trackers;
}

void completeCableData(cableLine *cables, UWCW_INT num,
	UWCW_FLOAT *RP_dBm, UWCW_FLOAT *Ang_plot, UWCW_INT M, UWCW_INT N,
	UWCW_INT *BW, UWCW_INT x0)
{
	UWCW_INT dataNum, i;
	UWCW_FLOAT lineData[512];
	UWCW_INT lineRange[512];
	UWCW_INT lineAngle[512];

	for (i=0; i<num; i++)
	{
		dataNum = getLineData(RP_dBm, Ang_plot, M, N, BW, x0, cables[i].rho, cables[i].theta, lineData, lineRange, lineAngle, true);
		cables[i].nPix = dataNum;
		
		cables[i].range = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_INT));
		cables[i].angle = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_INT));
		cables[i].data  = (UWCW_FLOAT*)calloc(dataNum, sizeof(UWCW_FLOAT));

		memcpy(cables[i].range, lineRange, dataNum*sizeof(UWCW_INT));
		memcpy(cables[i].angle, lineAngle, dataNum*sizeof(UWCW_INT));
		memcpy(cables[i].data, lineData, dataNum*sizeof(UWCW_FLOAT));
	}
}

UWCW_INT get_rtracker_rho(rtracker *rtracker_, UWCW_FLOAT *rho, UWCW_INT *is_cable)
{
	UWCW_INT num = 0;
	while (rtracker_)
	{
		if (fabs(rtracker_->rho)>0.01)
		{			
			*rho = rtracker_->rho;
			if (rtracker_->rstatus == healthy)
				*is_cable = 1;
			else if (rtracker_->rstatus == dormant)
				*is_cable = 2;
			else
				*is_cable = 0;
			num++;
			rho++;
			is_cable++;
		}
		rtracker_ = rtracker_->next_;
	}
	return num;
}


// identify peaks after HT
UWCW_INT houghPeaks(UWCW_INT *H,			// Hough data
				UWCW_INT rhoLen,			// # of rows in H
				UWCW_INT thetaLen,			// # of columns in H
				UWCW_INT numPeaks,			// # of peaks to find
				UWCW_FLOAT hGate,			// gate value for a peak
				UWCW_INT nHoodRowHalf,		// neighborhood size half row
				UWCW_INT nHoodColHalf,		// neighborhood size half column
				UWCW_INT *peaks)
{
	UWCW_INT num, done, maxRow, maxCol;
	UWCW_FLOAT maxVal;
	UWCW_INT m, n, nn, m1, m2;
	UWCW_INT M, N;
	M = rhoLen;
	N = thetaLen;

	done = 0;
	num = 0;
	while (!done)
	{
		getMaxVal(H, rhoLen, thetaLen, &maxVal, &maxRow, &maxCol);
		if (maxVal < hGate)
		{
			done = 1;
		}
		else
		{
			if (num==0 && maxVal>2*hGate)
			{
				hGate = maxVal/2;
			}

			peaks[2*num] = maxRow;
			peaks[2*num+1] = maxCol;
			num++;

			if (num == numPeaks)
			{
				done = 1;
				break;
			}
			
			// suppress neighborhood
			for (n=maxCol-nHoodColHalf; n<=maxCol+nHoodColHalf; n++)
			{
				if (n<0)
					nn = n+N;
				else if (n>=N)
					nn = n-N;
				else
					nn = n;

				m1 = maxRow-nHoodRowHalf;
				if (m1<0) m1=0;
				m2 = maxRow+nHoodRowHalf;
				if (m2>M-1) m2=M-1;

				for (m=m1; m<=m2; m++)
				{
					H[m*thetaLen+nn] = 0;
				}
			}
		}
	}

	return num;
}

// lines is allocated outside
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
				 UWCW_FLOAT ratio)
{
	UWCW_INT sliceStep;
	UWCW_INT sliceNum;
	UWCW_INT offset, n, y1, y2, i;
	UWCW_INT M, N;	// size of one slice
	UWCW_INT q;
	UWCW_FLOAT d;
	UWCW_INT *peaks;
	UWCW_INT totalNum, num;
	UWCW_INT nHoodRowHalf, nHoodColHalf;
	UWCW_FLOAT nHood1, nHood2, hGate;
	UWCW_INT maxR, maxC;

	UWCW_INT rhoLen;
	UWCW_FLOAT *rho;
	UWCW_INT *H;

	getMaxVal(BWcart, MC, NC, &hGate, &maxR, &maxC);
	hGate *= ratio;
	hGate *= HGATE_FACTOR;

	sliceStep = sliceHeight - sliceOverlap;
	sliceNum = (UWCW_INT)((UWCW_FLOAT)(MC-sliceStep)/(UWCW_FLOAT)sliceStep + 1.0);

	M = sliceHeight;
	N = NC;
	d = sqrt((UWCW_FLOAT)(M-1)*(M-1)+(UWCW_FLOAT)(N-1)*(N-1));
	q = (UWCW_INT)(d)+1;
	rhoLen = 2*q+1;
	rho = (UWCW_FLOAT*)calloc(rhoLen, sizeof(UWCW_FLOAT));
	H = (UWCW_INT*)calloc(rhoLen*thetaLen, sizeof(UWCW_INT));

	totalNum = 0;
	peaks = (UWCW_INT*)calloc(numPeaks*2, sizeof(UWCW_INT));
	for (n=0; n<sliceNum; n++)
	{
		offset = n*sliceStep;
		y1 = offset;
		y2 = offset + sliceHeight;
		if (y2>MC) y2=MC;

		// not in the region yet
		if (ystart>=y2 || yend<=y1) continue;

		// HT and related stuff, rho, theta, data allocation
		M = y2-y1;
		N = NC;
		d = sqrt((UWCW_FLOAT)(M-1)*(M-1)+(UWCW_FLOAT)(N-1)*(N-1));
		q = (UWCW_INT)(d)+1;
		rhoLen = 2*q+1;
		memset(H, 0, sizeof(UWCW_INT)*rhoLen*thetaLen);
		for (i=-q; i<=q; i++)
		{
			rho[i+q] = (UWCW_FLOAT)i;
		}

		houghTransform(BWcart, MC, NC, y1, y2, rho, rhoLen, theta, thetaLen, H);

		// debug code
		UWCW_FLOAT dumVal;
		UWCW_INT dum1, dum2;
		getMaxVal(H, rhoLen, thetaLen, &dumVal, &dum1, &dum2);

		// peaks identification, hGate and nHoodSize, return line
		nHood1 = (UWCW_FLOAT)rhoLen/100.;
		nHood2 = (UWCW_FLOAT)thetaLen/30.;
		nHoodRowHalf = (UWCW_INT)(nHood1);
		nHoodColHalf = __max((UWCW_INT)(nHood2), 1);
		//hGate = (double)nHoodRowHalf*10;

		num = houghPeaks(H, rhoLen, thetaLen, numPeaks, hGate, nHoodRowHalf, nHoodColHalf, peaks);

		// from peaks location to rho and theta value
		for (i=0; i<num; i++)
		{
			// rho
			lines[2*totalNum] = rho[ peaks[2*i] ];
			lines[2*totalNum] += (UWCW_FLOAT)(offset)*sin(theta[ peaks[2*i+1] ]*UW_PI/180);

			// theta
			lines[2*totalNum+1] = theta[ peaks[2*i+1] ];			
			totalNum++;
		}
	}

	if (peaks) free(peaks);
	if (rho) free(rho);
	if (H) free(H);

	return totalNum;
}

UWCW_INT removeRepLines(UWCW_FLOAT *cLs, UWCW_INT num)
{
	UWCW_INT *idx = (UWCW_INT*)calloc(num, sizeof(UWCW_INT));
	UWCW_FLOAT *cL2 = (UWCW_FLOAT*)calloc(num, sizeof(UWCW_FLOAT)*2);
	UWCW_INT i, j;
	UWCW_INT newNum;

	for (i=0; i<num; i++)
	{
		idx[i] = 1;
		cL2[2*i] = cLs[2*i];
		cL2[2*i+1] = cLs[2*i+1];
	}

	for (i=num-1; i>=0; i--)
	{
		for (j=i-1; j>=0; j--)
		{
			if (fabs(cLs[2*j+1]-cLs[2*i+1])<=1.0) // same line, theta
			{
				if (fabs(cLs[2*j]-cLs[2*i])<=5.0)				// rho
				{
					idx[i] = 0;
				}
			}
		}
	}

	newNum=0;
	for (i=0; i<num; i++)
	{
		if (idx[i])
		{
			cLs[2*newNum] = cL2[2*i];
			cLs[2*newNum+1] = cL2[2*i+1];
			newNum++;
		}
	}

	if (idx) free(idx);
	if (cL2) free(cL2);

	return newNum;
}

//detectCablesOld(RP_dBm, M, N, Ang_plot, svm, y1, y2, BW, BWcart, NC, x0);
void detectCablesOld(UWCW_FLOAT *RP_dBm, UWCW_INT M, UWCW_INT N,
	UWCW_FLOAT *Ang_plot, svmData *svm, UWCW_INT y1, UWCW_INT y2,
	UWCW_INT *BW, UWCW_INT *BWcart, UWCW_INT NC, UWCW_INT x0)
{
	UWCW_INT sliceSize, sliceOverlap, sliceLineNum, sliceStep, sliceNum, nLine, dataNum, nCable, i;
	UWCW_FLOAT *lines, linearVal, svmVal;
	cableLine *tmp;

	UWCW_FLOAT lineData[512];
	UWCW_INT lineRange[512];
	UWCW_INT lineAngle[512];
	UWCW_FLOAT f[FEAT_SIZE];

	sliceSize = (UWCW_INT)((double)M/24.0);
	sliceOverlap = (UWCW_INT)((double)sliceSize/4.0);
	sliceLineNum = 3*(UWCW_INT)(sliceSize/(double)(M/16.0)+0.5);
	sliceStep = sliceSize - sliceOverlap;
	sliceNum = (UWCW_INT)((double)(M-sliceStep)/(double)sliceStep + 1.0);
	lines = (UWCW_FLOAT*)malloc(sizeof(UWCW_FLOAT)*2*sliceNum*sliceLineNum);

	nLine = detectLines(BWcart, M, NC, y1, y2, THETA, THETA_LEN,
		sliceSize, sliceOverlap, sliceLineNum, lines, 1.0);
	nLine = removeRepLines(lines, nLine);

	tmp = (cableLine*)calloc(128, sizeof(cableLine));
	nCable = 0;
	for (i=0; i<nLine; i++)
	{
		dataNum = getLineData(RP_dBm, Ang_plot, M, N,
			BW, x0, lines[2*i], lines[2*i+1],
			lineData, lineRange, lineAngle, false);
		computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
		// svm
		linearVal = linearCalc(svm, f);
		svmVal = svmCalc(svm, f);

		if (dataNum>=15 && svmVal<0 && linearVal>=0.5)
		{
			tmp[nCable].isCable = true;
			tmp[nCable].score = linearVal;
			tmp[nCable].rho = lines[2*i];
			tmp[nCable].theta = lines[2*i+1];
			tmp[nCable].nPix = dataNum;
			tmp[nCable].data = (UWCW_FLOAT*)calloc(dataNum, sizeof(UWCW_FLOAT));
			tmp[nCable].range = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_INT));
			tmp[nCable].angle = (UWCW_INT*)calloc(dataNum, sizeof(UWCW_INT));
			memcpy(tmp[nCable].data, lineData, dataNum*sizeof(UWCW_FLOAT));
			memcpy(tmp[nCable].range, lineRange, dataNum*sizeof(UWCW_INT));
			memcpy(tmp[nCable].angle, lineAngle, dataNum*sizeof(UWCW_INT));
			nCable++;
		}
	}

	voteTheta(tmp, nCable);

	rhoNum = 0;
	for (i=0; i<nCable; i++)
	{
		if (tmp[i].isCable)
			rhoNum++;
	}
	if (rhoNum)
		rhoCableLine = (cableLine*)calloc(rhoNum, sizeof(cableLine));

	rhoNum = 0;
	for (i=0; i<nCable; i++)
	{
		if (tmp[i].isCable)
			rhoCableLine[rhoNum++] = tmp[i];
	}

	if (tmp) free(tmp);
	if (lines) free(lines);
}

void voteTheta(cableLine *cables, UWCW_INT cableNum)
{
	int *countTheta = NULL;
	countTheta = (int*)calloc(THETA_LEN, sizeof(int));
	int i;
	int idx;
	int maxIdx;
	int maxCount, maxCount2;
	double maxTheta, maxTheta2;
	bool mulMax, secondMax;

	// vote
	for (i=0; i<cableNum; i++)
	{
		if (cables[i].isCable && cables[i].score>.5)
		{
			// original position: +6
			idx = thetaToIdx2(cables[i].theta);
			countTheta[idx] += 6;

			// +/-1: +2
			idx = thetaToIdx2(cables[i].theta+1);
			countTheta[idx] += 2;
			idx = thetaToIdx2(cables[i].theta-1);
			countTheta[idx] += 2;

			// +/-2: +1
			idx = thetaToIdx2(cables[i].theta+2);
			countTheta[idx] += 1;
			idx = thetaToIdx2(cables[i].theta-2);
			countTheta[idx] += 1;
		}
	}

	// find out the maximum; test if multiple max
	mulMax = false;
	maxIdx = 0;
	maxCount = -1;
	for (i=0; i<THETA_LEN; i++)
	{
		if (countTheta[i] > maxCount)
		{
			maxCount = countTheta[i];
			maxIdx   = i;
			mulMax = false;
		}
		else if(countTheta[i] == maxCount)
		{
			mulMax = true;
		}
	}

	// if multiple max, do nothing
	if (!mulMax)
	{
		maxTheta = THETA[maxIdx];
		
		// get the second maximum count theta
		countTheta[maxIdx] = 0;
		countTheta[(maxIdx-1+THETA_LEN)%THETA_LEN] = 0;
		countTheta[(maxIdx+1)%THETA_LEN] = 0;

		maxIdx = 0;
		maxCount2 = -1;
		for (i=0; i<THETA_LEN; i++)
		{
			if (countTheta[i] > maxCount2)
			{
				maxCount2 = countTheta[i];
				maxIdx = i;
			}
		}

		// if second maximum is large enough, keep it
		maxTheta2 = 0;
		if (maxCount2 >= int(double(maxCount)*0.5))
		{
			secondMax = true;
			maxTheta2 = THETA[maxIdx];
		}
		else
			secondMax = false;

		// screen out thetas not in maxTheta or maxTheta2
		for (i=0; i<cableNum; i++)
		{
			if (cables[i].isCable && cables[i].score>.5)
			{
				if (fabs(cables[i].theta-maxTheta)>1)
					if (!secondMax || fabs(cables[i].theta-maxTheta2)>1)
						cables[i].isCable = false;
			}
		}
	}

	free(countTheta);
}

UWCW_INT thetaToIdx2(UWCW_FLOAT theta)
{
	UWCW_INT thetaI = (UWCW_INT)theta;
	UWCW_INT idx;

	if (thetaI<0)
		idx = thetaI + 90;
	else
		idx = thetaI - 29;

	if (idx<0 || idx>= THETA_LEN)
		idx = (idx+THETA_LEN) % THETA_LEN;

	return idx;
}

void saveCableData(cableLine *cls, UWCW_INT num)
{
	int i, j;
	char name[256];

	for (i=0; i<num; i++)
	{
		sprintf(name, "data_%d.txt", i);
		FILE *file = fopen(name, "w");
		fprintf(file, "rho = %f, theta = %f\n", cls[i].rho, cls[i].theta);
		for (j=0; j<cls[i].nPix; j++)
		{
			fprintf(file, "%f\n", cls[i].data[j]);
		}
		fclose(file);
	}
}

void saveLineData(UWCW_FLOAT *lineData, UWCW_INT dataNum)
{
	FILE *file = fopen("lineData.txt", "w");
	UWCW_INT i;

	for (i=0; i<dataNum; i++)
	{
		fprintf(file, "%f\n", lineData[i]);
	}

	fclose(file);
}