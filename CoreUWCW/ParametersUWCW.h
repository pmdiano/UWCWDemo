#include "DefUWCW.h"

// ********************************************************************************** //
// general parameters
// ********************************************************************************** //
#define UW_PI			3.14159265358979323846
#define CABLE_TIP_END 15	/** number of pixels to exclude on each tip */

// ********************************************************************************** //
// thresholding parameters
// ********************************************************************************** //
extern UWCW_FLOAT TH_PCT;
extern UWCW_INT HGATE_FACTOR;
extern UWCW_FLOAT TH_GATE;
extern UWCW_INT SLICE_HEIGHT;
extern UWCW_INT SLICE_WIDTH;
extern UWCW_INT wFactor;
extern UWCW_INT hFactor;

// ********************************************************************************** //
// feature parameters
// ********************************************************************************** //
#define SPEC_SIZE		2048			// power spectrum DFT size
#define AUTO_SIZE		128				// autocorrelation DFT size
#define FEAT_SIZE		14				// feature size
extern UWCW_FLOAT LOWBOUND_PCT;

// ********************************************************************************** //
// hough transform parameters
// ********************************************************************************** //
extern UWCW_FLOAT PTHETA[];
extern UWCW_INT PTHETA_LEN;
extern UWCW_INT PTHETA_END;
extern UWCW_FLOAT THETA[];
extern UWCW_INT THETA_LEN;

// ********************************************************************************** //
// theta particle filter parameters
// ********************************************************************************** //
extern UWCW_INT THETA_NUM_PARTICLES;
extern UWCW_FLOAT THETA_VARIANCE;
extern UWCW_INT THETA_SIDE;

// ********************************************************************************** //
// rho particle filter parameters
// ********************************************************************************** //
#define RHO_NUM_PARTICLES 30
#define SVM_BASE 0.3
#define ASSOCIATION_THRESHOLD 0.4
#define DORMANCY_THRESHOLD 1
#define HOUGH_SIDE 3
extern UWCW_FLOAT RHO_VARIANCE;
#define MAX_TRACKERS 8