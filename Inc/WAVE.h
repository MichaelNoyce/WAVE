/**
 * waveFunctions.h
 *
 * Created on: Feb 24, 2022
 *
 * Author: Michael Noyce
 */

#ifndef HAL_WAVE_H_
#define HAL_WAVE_H_


/*User Includes------------------------------------------------*/
#include "arm_math.h"
#include "string.h"
/*End User Includes-------------------------------------------*/

/*Definitions-------------------------------------------------*/
//Input Data Definitions
#define RAW_DATA_LENGTH 5120
#define FFT_LENGTH 1024
#define SEGMENT_NO 5
#define DECIMATION_CONSTANT 16
#define F_SAMPLE 100
#define IFFT 1
#define FFT 0
//LPF Decimate Defintions
#define SNR_THRESHOLD_F32    75.0f
#define BLOCK_SIZE            32
/* Must be a multiple of 16 */
#define NUM_TAPS_ARRAY_SIZE              32
#define NUM_TAPS              32

/*End Definitions ---------------------------------------------*/



/* Structs ---------------------------------------------------*/

/*
 * @brief: Wave data structure
 */
typedef struct
{
	float32_t Hm0; //Significant waveheight
	float32_t Hrms; //Root mean square waveheight
	float32_t T0; //Mean zero crossing period
	float32_t TM01; //Mean spectral wave period
	float32_t Tp; // Peak wave period
	float32_t M0; //0th Spectral moment
	float32_t M1; //1st spectral moment
	float32_t M2; //2nd spectral moment
	float32_t PSD[30]; //PSD Period in spectral bins

}Wave_Data_t;




/* Wave Parameter Extraction functions ------------------------------------------------------------*/
/**
 * @brief
 *
 * @param waveData
 */
void fullPipeline(Wave_Data_t waveData, uint32_t WaveDirNo);

/**
 *
 * @param rawData
 * @param decArray
 */
void singleSegmentPipeline(float32_t* singleSegment,  uint32_t waveDirNo, uint32_t waveLogNo);

/**
 * @brief Calibration Function for the ICM20689
 * @param rawData
 * @param calOutput
 */
void calibrate(int32_t* rawData, float32_t* calOutput);

/**
 * @brief Detrending function
 * @param rawData
 * @param detrendOutput
 */
void detrend(float32_t* calOutput, float32_t* detrendOutput);

/**
 * @brief Decimating input low pass filtering function
 * @param detrendOutput
 * @param decimatedOutput
 */
void LPFDecimate(float32_t* detrendOutput, float32_t* decimatedOutput);

/**
 * @brief Performs Four Step FFT method on time series greater than 1024 (allows use of CMSIS DSP)
 * @param timeSeries
 * @param frequencyOutput
 */
void embeddedFFT(float32_t* timeSeries, float32_t* frequencyOutput, uint32_t fftLength, uint32_t ifftFlag);

/**
 * @brief Double Integrating time series
 * @param rfftInput
 * @param waveAmplitude
 */
void HPFDoubleIntegration(float32_t* rfftInput, float32_t* waveAmplitude);


/**
 * @brief WelchMethod
 * @param waveAmplitudeArray
 * @param waveAmplitude
 */
void WelchMethod(float32_t* PSD, uint32_t waveDirNo);


/**
 * @brief waveParamExtract
 * @param PSD
 * @param Hm0
 * @param Hrms
 * @param T0
 * @param Tp
 * @param M0
 * @param M1
 * @param M2
 */
void waveParamExtract(float32_t* PSD, float32_t* Hm0, float32_t* Tm01, float32_t* Hrms, float32_t* T0, float32_t* Tp, float32_t*  M0, float32_t* M1, float32_t* M2);

#endif /* HAL_WAVE_H_ */
