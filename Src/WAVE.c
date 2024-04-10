/**
 * waveFunctions.c
 *
 * Created on: Feb 24, 2022
 *
 * Author: Michael Noyce
 */

/*User Includes-------------------------------------------------*/
#include "WAVE.h"
#include "HAL_SD.h"
#include "arm_math.h"
/*End User Includes--------------------------------------------*/

/*Private Variables---------------------------------------------*/

//LPF Filter variables
static float32_t firStateF32[BLOCK_SIZE + NUM_TAPS - 1];
const float32_t firCoeffs32[NUM_TAPS_ARRAY_SIZE] = {
		0.00251240944810394,0.00312707878616391,0.00451235694993659,0.00680393428856072,0.0100755541121319,0.0143279617567923,0.0194834494742694,0.0253866652999644,0.0318118288415239,0.0384759479928208,0.0450571065900547,0.051216441065083,0.0566220853690755,0.060973169244554,0.0640219240256733,0.0655920867552915,0.0655920867552915,0.0640219240256733,0.060973169244554,0.0566220853690755,0.051216441065083,0.0450571065900547,0.0384759479928208,0.0318118288415239,0.0253866652999644,0.0194834494742694,0.0143279617567923,0.0100755541121319,0.00680393428856072,0.00451235694993659,0.00312707878616391,0.00251240944810394
};

/*----------------------------------------------------------------*/

/**
 * @brief fullPipeline
 * @param freqSegments
 * @param PSDOutput
 */
void fullPipeline(Wave_Data_t waveData, uint32_t WaveDirNo)
{

	float32_t PSD[FFT_LENGTH];

	//Welch Method - calls single segment pipeline in function
	WelchMethod(PSD, WaveDirNo);
	/*for(int i=0; i<512; i++)
	  				    {
	  				    	printmsg("%.2f \r\n", PSD[i]);
	  				    	PSDOutput[i] = PSD[i];
	  				    	printmsg("PSD sample %d \r\n", i);
	  				    } */

	//Wave Extraction Algorithm
	//PSD Spectra input - wave parameter output
	//waveParamExtract(PSD, &Hm0, &Hrms, &T0, &Tp, &M0, &M1, &M2);
	waveParamExtract(PSD, &waveData.Hm0, &waveData.TM01, &waveData.Hrms, &waveData.T0, &waveData.Tp, &waveData.M0, &waveData.M1, &waveData.M2);


	/*for(int i = 0; i<(FFT_LENGTH/2); i++)
	{
		if(i<(FFT_LENGTH/2)

	} */

	printmsg("Hm0 %.8f \r\n", waveData.Hm0);
	printmsg("Hrms %.8f \r\n", waveData.Hrms);
	printmsg("M0 %.8f \r\n", waveData.M0);
	printmsg("M1 %.8f \r\n", waveData.M1);
	printmsg("M2 %.8f \r\n", waveData.M2);
	printmsg("Tm01 %.8f \r\n", waveData.TM01);
	printmsg("T0 %.8f \r\n", waveData.T0);
	printmsg("Tp %.8f \r\n", waveData.Tp);


}

/**
 * @brief singleSegmentPipeline
 * @param rawData
 * @param decArray
 */
void singleSegmentPipeline(float32_t* singleSegment, uint32_t waveDirNo, uint32_t waveLogNo)
{

    uint32_t count;
    uint32_t i;
    uint32_t j;

    int32_t rawData[FFT_LENGTH];
    float32_t smallBuf[FFT_LENGTH/DECIMATION_CONSTANT];
    float32_t tempBuf[FFT_LENGTH];
    float32_t tempBuf2[FFT_LENGTH];
    float32_t tempBuf3[FFT_LENGTH];
    float32_t tempBuf4[FFT_LENGTH/2];
    float32_t tempBuf5[FFT_LENGTH];
    float32_t frequencyOutput[FFT_LENGTH/2];

    //get FFT_LENGTH raw data and decimate to FFT_LENGTH/DECIMATION_CONSTANT
    //append until array equals FFT_LENGTH


    uint32_t fpointer = 0;


    for(count = 0; count<DECIMATION_CONSTANT; count++)
    {
    	/*----------Read in data here----------------------*/

    	SD_Wave_Read(&File, rawData, waveDirNo, waveLogNo, Z_ACC, &fpointer);

    	/*----------Read in data here----------------------*/

    	calibrate(rawData, tempBuf);

    	detrend(tempBuf, tempBuf2);

    	LPFDecimate(tempBuf2, smallBuf);

    	//Append detrended, calibrated and decimated array

        for (j = 0; j < FFT_LENGTH/DECIMATION_CONSTANT ; j++)
        {
        	tempBuf3[count*(FFT_LENGTH)/(DECIMATION_CONSTANT) + j] = smallBuf[j];
        }

    }

    embeddedFFT(tempBuf3, frequencyOutput, FFT_LENGTH, FFT);


    HPFDoubleIntegration(frequencyOutput, tempBuf4);


    //Output Inverse - Wave Amplitude series
    embeddedFFT(tempBuf4, tempBuf5, FFT_LENGTH, IFFT);


    //Apply Hanning Window to reduce spectral leakage
   float32_t tempWindow;

    for(i = 0; i<FFT_LENGTH; i++)
    {
    	tempWindow = 0.5 * (1 - cos(2*PI*i/FFT_LENGTH));
    	singleSegment[i] = tempBuf5[i]*tempWindow;
    }

}


/**
 * @brief Calibration function for the ICM20649
 * @param rawData
 * @param calOutput
 */
void calibrate(int32_t* rawData, float32_t* calOutput)
{
	uint32_t i;

	for (i = 0; i < FFT_LENGTH ; i++)
	{
		calOutput[i] = ((float32_t)rawData[i] - 92.22)/1686.50 -9.81; //remove bias offset, scale and subtract g
	}

}

/**
 * @brief Detrending function
 * @param rawData
 * @param detrendOutput
 */
void detrend(float32_t* calOutput, float32_t* detrendOutput)
{

	//Determining Mean of data
	uint32_t i;
	float32_t sumTotal, dataMean;

	for (i= 0; i < FFT_LENGTH; i++)
	{
		sumTotal = sumTotal + calOutput[i];
	}

	dataMean = sumTotal/FFT_LENGTH;

	for (i= 0; i < FFT_LENGTH; i++)
	{
		detrendOutput[i] = calOutput[i] - dataMean;
	}
	// RC Highpass function	for detrending (moving average filter)
	float32_t k = 0.9995;
	float32_t s[FFT_LENGTH];
	float32_t detrendTemp[FFT_LENGTH];


	s[0] = 0;

	for (i= 1; i < FFT_LENGTH; i++)
	{
		 s[i] = detrendTemp[i] + k*(s[i-1]);
		 detrendTemp[i] = detrendTemp[i] - (1 - k)*(s[i]);
	}

	detrendOutput = detrendTemp;

}

/**
 * @brief LPFDecimate
 * @param testInput
 * @param testOutput
 */
void LPFDecimate(float32_t* testInput, float32_t* testOutput)
{

	uint32_t blockSize = BLOCK_SIZE;
	uint32_t numBlocks = FFT_LENGTH/BLOCK_SIZE;
	arm_fir_decimate_instance_f32 S;
	uint32_t i;
	float32_t  *inputF32, *outputF32;

	/* Initialize input and output buffer pointers */
	inputF32 = &testInput[0];
	outputF32 = &testOutput[0];

    arm_fir_decimate_init_f32(&S, NUM_TAPS, DECIMATION_CONSTANT, (float32_t *)&firCoeffs32[0], &firStateF32[0], blockSize);

    for(i=0; i < numBlocks; i++)
    {
      arm_fir_decimate_f32(&S, inputF32 + (i * blockSize), outputF32 + (i * blockSize/DECIMATION_CONSTANT), blockSize);
    }

}

/**
 * @brief embeddedFFT
 * @param inputArray
 * @param frequencyOutput
 */
void embeddedFFT(float32_t* inputArray,float32_t* frequencyOutput, uint32_t FFTLength, uint32_t ifftFlag)
{

	arm_rfft_fast_instance_f32 varInstRfftF32;

	arm_rfft_fast_init_f32(&varInstRfftF32, FFTLength);

	 //Process the data through the CFFT/CIFFT module
	 arm_rfft_fast_f32(&varInstRfftF32, inputArray, frequencyOutput, ifftFlag);

	 //Process the data through the Complex Magnitude Module for calculating the magnitude at each bin
	 //arm_cmplx_mag_f32(frequencyOutput+2, magnitudeOutput+1, FFTLength/2-1);

	 // DC component - consider removing if still conservative

}

/**
 * @brief HPFDoubleIntegration
 * @param fftInput
 * @param waveAmplitude
 */
void HPFDoubleIntegration(float32_t* rfftInput, float32_t* waveAmplitude)
{
	float32_t freq1 = 0.02;
	float32_t freq2 = 0.03;
	float32_t fnyq = F_SAMPLE/DECIMATION_CONSTANT/2;
	uint32_t i;
	float32_t Fs = F_SAMPLE/DECIMATION_CONSTANT;
	float32_t freq[FFT_LENGTH/2];



	//Create frequency array
	for(i=0; i<FFT_LENGTH/2; i++)
	{
		freq[i] = Fs*((float32_t)i/(FFT_LENGTH/2));
	}


	//HPF and double integration in frequency domain
	//Based on Tucker and Pitt
	for(i=0;i<FFT_LENGTH/2; i++)
	{
		if(freq[i]<freq1)
		{
			waveAmplitude[i] = 0.0f;
		}

		if(freq[i]>=freq1 && freq[i]<freq2)
		{
			waveAmplitude[i] = (1/2)*(1-cos(PI*(freq[i]-freq1)/(freq2-freq1))*(-1/(2*PI*(freq[i])*(rfftInput[i]))));
		}

		if(freq[i]>=freq2 && freq[i]<fnyq)
		{
			waveAmplitude[i] = -1/(2*PI*((freq[i])*(freq[i])))*rfftInput[i];
		}

		if(isnan(waveAmplitude[i]))
		{
			waveAmplitude[i] = 0.0f;
		}

	}

}


//Welch Method

void WelchMethod(float32_t* PSD, uint32_t waveDirNo)
{
	uint32_t i;
	uint32_t j;
	float32_t singleSegment[FFT_LENGTH];
	float32_t singleSegmentFFT[FFT_LENGTH/2];
	float32_t singleSegmentMag[FFT_LENGTH/2];
	float32_t tempSum;
	float32_t waveAmplitudeArray[SEGMENT_NO][FFT_LENGTH/2];
	float32_t normalConst = (float32_t)(FFT_LENGTH * F_SAMPLE); //normalizes periodogram estimate

	//Populate wave amplitude array with windowed FFT functions

	uint32_t waveLogNo = 1;

	for(i = 1; i<SEGMENT_NO; i++)
	{
		//populate singleSegment array
		singleSegmentPipeline(singleSegment, waveDirNo, waveLogNo);
		waveLogNo++;
		//FFT each segment
		embeddedFFT(singleSegment, singleSegmentFFT, FFT_LENGTH, FFT);
		arm_cmplx_mag_f32(singleSegmentFFT+2, singleSegmentMag+1, FFT_LENGTH/2-1);
		//Remove frequencies above 2 Hz 
		int index_1Hz = (1 * FFT_LENGTH) / ( F_SAMPLE / DECIMATION_CONSTANT);
		// Set all frequency bins beyond index_1Hz to zero
		for (int i = index_1Hz + 1; i < FFT_LENGTH / 2; i++) {
			singleSegmentMag[i] = 0.0f;  // Assuming outputF32 is your RFFT result array
		}
		int index_0_03Hz = (0.05 * FFT_LENGTH) / ( F_SAMPLE / DECIMATION_CONSTANT);
		// Set all frequency bins bleow index_0_03Hz to zero
		for (int i = 0; i < index_0_03Hz; i++) {
			singleSegmentMag[i] = 0.0f;  // Assuming outputF32 is your RFFT result array
		}
		//Fill waveAmplitudeArray with segments
		for(j = 0; j<FFT_LENGTH/2; j++)
		{
			waveAmplitudeArray[i][j] = singleSegmentMag[j];

		}
	}

	//Periodogram

	for (i = 0; i < SEGMENT_NO ; i++)
	{
		for (j = 0; j < FFT_LENGTH/2; j++)
		{
			waveAmplitudeArray[i][j] = (waveAmplitudeArray[i][j]*waveAmplitudeArray[i][j])/(normalConst);
		}

	}

	// Averaging
	for (i = 0; i < FFT_LENGTH/2 ; i++)
	{
		tempSum = 0;
		for (j = 0; j < SEGMENT_NO ; j++)
		{
			tempSum += waveAmplitudeArray [j][i];
		}
		PSD [i] = tempSum / SEGMENT_NO * FFT_LENGTH ;
	}

	PSD [0] /= 2;

}

//Wave Parameter extraction
//Based on SWASH Model
void waveParamExtract(float32_t* PSD, float32_t* Hm0, float32_t* Tm01, float32_t* Hrms, float32_t* T0, float32_t* Tp, float32_t*  M0, float32_t* M1, float32_t* M2)
{

	//Iterators
	uint32_t i, j;
	//Length of PSD input
	float32_t sample_length = FFT_LENGTH/2;
	//Subint Length
	float32_t subInt_length;
	// Number of subintervals
	uint32_t subInt_no = sample_length - 1;
	//Integral
	float32_t Sum = 0;
	//Frequency array
	float32_t freq[FFT_LENGTH/2];
	//Max Amplitude freq
	float32_t fm;
	//Freq index of max amplitude
	uint32_t fi;
	//temp moment values
	float32_t m0, m1, m2;

	//Create frequency array
	for(i=0; i<FFT_LENGTH/2; i++)
	{
		freq[i] = ((float32_t)F_SAMPLE/(float32_t)DECIMATION_CONSTANT/2)*((float32_t)i/(FFT_LENGTH));
	}



	//Length of subinterval
	subInt_length=(freq[subInt_no+1]-freq[1])/subInt_no;

	//Newton Cotes Quadrature to calculate spectral moments
	//Note scaling by 10000 used to preserve resolution

	for(i=2; i<FFT_LENGTH/2; i++)
	{
		Sum = Sum + PSD[i];
	}

	m0 = 10000*0.5*subInt_length*(PSD[1]+2*Sum+PSD[subInt_no+1]);

	*M0 = m0/10000;

	printmsg("M0 calculated!");

	Sum = 0;

	for(i=2; i<FFT_LENGTH/2; i++)
	{
		Sum = Sum + (freq[i])*PSD[i];
	}

	printmsg("f^2 sum %.8f", Sum);

	m1 = 10000*0.5*subInt_length*(freq[1]*PSD[1]+2*Sum+(freq[subInt_no+1]*PSD[subInt_no+1]));

	*M1 = m1/10000;

	printmsg("M1 calculated!: m1 %.8f", m1);

	Sum = 0;

	for(i=2; i<FFT_LENGTH/2; i++)
		{
			Sum = Sum + (freq[i]*freq[i])*PSD[i];
		}

	 m2 = 10000*0.5*subInt_length*(freq[1]*freq[1]*PSD[1]+2*Sum+(freq[subInt_no+1]*freq[subInt_no+1]*PSD[subInt_no+1]));

	 *M2 = m2/10000;

	printmsg("M2 calculated!");

	//add M2, T0, Hrms and H0 calculations

	*Hm0 = 4*sqrt(fabsf(m0/10000)); //correct
	*Tm01 = (fabsf(m0/m1)); //correct
	*Hrms = sqrt(2)/2*(*Hm0); // Correct
	*T0 = sqrt(fabsf(m2)/fabsf(m0)); //Correct

	printmsg("Parameters calculated!");

	 for (i = 1; i < FFT_LENGTH/2; ++i) {
	    if (PSD[0] < PSD[i]) {
	    	fi = i;
	    }
	  }

	fm = freq[fi];

	*Tp = 1/fm;

	printmsg("Tp calculated!");

}





