/* 
 * Copyright (c) 2012.
 * This file is part of ALOE (http://flexnets.upc.edu/)
 * 
 * Eq.ZF a& Eq.MMSE designed by: Rubén Aldana, Sandra Dacosta & Berta García
 * ALOE++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ALOE++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ALOE++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <phal_sw_api.h>
#include "skeleton.h"
#include "params.h"

#include "UPLINK_EQUALIZER_interfaces.h"
#include "UPLINK_EQUALIZER_functions.h"
#include "UPLINK_EQUALIZER.h"

//ALOE Module Defined Parameters. Do not delete.
char mname[STR_LEN]="UPLINK_EQUALIZER";

//Module User Defined Parameters


//Global Variables
int numFFTs=14;
int DATAsize=156;
static int cont;
//CODIGO AÑADIDO----------------------------------------
//  ZC sequence Añadido
#define BUFFER_SZ	2048
_Complex float DMRS_SEQ0[156];
_Complex float DMRS_SEQ1[156];
_Complex float COEFF0[156];
_Complex float COEFF1[156];
_Complex float H[156];
static _Complex float Htot[156];
_Complex float F[156];
int DMRS_length=0;

float SNRm;

float pSignal, pNoise;

float pSignal0, pSignal1, Pnoise0, Pnoise1, SNR; 

int M_RS_SC = 156;
//CODIGO AÑADIDO-------------------------------------------

/*
 * Function documentation
 *
 * @returns 0 on success, -1 on error
 */
int initialize() {

	/* Get control parameters*/

	/* Verify control parameters */

	/* Print Module Init Parameters */
	strcpy(mname, GetObjectName());
	printf("O--------------------------------------- Eq.MMSE-SNR ----------------------------------------O\n");
	printf("O    SPECIFIC PARAMETERS SETUP: \033[1;34m%s\033[0m\n", mname);
	printf("O      Nof Inputs=%d, DataTypeIN=%s, Nof Outputs=%d, DataTypeOUT=%s\n", 
		       NOF_INPUT_ITF, IN_TYPE, NOF_OUTPUT_ITF, OUT_TYPE);
	printf("O--------------------------------------------------------------------------------------------O\n");

	/* do some other initialization stuff */
	//CODIGO AÑADIDO-------------------------------------------
	// CALCULATE DMRS sequence 
	M_RS_SC = 156;
	DMRS_length=genRSsignalargerThan3RB(1, 0, 10, M_RS_SC, DMRS_SEQ0, 0);
	DMRS_length=genRSsignalargerThan3RB(0, 1, 10, M_RS_SC, DMRS_SEQ1, 0);
	int k;
	cont=1;
	//CODIGO AÑADIDO-------------------------------------------	


	return 0;
}



/**
 * @brief Function documentation
 *
 * @param inp Input interface buffers. Data from other interfaces is stacked in the buffer.
 * Use in(ptr,idx) to access the address. To obtain the number of received samples use the function
 * int get_input_samples(int idx) where idx is the interface index.
 *
 * @param out Input interface buffers. Data to other interfaces must be stacked in the buffer.
 * Use out(ptr,idx) to access the address.
 *
 * @return On success, returns a non-negative number indicating the output
 * samples that should be transmitted through all output interface. To specify a different length
 * for certain interface, use the function set_output_samples(int idx, int len)
 * On error returns -1.
 *
 * @code
 * 	input_t *first_interface = inp;
	input_t *second_interface = in(inp,1);
	output_t *first_output_interface = out;
	output_t *second_output_interface = out(out,1);
 *
 */
int work(input_t *inp, output_t *out) {

//añadido para graph//----------------------------
	output_t *output1 = out(out,1);
	output_t *output2 = out(out,2);
//añadido para graph//------------------------------
	int rcv_samples = get_input_samples(0); /** number of samples at itf 0 buffer */
	int snd_samples=0;
	int i,j,k;
	

	if(rcv_samples == 0)return(0);

	// OUTPUT 1
	memcpy(output1, DMRS_SEQ0, sizeof(_Complex float)*DATAsize);
	int snd_samples1=156;

	// OUTPUT 2
	memcpy(output2, &inp[DATAsize*3], sizeof(_Complex float)*DATAsize);
	int snd_samples2=156;

	// Indicate the number of samples at output number N
	set_output_samples(1, snd_samples1);
	set_output_samples(2, snd_samples2);
						//EQUALIZER MMS-SNR
	for(k=0;k<DATAsize;k++){
		COEFF0[k] = inp[(DATAsize*3)+k]/DMRS_SEQ0[k]; // [3.1.1.1] 
		COEFF1[k] = inp[(DATAsize*10)+k]/DMRS_SEQ1[k];// [3.1.1.2] 
		F[k] = (COEFF0[k] + COEFF1[k])/2; // [3.1.1.3] 
		pSignal0 = (float)((inp[(DATAsize*3)+k]) * (conj(inp[(DATAsize*3)+k]))); // [3.1.3.3]
		pSignal1 = (float)((inp[(DATAsize*10)+k]) * (conj(inp[(DATAsize*10)+k]))); //[3.1.3.4]
		Pnoise0 = (float)((inp[(DATAsize*3)+k]-(DMRS_SEQ0[k]*F[k]))*conj(inp[(DATAsize*3)+k]-(DMRS_SEQ0[k]*F[k]))); // [3.1.2.1] & [3.1.3.1]
		Pnoise1 = (float)((inp[(DATAsize*10)+k]-(DMRS_SEQ1[k]*F[k]))*conj(inp[(DATAsize*10)+k]-(DMRS_SEQ1[k]*F[k]))); // [3.1.2.2] & [3.1.3.2]
		pSignal = ((pSignal0+pSignal1)/2);//[3.1.3.6]
		pNoise = ((Pnoise0+Pnoise1)/2); // [3.1.3.5]
		SNR = (((pSignal0-Pnoise0)/Pnoise0)+((pSignal1-Pnoise1)/Pnoise1))/2; // [3.1.4.1]
		SNRm = SNRm + SNR; // [3.1.4.2]
	}
		SNRm = SNRm/156; // [3.1.4.2]
	for(k=0;k<DATAsize;k++) //para cada portadora
	{		
		COEFF0[k] = inp[(DATAsize*3)+k]/DMRS_SEQ0[k]; // [3.1.1.1] 
		COEFF1[k] = inp[(DATAsize*10)+k]/DMRS_SEQ1[k];// [3.1.1.2] 
		F[k] = (COEFF0[k] + COEFF1[k])/2; // [3.1.1.3] 
		H[k] = (conj(F[k])/(pow(F[k],2)+(1/SNRm))); // [3.1.5.1]
		Htot[k] = Htot[k] + H[k]; //Hacemos la suma de los factores para ver su comportamiento general
	}
	memcpy(output1, Htot, sizeof(_Complex float)*DATAsize); //output para enviar los coefs de corr. a la graph
	cont++;
	printf("\nSNR media medida = %fdB\n",(10*log(SNRm/DATAsize)));

	// ELIMINATE PUSCH REFERENCE SIGNAL
	j=0;
	while(j<DATAsize)
	{
		for(i=0;i<numFFTs;i++)
		{
			inp[(DATAsize*i)+j] = inp[(DATAsize*i)+j]*(Htot[j]/cont);
		}
		j++;
	}
	j=0;
	for(i=0; i<numFFTs; i++){
		if(i!=3 && i!=10){
			memcpy(&out[DATAsize*j], &inp[DATAsize*i], sizeof(_Complex float)*DATAsize);
			snd_samples += DATAsize;
			j++;
		}
	}
	// Indicate the number of samples at output 0 with return value
	return snd_samples;
}

/** @brief Deallocates resources created during initialize().
 * @return 0 on success -1 on error
 */
int stop() {
	return 0;
}


