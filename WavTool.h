/*
 * WavTool.h
 *
 *  Created on: Jan 11, 2013
 *      Author: panda
 */

#ifndef WAVTOOL_H_
#define WAVTOOL_H_
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

typedef struct  WAV_HEADER
   {
     char                RIFF[4];        /* RIFF Header      */ //Magic header
     unsigned long       ChunkSize;      /* RIFF Chunk Size  */
     char                WAVE[4];        /* WAVE Header      */
     char                fmt[4];         /* FMT header       */
     unsigned long       Subchunk1Size;  /* Size of the fmt chunk                                */
     unsigned short      AudioFormat;    /* Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM */
     unsigned short      NumOfChan;      /* Number of channels 1=Mono 2=Sterio                   */
     unsigned long       SamplesPerSec;  /* Sampling Frequency in Hz                             */
     unsigned long       bytesPerSec;    /* bytes per second */
     unsigned short      blockAlign;     /* 2=16-bit mono, 4=16-bit stereo */
     unsigned short      bitsPerSample;  /* Number of bits per sample      */
     char                Subchunk2ID[4]; /* "data"  string   */
     unsigned long       Subchunk2Size;  /* Sampled data length    */
   }wav_hdr;
typedef struct  WAV_INFO
{
	wav_hdr header;
	short *samples;
	long nsample;
}wav_inf;


int getFileSize(FILE *inFile);


void Wavtool_printinfo(char *filepath);
wav_hdr* Wavtool_read(char *filepath);
wav_inf* Wavtool_readSam(char *filepath);

#endif /* WAVTOOL_H_ */
