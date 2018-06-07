/*
 * Wav_tool.c
 *
 *  Created on: Jan 11, 2013
 *      Author: panda
 */
#include "WavTool.h"
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <QFile>
int getFileSize(FILE *inFile)
{
 int fileSize = 0;
    fseek(inFile,0,SEEK_END);
    fileSize=ftell(inFile);
    fseek(inFile,0,SEEK_SET);
 return fileSize;
}
void Wavtool_printinfo(char *filepath)
{
    wav_hdr wavHeader;
    FILE *wavFile;
    int headerSize = sizeof(wav_hdr),filelength = 0;
    wavFile = fopen(filepath,"r");
    if(wavFile == NULL)
    {
        printf("Can not able to open wave file\n");
        exit(EXIT_FAILURE);
    }
    fread(&wavHeader,headerSize,1,wavFile);
    filelength = getFileSize(wavFile);
    fclose(wavFile);
    printf("File is %d bytes.\n",filelength);
    printf("RIFF header                           :%c%c%c%c\n",wavHeader.RIFF[0],wavHeader.RIFF[1],wavHeader.RIFF[2],wavHeader.RIFF[3]);
    printf("WAVE header                           :%c%c%c%c\n",wavHeader.WAVE[0],wavHeader.WAVE[1],wavHeader.WAVE[2],wavHeader.WAVE[3]);
    printf("FMT                                   :%c%c%c%c\n",wavHeader.fmt[0],wavHeader.fmt[1],wavHeader.fmt[2],wavHeader.fmt[3]);
    printf("Data size (based on bits used)        :%ld\n",wavHeader.ChunkSize);

    // Display the sampling Rate form the header
    printf("Sampling Rate                         :%ld\n",wavHeader.SamplesPerSec); //Sampling frequency of the wav file
    printf("Number of bits used                   :%d\n",wavHeader.bitsPerSample); //Number of bits used per sample
    printf("Number of channels                    :%d\n",wavHeader.bitsPerSample);     //Number of channels (mono=1/sterio=2)
    printf("Number of bytes per second            :%ld\n",wavHeader.bytesPerSec);   //Number of bytes per second
    printf("data length                           :%ld\n",wavHeader.Subchunk2Size);

}
wav_hdr* Wavtool_read(char *filepath)
{
	wav_hdr* wavHeader;
    FILE *wavFile;
    wavHeader = (wav_hdr*)malloc(sizeof(wav_hdr));
    int headerSize = sizeof(wav_hdr);
    wavFile = fopen(filepath,"r");
    if(wavFile == NULL)
    {
        printf("Can not able to open wave file\n");
        exit(EXIT_FAILURE);
    }
    fread(wavHeader,headerSize,1,wavFile);
    fclose(wavFile);
    return wavHeader;
}
/* find the file size */
wav_inf* Wavtool_readSam(char *filepath)
{
	wav_inf *re;
    int headerSize = sizeof(wav_hdr);
    re = (wav_inf*)malloc(sizeof(wav_inf));


    QFile file(filepath);
    if (!file.open(QIODevice::ReadOnly)) {
       return NULL;
    }
    QByteArray blob = file.readAll();


    //fread(&(re->header),headerSize,1,wavFile);
    memcpy(&(re->header),blob.data(),headerSize);
    re->samples=(short*)malloc(re->header.Subchunk2Size);
    memcpy(re->samples,blob.data()+headerSize,re->header.Subchunk2Size);
    re->nsample = re->header.Subchunk2Size/(re->header.bitsPerSample/8);
    return re;
}
