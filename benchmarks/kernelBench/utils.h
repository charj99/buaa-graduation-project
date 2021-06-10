#ifndef __UTILS_H
#define __UTILS_H

#define MAXN 256
#include <stdio.h>
#include <time.h>

typedef struct timespec time_measure;

void getTime(time_measure* t) {
    clock_gettime(CLOCK_MONOTONIC, t);
}

double measuring_difftime(struct timespec t0, struct timespec t1) {
  double secdiff = difftime(t1.tv_sec, t0.tv_sec);
  if (t1.tv_nsec < t0.tv_nsec) {
    long val = 1000000000l - t0.tv_nsec + t1.tv_nsec;
    secdiff += (double)val / 1e9 - 1.;
  } else {
    long val = t1.tv_nsec - t0.tv_nsec;
    secdiff += (double)val / 1e9;
  }
  return secdiff;
}

void print2file_2D(char* fileName, double (*u)[MAXN]) {
    FILE* outFile = fopen(fileName, "w");
    if (outFile == NULL) {
        printf("Can not open output file %s\n", fileName);
        return;
    }
    for (int i = 0; i < MAXN; i++) {
        for (int j = 0; j < MAXN; j++)
            fprintf(outFile, "%lf\t", u[i][j]);
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}

void print2file_3D(char* fileName, double (*u)[MAXN][MAXN]) {
    FILE* outFile = fopen(fileName, "w");
    if (outFile == NULL) {
        printf("Can not open output file %s\n", fileName);
        return;
    }
    for (int i = 0; i < MAXN; i++) {
        for (int j = 0; j < MAXN; j++) {
            for (int k = 0; k < MAXN; k++)
                fprintf(outFile, "%lf\t", u[i][j][k]);
            fprintf(outFile, "\n");    
        }
        fprintf(outFile, "\n\n\n");
    }
    fclose(outFile);
}
#endif // __UTILS_H
