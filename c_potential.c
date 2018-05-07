#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define GCGS 6.674e-8
#define KPCTOCM 3.08e21
#define CODETOGRAMS 2e43

#define ALLTOGETHER 1.407e-7 // GCGS*CODETOGRAMS/(KPCTOCM*KPCTOCM)

void printArray(float * arr,int Narr){
    for (int i=0; i< Narr; i++){
        printf("%.2f\t",arr[i]);
    }
    printf("\n");
}

float calculateGravityZ(float * point, float * masses, float * xs, float * ys, float * zs,int Narr){
    float summ = 0;
    for (int i =0; i<Narr; i++){
        float dx = point[0]-xs[i];
        float dy = point[1]-ys[i];
        float dz = point[2]-zs[i];

        float dr = sqrt(dx*dx + dy*dy + dz*dz); // RIP speed
        summ += masses[i]/(dr*dr*dr)*dz; 
    }
    return ALLTOGETHER*summ;//cgs
}

int calcDists(
    int Narr,
    float * xs, float * ys, float * zs,
    float * masses,
    int Ntest,
    float * test_xs, float * test_ys, float * test_zs,
    float * H_OUT ){

    float point[3];
    for (int i=0; i<Ntest; i++){

        point[0]=test_xs[i]; point[1]=test_ys[i]; point[2]=test_zs[i];

        H_OUT[i]=calculateGravityZ(point,masses,xs,ys,zs,Narr);

    }

    return 1;
}


