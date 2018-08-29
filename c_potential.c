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
        float dx = xs[i]-point[0];
        float dy = ys[i]-point[1];
        float dz = zs[i]-point[2];

        float dr = sqrt(dx*dx + dy*dy + dz*dz); // RIP speed
        summ += masses[i]/(dr*dr*dr)*dz; 
    }
    return ALLTOGETHER*summ;//cgs
}

float calculateGravityR(float * point, float * masses, float * xs, float * ys, float * zs,int Narr){
    float summ = 0;
    for (int i =0; i<Narr; i++){
        float dx = xs[i]-point[0];
        float dy = ys[i]-point[1];
        float dz = zs[i]-point[2];

        float dr2 = dx*dx + dy*dy + dz*dz;
        summ += masses[i]/dr2; 
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

int calcPairwiseDistsSquared(
    int Narr,
    float * xs, float * ys, float * zs,
    float * H_OUT ){

    float x,y,z;
    float dx,dy,dz;
    int count = 0;
    // loop over each point in the array
    for (int i=0; i<Narr; i++){
        // pick out this point
        x = xs[i];
        y = ys[i];
        z = zs[i];
        // loop over the remaining points, starting with the next one
        for (int j=i+1; j<Narr; j++){
            // calculate the separation
            dx = x - xs[j];
            dy = y - ys[j];
            dz = z - zs[j];

            // put this square distance into the output array and move onto the next
            H_OUT[count]=dx*dx + dy*dy + dz*dz;
            count++;
        }
    }
}
