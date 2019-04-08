#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "include.h"

int calcTestPairwiseDistsSquared(
    int Ntest,
    float * test_xs, float * test_ys, float * test_zs,
    int Narr,
    float * xs, float * ys, float * zs,
    int count_skip,
    float * H_OUT ){

    float x,y,z;
    float dx,dy,dz;
    int count = 0;
    // loop over each point in the array
    for (int i=0; i<Ntest; i++){
        // pick out this point
        // if test_xs is xs[:Ntest] then you can iteratively call this function 
        //  externally to save memory
        x = test_xs[i];
        y = test_ys[i];
        z = test_zs[i];
        // loop over the remaining points, starting with the next one
        // NOTE the fact that j = i+1 to start **assumes** that test_xs = xs[:Ntest]
        for (int j=i+1; j<Narr; j++){
            // calculate the separation
            dx = x - xs[j];
            dy = y - ys[j];
            dz = z - zs[j];
            // put this square distance into the output array and move onto the next
            H_OUT[count+count_skip]=dx*dx + dy*dy + dz*dz;
            count++;
        }
    }
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
