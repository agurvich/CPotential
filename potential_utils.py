
import os
import numpy as np 
import ctypes
import copy
import time
from memory_profiler import profile,memory_usage
import itertools
import warnings

import pdb

def enter(sleep_time=1):
    init_time = time.time()
    while time.time() - init_time < sleep_time:
        pass
    pass

def exit(sleep_time=1):
    init_time = time.time()
    while time.time() - init_time < sleep_time:
        pass
    pass
    
def subCalculateCPairwiseDistancesSquared(
    test_xs,test_ys,test_zs,
    xs,ys,zs,
    fin_array=None,
    count=0,tot_count=0):

    ## what are the lengths of the test and full arrays?
    Ntest = test_xs.shape[0]
    Narr = xs.shape[0]

    ## how many pairwise distances will we be finding?
    ##  use partial arithmetic series sum formula 
    ##  n/2 (a1 + an) 
    arr_size = int(Ntest/2.0 * (Narr - 1 - count + Narr - count - Ntest))

    ## do we have a place to put the distances?
    if fin_array is None:
        fin_array = np.zeros(arr_size,dtype=np.float32)

    ## if we're sub-stepping through the list, don't need to calculate distances
    ##  for points we've already done
    new_xs = xs[count:]
    new_ys = ys[count:]
    new_zs = zs[count:]

    print(Ntest,'this run -> ',arr_size,'many distances')

    ## call the c executable
    exec_call = os.path.join(os.environ['HOME'],"python/CPotential/src/c_potential.so")
    c_obj = ctypes.CDLL(exec_call)

    print("Calling the C executable...",)
    c_obj.calcTestPairwiseDistsSquared(
        ctypes.c_int(Ntest),
        test_xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        test_ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        test_zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        ctypes.c_int(Narr-count),
        new_xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        new_ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        new_zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        ctypes.c_int(tot_count),
        fin_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
    print("... done!",)

    return arr_size,fin_array

def calculateCPairwiseDistancesSquared(all_pos,ncutoff=2e4):
    ncutoff = int(ncutoff)

    ## recast and copy to safe memory
    all_pos = all_pos.astype('f')
    xs,ys,zs = all_pos.T

    Narr = all_pos.shape[0]

    ## arithmetic series sum formula
    ##  sum of arithmetic series = n*(n+1)/2 with n = Narr -1
    big_dist_count = (Narr-1)*(Narr+1-1)//2

    if Narr > ncutoff:
        nruns = int(Narr/ncutoff) + ((Narr%ncutoff)!=0)

        fin_array = np.zeros(big_dist_count,dtype=np.float32)

        ## accumulation variables
        count = 0
        tot_count = 0

        for run in range(nruns):

            test_xs = xs[run*ncutoff:(run+1)*ncutoff]
            test_ys = ys[run*ncutoff:(run+1)*ncutoff]
            test_zs = zs[run*ncutoff:(run+1)*ncutoff]

            Ntest = test_xs.shape[0]

            print("Working on subset",run+1,"out of",nruns)
            arr_size,fin_array = subCalculateCPairwiseDistancesSquared(
            test_xs,test_ys,test_zs,
            xs,ys,zs,
            fin_array,
            count,tot_count)

            count += Ntest
            tot_count += arr_size

        ## make sure we calculated a distance to all the points
        assert tot_count==big_dist_count

        return fin_array
    else:
        xs = copy.copy(xs)
        ys = copy.copy(ys)
        zs = copy.copy(zs)

        ## prepare the output pointer
        ##  sum of arithmetic series = n*(n+1)/2 with n = Narr -1
        h_out_cast=ctypes.c_float*(big_dist_count)
        H_OUT=h_out_cast()

        ## call the c executable
        exec_call = os.path.join(os.environ['HOME'],"python/CPotential/src/c_potential.so")
        c_obj = ctypes.CDLL(exec_call)

        print("Calling the C executable...",)
        c_obj.calcPairwiseDistsSquared(
            ctypes.c_int(Narr),
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

            H_OUT)
        print("... done!",)

        h=np.ctypeslib.as_array(H_OUT)
        print(h.shape[0],'many distances')
        return h

def calculateCPotential(all_pos,all_masses,test_pos,print_flag = 1,nthreads =1):
    ## recast 
    all_pos = all_pos.astype('f')
    all_masses = all_masses.astype('f')

    test_pos = test_pos.astype('f')

    ## unpack positions
    xs,ys,zs = all_pos.T
     ## need to reallocate arrays in memory 
    xs = np.ascontiguousarray(xs)#copy.copy(xs)
    ys = np.ascontiguousarray(ys)#copy.copy(ys)
    zs = np.ascontiguousarray(zs)#copy.copy(zs)

    test_xs,test_ys,test_zs = test_pos.T
    ## need to reallocate arrays in memory 
    test_xs = np.ascontiguousarray(test_xs)#copy.copy(test_xs)
    test_ys = np.ascontiguousarray(test_ys)#copy.copy(test_ys)
    test_zs = np.ascontiguousarray(test_zs)#copy.copy(test_zs)

    Narr = all_pos.shape[0]
    Ntest = test_pos.shape[0]

    ## prepare the output pointer
    h_out_cast=ctypes.c_float*Ntest
    H_OUT=h_out_cast()

    ## call the c executable

    ## alternatively can link to this __file__'s location
    exec_call = os.path.join(os.environ['HOME'],"python/CPotential/src/c_potential.so")
    c_obj = ctypes.CDLL(exec_call)


    ## from "extensive" testing this seems to be around the point where you get
    ##  consistent speedup. interestingly increasing either Ntest or Narr will
    ##  do it (even though only the Narr portion is actually multi-threaded). 
    if np.log10(Narr*Ntest) < 8:
        warnings.warning("Not enough points for multi-threading to be useful, setting threads to 1")
        nthreads = 1

    if nthreads >1:
        if print_flag:
            print("Calling the Multi-C executable...",)
        c_obj.multiCalculateZGravityAtLocations(
            ctypes.c_int(nthreads),
            ctypes.c_int(Narr),
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

            all_masses.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), #masses

            ctypes.c_int(Ntest),
            test_xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            test_ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            test_zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

            ctypes.byref(H_OUT))
    else:
        if print_flag:
            print("Calling the C executable...",)
        c_obj.calculateZGravityAtLocations(
            ctypes.c_int(Narr),
            xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

            all_masses.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), #masses

            ctypes.c_int(Ntest),
            test_xs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            test_ys.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            test_zs.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),

            ctypes.byref(H_OUT))
    if print_flag:
        print("...done!")

    h=np.ctypeslib.as_array(H_OUT)
    return h


def test_pot(all_pos,all_masses,test_pos,nthreads =1):
    ## py-calculate
    from scipy.spatial.distance import cdist as cdist
    distss = cdist(test_pos,all_pos)

    rvects = all_pos-test_pos[:,None]
    zcomponents = (rvects/distss[:,:,None])[:,:,2]

    ALLTOGETHER=1.407e-7
    py_grav_cgs = np.sum(all_masses/distss**2*zcomponents,axis=1)*ALLTOGETHER

    ## c-calculate
    pos_cgs_zgravs = calculateCPotential(all_pos,all_masses,test_pos,print_flag = True,nthreads=nthreads)

    rats = py_grav_cgs / pos_cgs_zgravs

    print('Test passed:',np.all((rats-1)**2 < 1e8))

    return rats

def test_pair(all_pos):
    print(calculateCPairwiseDistancesSquared(all_pos))

if __name__ == '__main__':


    nmasses = int(5e5)
    ntest = 50*4
    nthreads = 10
    dec_pos = np.random.random(size=(nmasses,3))
    dec_mass = np.ones(nmasses)
    test_pos = np.random.random(size=(ntest,3))
    init_time = time.time()
    print('%.1e elements'%nmasses)
    #test_pot(dec_pos,dec_mass,test_pos,nthreads = nthreads)
    calculateCPotential(dec_pos,dec_mass,test_pos,print_flag = True,nthreads=nthreads)
    multi_dur = time.time()-init_time
    calculateCPotential(dec_pos,dec_mass,test_pos,print_flag = True,nthreads=1)
    #test_pot(dec_pos,dec_mass,test_pos,nthreads = 1)
    single_dur = time.time()-init_time - multi_dur
    print('1 thread: %.4f s'%single_dur,'%d threads: %.4f s'%(nthreads,multi_dur))
    print(single_dur/multi_dur,'X speedup')

    """
    all_pos = np.array([[1,1,1],[2,2,2],[1,1,1],[2,2,2]])
    test_pair(all_pos)

    coords = np.random.random(size=(int(5e4),3))

    calculateCPairwiseDistancesSquared(coords,ncutoff=1e4)
    print(time.time()-init_time,'s elapsed calculating pairwise distances')
    """
    

