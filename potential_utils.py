
import os
import numpy as np 
import ctypes
import copy


def calculateCPairwiseDistancesSquared(all_pos):
    ## recast and copy to safe memory
    all_pos = all_pos.astype('f')
    xs,ys,zs = all_pos.T
    xs = copy.copy(xs)
    ys = copy.copy(ys)
    zs = copy.copy(zs)

    Narr = all_pos.shape[0]

    ## prepare the output pointer
    ##  sum of arithmetic series = n*(n+1)/2 with n = Narr -1
    h_out_cast=ctypes.c_float*( (Narr-1)*(Narr+1-1)/2 )
    H_OUT=h_out_cast()

    ## call the c executable
    exec_call = os.path.join(os.environ['HOME'],"python/CPotential/c_potential.so")
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
    return h

def calculateCPotential(all_pos,all_masses,test_pos):
    ## recast 
    all_pos = all_pos.astype('f')
    all_masses = all_masses.astype('f')

    test_pos = test_pos.astype('f')

    ## unpack positions
    xs,ys,zs = all_pos.T
     ## need to reallocate arrays in memory 
    xs = copy.copy(xs)
    ys = copy.copy(ys)
    zs = copy.copy(zs)

    test_xs,test_ys,test_zs = test_pos.T
     ## need to reallocate arrays in memory 
    test_xs = copy.copy(test_xs)
    test_ys = copy.copy(test_ys)
    test_zs = copy.copy(test_zs)

    Narr = all_pos.shape[0]
    Ntest = test_pos.shape[0]

    ## prepare the output pointer
    h_out_cast=ctypes.c_float*Ntest
    H_OUT=h_out_cast()

    ## call the c executable

    ## alternatively can link to this __file__'s location
    exec_call = os.path.join(os.environ['HOME'],"python/CPotential/c_potential.so")
    c_obj = ctypes.CDLL(exec_call)

    print("Calling the C executable...",)
    c_obj.calcDists(
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
    print("...done!")

    h=np.ctypeslib.as_array(H_OUT)
    return h


def test_pot(all_pos,all_masses,test_pos):
    ## py-calculate
    from scipy.spatial.distance import cdist as cdist
    distss = cdist(test_pos,all_pos)

    rvects = all_pos-test_pos[:,None]
    zcomponents = (rvects/distss[:,:,None])[:,:,2]

    ALLTOGETHER=1.407e-7
    py_grav_cgs = np.sum(all_masses/distss**2*zcomponents,axis=1)*ALLTOGETHER

    ## c-calculate
    pos_cgs_zgravs = calculateCPotential(all_pos,all_masses,test_pos)

    rats = py_grav_cgs / pos_cgs_zgravs

    print(rats)

    return rats

def test_pair(all_pos):
    print(calculateCPairwiseDistancesSquared(all_pos))

if __name__ == '__main__':

    dec_pos = np.array([[0,0,0]]*5)
    dec_mass = np.array([1,1,1,1,1])
    test_pos = np.array([[0,0,1]]*3)
    #print(test_pos)
    #test_pot(dec_pos,dec_mass,test_pos)

    all_pos = np.array([[1,1,1],[2,2,2],[1,1,1],[2,2,2]])
    test_pair(all_pos)

