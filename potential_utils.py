
import os
import numpy as np 
import ctypes

def calculateCPotential(all_pos,all_masses,test_pos):
    ## recast 
    all_pos = all_pos.astype('f')
    all_masses = all_masses.astype('f')

    test_pos = test_pos.astype('f')

    ## unpack positions
    xs,ys,zs = all_pos.T
    test_xs,test_ys,test_zs = test_pos.T

    Narr = all_pos.shape[0]
    Ntest = test_pos.shape[0]

    ## prepare the output pointer
    h_out_cast=ctypes.c_float*Ntest
    H_OUT=h_out_cast()

    ## call the c executable

    ## alternatively can link to this __file__'s location
    exec_call = os.path.join(os.environ['HOME'],"python/CPotential/c_potential.so")
    c_obj = ctypes.CDLL(exec_call)

    print "Calling the C executable...",
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
    print "Finished!"

    h=np.ctypeslib.as_array(H_OUT)
    return h

