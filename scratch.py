
# coding: utf-8

# In[27]:

import os
import numpy as np 
import time
from potential_utils import calculateCPotential


# In[28]:



np.random.seed(516)

## set-up some random positions
Narr=10**7
x,y,z = ((np.random.rand(3,Narr)-0.5)*15).astype('f')
masses = np.ones(Narr,dtype='f')

Ntest=1000
test_x,test_y,test_z = ((np.random.rand(3,Ntest)-0.5)*15).astype('f')

init_time = time.time()
h = calculateCPotential(np.array([x,y,z]).T,masses,np.array([test_x,test_y,test_z]).T)
print time.time()-init_time,'s elapsed'
