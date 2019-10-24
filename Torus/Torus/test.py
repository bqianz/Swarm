# from __future__ import division # true divide of arrays

import numpy as np
import math
from numpy import linalg as LA

x = 4
y = 1


# case of not crossing half theta
def general_segment(a,b,n):
    half = int(n/2)
    s = np.sign(b-a)
    temp = np.arange(a,b,s)
    if temp.size > half:
        return np.arange(a,b-s*n,-s) % n
    else:
        return temp % n

if __name__ == "__main__":
    print(general_segment(x,y,10))
    print(x+y)
    arr1 = np.ones(7)
    arr2 = np.ones(7) * 2
    print(arr1/arr2)
    print(arr1 * arr2)