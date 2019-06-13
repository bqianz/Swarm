from __future__ import division # true divide of arrays
import numpy as np
import math
from numpy import linalg as LA

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
    print(generate_array(4,1,10))