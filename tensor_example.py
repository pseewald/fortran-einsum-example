#!/usr/bin/env python
from __future__ import print_function

import numpy as np

a = np.arange(60.).reshape(3,4,5)
b = np.arange(24.).reshape(4,3,2)
c = np.einsum(a, [1,2,3], b, [2,1,4], [3,4])

for index, val in np.ndenumerate(c):
    print(index, val)

a = np.arange(10.).reshape(5,2)
b = np.arange(30.).reshape(3,5,2)
c = np.einsum(a,[1,2], b, [3,1,4], [4,2,3])

for index, val in np.ndenumerate(c):
    print(index, val)

a = np.arange(24.).reshape(3,2,4)
b = np.arange(8.).reshape(4,2)
c = np.einsum(a,[1,2,3], b, [3,2], [1])

for index, val in np.ndenumerate(c):
    print(index, val)

a = np.arange(6.).reshape(3,2)
b = np.arange(5.).reshape(5)
c = np.einsum(a, [1,2], b, [3], [1,3,2])

for index, val in np.ndenumerate(c):
    print(index, val)

a=np.arange(6.).reshape(3,2)
b = np.arange(6.).reshape(2,3)
c=np.einsum(a, [1,2], b, [2,1])

for index, val in np.ndenumerate(c):
    print(index, val)
