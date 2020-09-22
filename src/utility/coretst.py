import numpy as np
import time

size = 10000

start = time.time()
a = np.random.random_sample((size, size))
b = np.random.random_sample((size, size))
n = np.dot(a,b)

print(time.time()-start)