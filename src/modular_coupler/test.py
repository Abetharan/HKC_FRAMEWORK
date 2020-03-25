import numpy as np

from scipy.signal import argrelextrema
# for local maxima
argrelextrema(xNorm, np.greater)

# for local minima
argrelextrema(xNorm, np.less)
