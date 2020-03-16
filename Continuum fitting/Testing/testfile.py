import numpy as np
#percent selction algorithm

array = np.array([7,3,4,2,1,5,6,10,8,9])
pct = 0.3

def findpct(array,pct):
    sortedarray = np.argsort(array)
    selection = int(len(array)*pct)
    selectedvals = sortedarray[-selection:]
    return selectedvals
