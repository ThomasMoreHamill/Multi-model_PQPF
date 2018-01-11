import numpy as np

gamma_threshes = \
    [0.0, 0.05, 0.1, 0.2, 0.3,           0.4, 0.5, 0.6, 0.7, 0.8, \
    0.9, 1.0, 1.2, 1.4, 1.6,            1.8, 2.0, 2.3, 2.6, 3.0, \
    3.3, 3.6, 4.0, 4.5, 5.0,            5.5, 6.0, 6.5, 7.0, 7.5, \
    8.0, 8.5, 9.0, 10.0, 11.0,          12.0, 13.0, 14.0, 15.0, 16.0, \
	18.0, 20.0, 22.5, 25.0, 27.5, 	    30.0, 33.0, 36.0, 40.0, 45.0, \
	50.0, 55.0, 60.0, 65.0, 70.0, 	    75.0, 80.0, 85.0, 90.0, 95.0, \
    100.0, 120.0, 140.0, 160.0, 180.0,  200.0, 250.0, 300.0 ]
nthreshes = len(gamma_threshes)
ilower_bound = np.zeros((nthreshes), dtype= np.int32)
iupper_bound = np.zeros((nthreshes), dtype= np.int32)
weightfn = np.zeros((nthreshes, nthreshes), dtype=np.float32)
gamma_threshes_alog = np.zeros((nthreshes), dtype= np.float32)


gamma_threshes_alog[0] = -10.0
gamma_threshes_alog[1:nthreshes-1] = np.log(gamma_threshes[1:nthreshes-1])
print 'gamma_threshes_alog = ', gamma_threshes_alog 


# ---- process each index

weightfn[0,0] = 1.0  # effectively no weight for 0 on other positive

print 'ithresh, ilower_bound, iupper_bound, gamma_threshes[ithresh] = ', \
    0, ilower_bound[0], iupper_bound[0], gamma_threshes[0]
print '  weightfn ==', weightfn[0,ilower_bound[0]:iupper_bound[0]+1]

for ithresh in range(1, nthreshes):
    if gamma_threshes[ithresh] < 1.:
        gmin = 0.0001
        gmax = 3.0
    else:
        gmin = gamma_threshes[ithresh]/2.
        gmax = gamma_threshes[ithresh]*2.
    nofindit = True
    for it2 in range(nthreshes):
        if gamma_threshes[it2] >= gmin and nofindit:
            ilower_bound[ithresh] = it2
            nofindit = False

    nofindit = True
    for it2 in range(nthreshes-1,0,-1):
        if gamma_threshes[it2] <= gmax and nofindit:
            iupper_bound[ithresh] = it2
            nofindit = False
    
    #  --- define the weights for indices inside bounds. 
    
    for iw1 in range(ilower_bound[ithresh], iupper_bound[ithresh]+1):
        denom = np.max([2. - gamma_threshes[ithresh], 0.6931])
        weightfn[ithresh,iw1] = 1.0 - \
            np.abs(gamma_threshes_alog[iw1] - gamma_threshes_alog[ithresh]) / denom
        if weightfn[ithresh,iw1] < 0.0: weightfn[ithresh,iw1] = 0.0

    print 'ithresh, ilower_bound, iupper_bound, gamma_threshes[ithresh] = ', \
        ithresh, ilower_bound[ithresh], iupper_bound[ithresh], gamma_threshes[ithresh]
    print '  weightfn ==', weightfn[ithresh,ilower_bound[ithresh]:iupper_bound[ithresh]+1]