def calculate_quantile_for_precip_values(npthreshes, nya, nxa, \
    rthreshes, fraction_zero, alphahat, betahat, conusmask_in):

    """ 
    for input precip amount thresholds in rthreshes array, and given 
    fraction zero and gamma distribution alpha and beta parameters,
    determine the cumulative non-exceedance probability associated
    with this input threshold.
    
    coded by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060
    
    compile with python setup_calculate_quantile_for_precip_values.py build_ext --inplace
    
    """
    
    import numpy as np
    import scipy.stats as stats
    
    quant_fpamt = np.zeros((npthreshes,nya,nxa), dtype=np.float64)
    for ixa in range(nxa):
        if ixa%20 == 0: print 'processing ',ixa,' of ', nxa
        for jya in range(nya):
            if conusmask_in[jya,ixa] == 1 and alphahat[jya,ixa] > 0.0 :
                for ithresh in range(npthreshes):
                    if fraction_zero[jya,ixa] < 1.0:
                        quant_fpamt[ithresh,jya,ixa] = fraction_zero[jya,ixa] + \
                            (1.0 - fraction_zero[jya,ixa])*\
                            stats.gamma.cdf(rthreshes[ithresh],\
                            alphahat[jya,ixa], loc=0, scale=betahat[jya,ixa])
                    else:
                        quant_fpamt[ithresh,jya,ixa] = 1.0
            else:
                quant_fpamt[:,jya,ixa] = -99.99
    
    return quant_fpamt

