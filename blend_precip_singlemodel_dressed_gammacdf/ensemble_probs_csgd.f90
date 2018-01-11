SUBROUTINE ensemble_probs_csgd(nxa, nya, nthreshes, ncsgd_params, pthreshes, &
    conusmask, ensmean, stddev, POP, CSGD_climo_mean, CSGD_climo_mu, &
    CSGD_climo_sigma, CSGD_climo_shift, csgd_parameters, rho, &
    prob_forecast_CSGD)
    
INTEGER, INTENT(IN) :: nxa, nya, nthreshes, ncsgd_params
INTEGER*2, INTENT(IN), DIMENSION(nxa, nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes
REAL, INTENT(IN), DIMENSION(nxa, nya) :: ensmean
REAL, INTENT(IN), DIMENSION(nxa, nya) :: stddev 
REAL, INTENT(IN), DIMENSION(nxa, nya) :: POP
REAL, INTENT(IN), DIMENSION(nxa, nya) :: CSGD_climo_mean
REAL, INTENT(IN), DIMENSION(nxa, nya) :: CSGD_climo_mu
REAL, INTENT(IN), DIMENSION(nxa, nya) :: CSGD_climo_sigma
REAL, INTENT(IN), DIMENSION(nxa, nya) :: CSGD_climo_shift
REAL, INTENT(IN), DIMENSION(nxa, nya) :: rho
!REAL, INTENT(IN), DIMENSION(nxa, nya) :: shift
REAL*8, INTENT(IN), DIMENSION(ncsgd_params) :: csgd_parameters
REAL, INTENT(OUT), DIMENSION(nxa, nya, nthreshes) :: prob_forecast_CSGD

REAL*8 logarg, ensmean_ano, mu, muratio, sigma, shape, scale
REAL*8 log1p, expm1, pgamma

PRINT *,'CSGD_climo_mean(1:nxa:10,nya/2) = ',CSGD_climo_mean(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_mu(1:nxa:10,nya/2) = ',CSGD_climo_mu(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_sigma(1:nxa:10,nya/2) = ',CSGD_climo_sigma(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_shift(1:nxa:10,nya/2) = ',CSGD_climo_shift(1:nxa:10,nya/2)
PRINT *,'rho(1:nxa:10,nya/2) = ',rho(1:nxa:10,nya/2)
PRINT *,'pthreshes = ', pthreshes
PRINT *,'conusmask(1:nxa:10,nya/2) = ', conusmask(1:nxa:10,nya/2)
PRINT *,'ensmean(1:nxa:20, nya/2) = ', ensmean(1:nxa:20, nya/2) 
PRINT *,'stddev(1:nxa:20, nya/2) = ', stddev(1:nxa:20, nya/2)
PRINT *,'POP(1:nxa:20, nya/2) = ', POP(1:nxa:20, nya/2)
PRINT *,'csgd_parameters = ', csgd_parameters

DO ithresh = 1, nthreshes
    rthresh = pthreshes(ithresh)
    DO jya = 1, nya
        DO ixa = 1, nxa
            IF (conusmask(ixa,jya) .eq. 1) THEN
                                
                ensmean_ano = ensmean(ixa,jya) / CSGD_climo_mean(ixa,jya)
                
                !OLD: logarg = par[1] + par[2]*enspop + par[3]*ensmean_ano
                !OLD: mu = mu_cl * np.log1p(np.expm1(par[0])*logarg) / par[0]
                !OLD: sigma = par[4]*sigma_cl*np.sqrt(mu/mu_cl) + par[5]*ensstdev   
                                           
                !logarg = csgd_parameters(2) + &
                !   csgd_parameters(3)*POP(ixa,jya) + &
                !    csgd_parameters(4)*ensmean_ano
                !mu = CSGD_climo_mu(ixa,jya) * &
                !    log1p(expm1(csgd_parameters(1))*logarg) / csgd_parameters(1)
                !muratio = mu / CSGD_climo_mu(ixa,jya)
                !sigma = csgd_parameters(5)*CSGD_climo_sigma(ixa,jya)* &
                !    SQRT(muratio) + csgd_parameters(6)*stddev(ixa,jya)
                              
                !NEW:  logarg = np.exp(-rho/par[1]) + par[2]*enspop + par[3]*rho*ensmean_ano 
                !NEW:  mu = mu_cl * np.log1p(np.expm1(par[0])*logarg) / par[0]
                !NEW:  sigma = par[4]*sigma_cl*np.sqrt(1.-np.square(rho))*np.sqrt(mu/mu_cl) + par[5]*rho*ensstdev 
                
                logarg = exp(-rho(ixa,jya)/csgd_parameters(2)) + &
                    csgd_parameters(3)*POP(ixa,jya) + &
                    csgd_parameters(4)*rho(ixa,jya)*ensmean_ano             
                mu = CSGD_climo_mu(ixa,jya) * &
                    log1p(expm1(csgd_parameters(1))*logarg) / csgd_parameters(1)
                muratio = mu / CSGD_climo_mu(ixa,jya)
                omrho2 = 1.-rho(ixa,jya)**2
                sigma = csgd_parameters(5)*CSGD_climo_sigma(ixa,jya)* &
                    SQRT(omrho2)*SQRT(muratio) + csgd_parameters(6)*rho(ixa,jya)*stddev(ixa,jya)             
                                     
                !shift  = shift_cl
                !shape  = (mu/sigma)**2
                !scale  = mu/shape
                !prob_pop(ix,jy)           = REAL(pgamma( 0.25d0-shift, shape, scale, 0, 0))
                !prob_1mm(ix,jy)           = REAL(pgamma( 1.0d0-shift, shape, scale, 0, 0))
                            
                shape = (mu/sigma)**2
                scale = mu/shape

                prob_forecast_CSGD(ixa,jya,ithresh) = &
                    REAL(pgamma( DBLE(rthresh)-CSGD_climo_shift(ixa,jya), shape, scale, 0, 0))
            ELSE
                prob_forecast_CSGD(ixa,jya,ithresh) = -99.99
            ENDIF
        END DO ! ixa
    END DO ! jya

    PRINT *,'rthresh = ', rthresh
    PRINT *,'prob_forecast_CSGD(1:nxa:10, nya/2, ithresh) = ', &
        prob_forecast_CSGD(1:nxa:10, nya/2, ithresh) 
END DO ! ithresh


RETURN
END SUBROUTINE ensemble_probs_csgd

