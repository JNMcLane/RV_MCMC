; Fit all RVs for a single target from a given observing run
;
; name = name of target
; savename = name of file for saving RVs
; nmc = number of monte carlo iterations to run (for error
;determination)
; plot = set to plot model fits as they are completed
;dec 8 15, added mcmc passthrough

pro fit_spec_all, name, savename=savename, $
                  nmc=nmc, plot=plot, mcmc=mcmc

  ; Set default for output filename
  if n_elements(savename) eq 0 then $
     savename = name+'_rvfit.sav'

  ; Read in parameter file to get filenames and number of observations
  fit_spec_readparfile, name+'_pars.dat', mbound, ids, par, dpar, $
                       skyfile, starfile
  n = n_elements(ids) ; Number of observations

  first = 1 ; Flag the first pass through the loop
  for i=0, n-1 do begin ; Loop to fit each observation

     ; Read header info from ech file
     rdech, sp, name+'_'+ids[i], /silent ; Reach ech file
     oname = sxpar(sp.head, 'OBJECT')    ; Object name
     ra = sxpar(sp.head, 'RA')           ; R.A.
     dec = sxpar(sp.head, 'DEC')         ; Declination
     mjd = sxpar(sp.head, 'MJD')         ; MJD of observation
     hcor = sxpar(sp.head, 'HELCOR')     ; Heliocentric RV correction
     beam = sxpar(sp.head, 'BEAM')       ; Telescope beam position (A/B)
     
     ; Fit spectrum, get RV
     fit_spec, name, ids[i], v, dv, p, dp, chi, plot=plot, $
               wavl = w, spec = s, unc = u, fit = fit, $
               nmc = nmc, mcmc=mcmc
     
     ; Initialize data structure first time through the loop
     if first then begin
        first = 0
        np = n_elements(p)      ; Number of model parameters
        npix = n_elements(w)    ; Number of spectrum elements
        f = { $
            name : oname, $        ; Target name
            ra : ra, $             ; RA
            dec : dec, $           ; Dec
            beam : strarr(n), $    ; Telescope beam
            id : lonarr(n), $      ; Observation ID
            t : dblarr(n), $       ; Time of observation (MJD)
            hcor : fltarr(n), $    ; Heliocentric RV correction
            v : fltarr(n), $       ; Heliocentric RV of target
            dv : fltarr(n), $      ; RV uncertainty
            p : fltarr(np,n), $    ; Best fit model parameters
            dp : fltarr(np,n), $   ; Parameter uncertainties 
            chi : fltarr(n), $     ; Chi^2 of fit
            snr : fltarr(n), $     ; Signal-to-noise of observation
            w : fltarr(npix, n), $ ; Wavelength solution
            s : fltarr(npix, n), $ ; Observed spectrum
            u : fltarr(npix, n), $ ; Spectrum uncertainty
            f : fltarr(npix, n), $ ; Model spectrum 
            sky_fn : skyfile, $    ; Path to telluric model
            star_fn : starfile $   ; Path to stellar model
            }
     endif
     
     ; Save outcome of model fit
     f.beam[i] = strtrim(beam, 2)
     f.id[i] = ids[i]
     f.t[i] = mjd
     f.hcor[i] = hcor
     f.v[i] = v+hcor            ; NOTE: RV includes heliocentric correction!
     f.dv[i] = dv
     f.p[*,i] = p
     f.dp[*,i] = dp
     f.chi[i] = chi
     f.snr[i] = median(sp.spec/sp.sig)
     f.w[*,i] = w
     f.s[*,i] = s
     f.u[*,i] = u
     f.f[*,i] = fit
  endfor
  
  ; Print summary info to terminal
  print, ''
  print, '---------------------------------'
  print, 'FIT_SPEC_ALL summary for '+name+':'
  print, ''
  print, 'Mean RV = '+strcompress(string(mean(f.v)), /remove_all)+' km/s'
  if n_elements(f.v) gt 1 then $
     print, 'Std. Dev. RV = '+strcompress(string(1000*stddev(f.v)), /remove_all)+' m/s'
  print, ''
  save, f, filename=savename ; Save results to file
  print, 'Output saved as '+savename
  print, '---------------------------------'
  print, ''
  
end
