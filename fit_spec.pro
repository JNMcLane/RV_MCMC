;
; This is a driver program to fit an IRTF CSHELL spectrum to
; determine its radial velocity with respect to the telluric lines in
; the spectrum.
;
; INPUTS:
;   objname = Name of target
;   obsid   = Observation ID
;
; OUTPUTS:
;   vel     = The radial velocity of the spectrum.
;   dvel    = The uncertainty in the velocity estimated by the fitting 
;              software.
;   parfit  = The best fit model parameters
;   sigpar  = The uncertainty in the model parameters
;   chi     = The reduced chi-square value for the best fit model
;
; OPTIONAL INPUTS:
;   nmc     = Number of Monte Carlo iterations to run (for error
;              estimation)
;   ps      = Name of postscript file to save plot
;
; KEYWORDS:
;   plot    = Plots the final spectrum and fit if set.
;   init    = Set up everything for fitting, but don't fit - useful for
;              checking initial guess.
;
; OPTIONAL OUTPUTS:
;   wavl    = Wavelength solution
;   spec    = Observed spectrum
;   unc     = Spectrum uncertainty
;   fit     = Model spectrum
;
; HISTORY:
;   20-March-2008 Written by CMJ.
;   08-December-2015  -Added MCMC fitting
;                     -Changed Chi square fitting, 
;                     -Added ability to auto update param files (needs improvement)
;                     JNM

function fs_chi2, x, y, wts, p, dp
; FitSpec Chi2 - Determine reduce chi^2 of model
; INPUTS:
;  x   = Independent variable
;  y   = Dependent variable
;  wts = Weights on dependent variable
;  p   = Model parameters
;  dp  = Parameter step values (0 = fixed)
;
; RETURN VALUE
;  Reduced chi^2 of model
;

  ifree = where(dp ne 0,nfree)  ; Number of free parameters
  iwhr = where(wts ne 0, ncon)  ; Number of free dependent variables

  degf = ncon - nfree           ; Number of degrees of freedom
  if degf lt 0 then message,'Too few data points - not enough constraints.'
  degf = degf > 1               ;no div by zero, if degf=0

  fit = spec_model(x, p)               ; Calculate model
 ; chi2 = total(wts*(y - fit)^2.0)/degf ; Reduced chi-square
  chi2 = total(((y-fit)^2.0)/y)

  return, chi2
end


pro fs_fit, x, s, wts, pinit, dp, pfit, wfit, fit, chi, markov = markov
; FitSpec Fitter
; INPUTS:
; x     = array of pixel numbers
; s     = spectrum
; wts   = spectrum weights
; pinit = initial guess for model parameters
; dp    = parameter step values
; 
; OUTPUTS:
; pfit = Best fit model parameters
; wfit = Best fit wavelength solution
; fit  = Model spectrum
; chi  = Reduced chi^2 of model
;

  npar = n_elements(pinit)      ; Number of model parameters

  ; Initialize parinfo structure for MPFIT
  parinfo = replicate({value:0d, step:0d, fixed:0, limited:[0,0], $
                       limits:[0d,0d]}, npar)
  if keyword_set(markov) then begin
    parinfo[0].fixed=1
    parinfo[1].fixed=1
    parinfo[2].fixed=1
    parinfo[3].fixed=1
    parinfo[4].fixed=1
    parinfo[*].step = dp
    parinfo[*].value = pinit
  
  endif else begin
    parinfo[1].limited[0] = 1
    parinfo[1].limits[0] = 0.1
    parinfo[3].limited[0]=1
    parinfo[3].limits[0]=0.1
    parinfo[2].fixed = 1          ; Fix telluric velocity shift
    parinfo[4].fixed = 1          ; Fix vsini
    parinfo[*].step = dp          ; Set parameter step size
    parinfo[*].value = pinit      ; Set initial values
  endelse

  ; Least-squares fitting of model
  pfit = mpfitfun('SPEC_MODEL', x, s, weights = wts, $
                  yfit=fit, parinfo=parinfo, /quiet)

  wfit = poly(x, pfit[9:*])          ; Wavelength solution
  chi = fs_chi2(x, s, wts, pfit, dp) ; Reduced chi^2
end


pro accept, mcpar, newpar, x, s, wts, dpar, parout
  chi_1=fs_chi2(x, s, wts, mcpar, dpar)
  chi_2=fs_chi2(x, s, wts, newpar, dpar)
  ratio=chi_1/chi_2
  if ratio ge 1 then parout=newpar
  if ratio lt 1 then begin
    pass=randomu(seed)
    if ratio le pass then parout=newpar else parout=mcpar
  endif
end


pro gibbs, val, min, max, step, v_new
  v_new = val
  for itt=1, 100 do begin
    v_new = v_new - step + 2*step*randomu(seed)
    if v_new le min then v_new = min
    if v_new ge max then v_new = max
  endfor
end


pro mcmc, parfit, x, s, wts, dpar, bestpar
  mcpar=parfit
  bestpar=parfit
  chibest=fs_chi2(x, s, wts, bestpar, dpar)
  print,'Starting MCMC'
  for j=1, 1000 do begin
    newpar=mcpar
    gibbs, newpar[0], -100, 100, 1, v_new
    newpar[0] = v_new
    accept,mcpar, newpar, x, s, wts, dpar, parout
    mcpar=parout
    chimc = fs_chi2(x, s, wts, mcpar, dpar)
    if chimc lt chibest then begin
      chibest = chimc
      bestpar = mcpar
    endif
    newpar=mcpar
    gibbs,newpar[1], 0.1, 2, 0.01, v_new
    newpar[1] = v_new
    accept,mcpar, newpar, x, s, wts, dpar, parout
    mcpar=parout
    chimc = fs_chi2(x, s, wts, mcpar, dpar)
    if chimc lt chibest then begin
      chibest = chimc
      bestpar = mcpar
    endif
    newpar=mcpar
    gibbs,newpar[3], 0.1, 2, 0.01, v_new
    newpar[3] = v_new
    accept,mcpar, newpar, x, s, wts, dpar, parout
    mcpar=parout
    chimc = fs_chi2(x, s, wts, mcpar, dpar)
    if chimc lt chibest then begin
      chibest = chimc
      bestpar = mcpar
    endif
  endfor   
  print,'Finished MCMC'
end


; ------- MAIN PROGRAM -----------
pro fit_spec, objname, obsid, $                 ; Inputs
              vel, dvel, parfit, sigpar, chi, $ ; Outputs
              nmc = nmc, ps=ps, $               ; Optional Inputs
              mcmc = mcmc, $                    ; MCMC input, optional
              plot = plot, init = init, $       ; Keywords
              wavl = wout, spec = sout, $       ; Optional Outputs
              unc = uout, fit = fit, write = write

 common share2, skymodel, starmodel ; Common block for passing around input models

  ; Check usage
  if n_params() lt 2 then begin
     print,'fit_spec, name, id, v, dv, par, sigpar, chi, /plot, /init'
     retall
  endif

  print, '--------------------------------------'

  ; Set default number of Monte Carlo iterations
  if n_elements(nmc) eq 0 then nmc = 100

  ; Read parameter file.
  fit_spec_readparfile, objname+'_pars.dat', mask_boundaries, ids, allpar, dpar, $
                        skyfile, starfile

  ; Isolate model parameters for requested night:
  k = where(strcmp(ids, obsid), kcount) ; Find obsid in parameter file
  if kcount eq 0 then begin             ; Can't find it, return
     print, 'Invalid observation ID: '+string(obsid)
     retall
  endif
  parinit = allpar[*,k]             ; Extract parameter guesses for requested obsid
  sigpar = parinit*0.0              ; Initialize parameter uncertainty storage
  npar = n_elements(parinit)        ; Number of model parameters

  ; Get observed spectrum:
  fname = objname+'_'+strcompress(string(obsid), /remove_all) ; Filename of observation
  print,'Restoring '+fname
  print, ''
  rdech, obj, fname, /silent    ; Read ECH format

  s = obj.spec                  ; Observed spectrum
  u = obj.sig                   ; Data uncertainty

  ; Apply mask
  ; Fit procedure ignores datapoints where uncertainty = 0
  ; Use mask boundaries in par file to set appropriate unc to zero
  uout = u                                                      ; Save uncertainty for output before masking
  wts = 1./u^2                                                  ; Weights for MPFIT
  nmaskpairs = n_elements(mask_boundaries)/2                    ; Number of mask regions
  mask = intarr(n_elements(s))+1                                ; Initialize mask array
  if mask_boundaries[0] ne mask_boundaries[1] then $            ; If first two mask boundaries are equal, 
                                                                ;  then skip
     for i=0, nmaskpairs-1 do $                                 ; Step through each mask region
        mask[mask_boundaries[i*2]:mask_boundaries[i*2 + 1]] = 0 ; Set mask to zero in region
  wts = wts*mask                                                ; Apply mask to weights

  ; Read in model templates:
  print, 'Using following templates:'
  print, '  Telluric = '+skyfile
  print, '  Photosphere = '+starfile
  print, ''
  if ~file_test(skyfile) then begin ; Check for existence of telluric model
     print, 'Telluric template: '+skyfile+'.  Does not exist!'
     retall
  endif
  if ~file_test(starfile) then begin ; Check for existence of stellar model
     print, 'Photosphere template: '+starfile+'.  Does not exist!'
     retall
  endif
  ; Read telluric model
  restore, skyfile
  skymodel = {w : watm, s : satm}
  ; Restore stellar model, force all values to be positive
  restore, starfile
  sm.s = sm.s > 0
  starmodel = sm

  ; Initialize pixel number array
  npix = n_elements(s)          ; Number of spectrum elements
  x = dindgen(npix) - npix/2    ; Pixel count, set zero to middle of spectrum

  ; Start the model fitting!
  print, 'Fitting the model...'
  if keyword_set(mcmc) then begin

    if keyword_set(init) then begin          ; Init only, generate model with parameter guesses
     fit = spec_model(x,parinit)           ; Create model spectrum
     wfit = poly(x, parinit[9:*])          ; Create wavelength solution
     ;print, wfit
     parfit = parinit                      ; Set best fit values equal to guesses
     chi = fs_chi2(x, s, wts, parfit, dpar)  ; Chi^2 of "fit"

    endif else begin
      ; Create nmc randomized spectra for Monte Carlo
      b = rotate(replicate(1, nmc), 1)
      smc = s # b + randomn(seed, n_elements(s), nmc) * (u # b)

      ; Initialize parameter storage for all Monte Carlo runs
      parfit_mc = fltarr(npar, nmc)

      ; Here's where the actual fitting takes place....
      print, format='($,"MC Iteration =")'
      for i = 0, nmc-1 do begin  ; Loop through Monte Carlo runs
        if i mod 5 eq 0 then $  ; Print update every 5th run
           print, format = '($," ",I3)', i
        si = smc[*,i]                               ; Extract spectrum
        fs_fit, x, si, wts, parinit, dpar, parfiti, $ ; Fit the model!
                   wfiti, fiti, chii;, /markov
        parfit_mc[*,i] = parfiti ; Store parameters from this iteration
      endfor

      print, ''
      parfit = mean(parfit_mc, dim=2)   ; Calculate mean parameter values
      sigpar = stddev(parfit_mc, dim=2) ; Calculate variation in parameters
      
      mcmc,parfit, x, s, wts, dpar,bestpar
      parfit=bestpar

      vel = parfit[0]            ; Extract radial velocity from parameters
      dvel = sigpar[0]           ; Extract RV uncertainty from parameter variation
     
      fit = spec_model(x, parfit)          ; Produce a model from mean best-fit parameters
      wfit = poly(x, parfit[9:*])          ; Wavelength solution
      chi = fs_chi2(x, s, wts, parfit, dpar) ; Reduced chi^2 of mean best-fit  
    endelse

   endif else begin
    if keyword_set(init) then begin          ; Init only, generate model with parameter guesses
     fit = spec_model(x,parinit)           ; Create model spectrum
     wfit = poly(x, parinit[9:*])          ; Create wavelength solution
     ;print, wfit
     parfit = parinit                      ; Set best fit values equal to guesses
     chi = fs_chi2(x, s, wts, parfit, dpar)  ; Chi^2 of "fit"

    endif else begin              ; Else, actually fit the model...
     ; Create nmc randomized spectra for Monte Carlo
     b = rotate(replicate(1, nmc), 1)
     smc = s # b + randomn(seed, n_elements(s), nmc) * (u # b)

     ; Initialize parameter storage for all Monte Carlo runs
     parfit_mc = fltarr(npar, nmc)

     ; Here's where the actual fitting takes place....
     print, format='($,"MC Iteration =")'
     for i = 0, nmc-1 do begin  ; Loop through Monte Carlo runs
        if i mod 5 eq 0 then $  ; Print update every 5th run
           print, format = '($," ",I3)', i
        si = smc[*,i]                               ; Extract spectrum
        fs_fit, x, si, wts, parinit, dpar, parfiti, $ ; Fit the model!
                   wfiti, fiti, chii
        parfit_mc[*,i] = parfiti ; Store parameters from this iteration
     endfor

     print, ''
     parfit = mean(parfit_mc, dim=2)   ; Calculate mean parameter values
     sigpar = stddev(parfit_mc, dim=2) ; Calculate variation in parameters
     
     vel = parfit[0]            ; Extract radial velocity from parameters
     dvel = sigpar[0]           ; Extract RV uncertainty from parameter variation
     
     fit = spec_model(x, parfit)          ; Produce a model from mean best-fit parameters
     wfit = poly(x, parfit[9:*])          ; Wavelength solution
     chi = fs_chi2(x, s, wts, parfit, dpar) ; Reduced chi^2 of mean best-fit
    endelse
  endelse

  r = s - fit                   ; Observed-Model residuals
  
  ; Plot best fit model, if requested
  if keyword_set(plot) then $
     plot_specfit, wfit, s, u, fit, parfit, $
                   objname = objname, jd = sxpar(obj.head, 'MJD'), ps=ps
     ;writecol,'fit.txt',wfit,s,fit,fmt='(d,d,d)' ;write fit to out file
  
  ; Print summary info to terminal
  print, ''
  print, 'S/N = ', median(s/u)
  print, 'Mean of residuals = ', mean(r)
  print, 'Std. Dev of resdiduals = ', stddev(r)
  print, 'chi^2 = ', chi
  print, 'parfit = ',parfit
  print, ''

  ; Write to procedure outputs
  wout = wfit
  sout = s

  if keyword_set(write) then begin
    file_delete,objname+'_pars.dat'
    wigrparfile, objname, parfit
  endif

end
