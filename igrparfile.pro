pro igrparfile, name, rv

  ; Find all ECH files for this target
  files = findfile(name+'_*.ech')
  
  openw, 1, name+'_pars.dat'    ; Open parameter file for writing
  printf, 1, '0 15'             ; Default mask boundaries
  printf, 1, strcompress(string(n_elements(files)), /remove_all) ; Number of observations
 
  ; Default parameter step sizes:
  printf, 1, '              1.    0.05   0.0   0.05  0.0   0.1   0.1  0.01  0.001     0.2  0.001  0.000001'

  ; Loop through ECH files and print line of default guesses for 
  ;  each obs id:
  for i=0, n_elements(files)-1 do begin
     id = strmid(files[i], strlen(files[i])-14, 10) ; Parse ID from file name
     printf, 1, id+'    '+rv+'   0.4    0.0   0.5   0.0   7.5   1.   0.0   0.0   23031.0  0.1625  0.0'
  endfor
  
  ; Print default paths for telluric and photosphere models:
  ;printf, 1, '~/data/spec_templates/spectra_nso/extracted/irtf/nso_sky_22940_23020.sav'
  ;printf, 1, '~/data/spec_templates/synthmag/spec/ng385000_b000000_vs00.0.spec'
  printf, 1, '~/data/spec_templates/spectra_nso/extracted/nirspec/nso_sky_22800_23300.sav'
  printf, 1, '~/data/spec_templates/synthmag/spec/ng385000_b000000_w22800_23300.spec'
  printf, 1, ''
  close, 1
end