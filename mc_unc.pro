pro mc_unc, param_file, mcmc=mcmc

 
 readcol, param_file, par, format='a'

 date=par[2]
 id=date+'00'

 mjd=par[7]-2400000.5 
 ra=par[3]
 dec=par[4]
 hc_mjd=HELIO_JD( mjd, ra, dec)


 get_igr_data2,par[0],wav,spec,con,sig

 a=size(spec)
 
 ;openw, 3, 'igr_unc.txt'
     openw, 3, 'MCMC_unc.txt'
 
 FOR k = 1, 100 DO BEGIN
 
    r=randomn(seed,a[1],1)

    specm=r*sig+spec
;Program desined to do MC simulations for IGRINS data
;rough version, written by JNM, May 2015

    mkhdr, hdr, spec

        sxaddpar, hdr, 'INSTRUME', 'IGRINS'
        sxaddpar, hdr, 'OBJECT', par[1]
        sxaddpar, hdr, 'RA', ra, 'Right Ascension of object in hours (2000.0)'
        sxaddpar, hdr, 'DEC', dec, 'Declination of objet in degrees (2000.0)'
        sxaddpar, hdr, 'MJD', mjd, 'Modified Julian Date (JD-2400000) at middle of exposure', $
           form='(F13.6)'
        sxaddpar, hdr, 'ITIME', par[5], 'Integration time (seconds)'
        sxaddpar, hdr, 'OBSLON', 104.0225, 'Longitude of observatory (degrees)'
        sxaddpar, hdr, 'OBSLAT', 30.6714, 'Latitude of observatory (degrees)'
        sxaddpar, hdr, 'OBSALT', 2070, 'Altitude of observatory (meters)'
        sxaddpar, hdr, 'HC_MJD', hc_mjd, 'Heliocentric MJD at middle of exposure', form='(F13.6)'
        sxaddpar, hdr, 'HELCOR', par[6], 'Barycentric correction (km/s)'
        sxaddpar, hdr, 'BEAM', 'a', 'Telescope nod position (A/B)'
        sxaddpar, head, 'GROUP', 1
        sxaddpar, hdr, 'ORDPWR', 2 
        sxaddpar, hdr, 'MASH', 0
        sxaddpar, hdr, 'HORNE', 1
        sxaddpar, hdr, 'SFPOLY', 2
 
        fileout=par[1]+'_'+id+'.ech'
 
        wdech, fileout, hdr, spec , $
                 sig=sig, $
                 wave=wav, $
                 cont=con, $
                 /overwrite
    
    print,' '
    print,'----------------------------'
    print,'Beginning iteration '+strcompress(string(mean(k)),/remove_all)
    print,'----------------------------'
    print,' '

    if keyword_set(mcmc) then fit_spec_all,par[1],/mcmc else fit_spec_all,par[1]
 
    restore,par[1]+'_rvfit.sav'
    ;printf, 3, f.v
    printf, 3, f.p(1), f.p(3), f.v

    print,'Output written'

    file_delete,par[1]+'_'+id+'.ech'
    file_delete,par[1]+'_rvfit.sav'
ENDFOR

close,3
END