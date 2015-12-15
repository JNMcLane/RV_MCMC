;Program designed to adapt reduced IGRINS data
;to be used with fitspec

;Param File Format:
;
;Data file name
;Target Name
;Local Date
;RA
;Dec
;Integration Time
;BVC
;JD at middle of exposure


;___Written by JNM, Dec 2015___

pro igr_adapt, param_file

 
 readcol,param_file,par,format='a'

 date=par[2]
 id=date+'00'

 mjd=par[7]-2400000.5 
 ra=par[3]
 dec=par[4]
 hc_mjd=HELIO_JD( mjd, ra, dec)
 ;4.56448 citau ra
 ;22.8393 citau dec
 ;7.6564 gj281 ra
 ;2.2032 gj281 ra

 get_igr_data2,par[0],wav,spec,con,sig

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
 END

END