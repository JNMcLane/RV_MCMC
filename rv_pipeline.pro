pro rv_pipeline, param_file,rv_guess

	readcol,param_file,par,format='a'

	igr_adapt, param_file

	igrparfile,par[1],rv_guess

	fit_spec, par[1], par[2]+'00', /plot, /write

	fit_spec, par[1], par[2]+'00', /mcmc, /write, /plot

	for cycles=1, 4 do begin
		fit_spec, par[1], par[2]+'00', /mcmc, /write
	end

	fit_spec_all, par[1], /mcmc

	file_mkdir, 'UNC'
	file_copy, [par[0], param_file, '*pars.dat'], 'UNC'
	cd, 'UNC'

	mc_unc,param_file,/mcmc
	
	readcol,'MCMC_unc.txt',a,b,c,format='d,d,d'
	RV=mean(c)
	dev=stddev(c)*1000

	print,' '
	print, 'Radial velocity is '+string(RV)+' km/s.'
	print, 'Standard deviation is '+string(dev)+' m/s.'
	print,' '

end