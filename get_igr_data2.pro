pro get_igr_data2, file, wav, spec, con, sig, out = out

 readcol,file,ind,wav,spec,snr,format='d,d,d,d'
 order=where(ind GE 15)
 wav[0:order[0]-1]=0
 start=where(wav GE 2.295)
 if start[0]+1023 NE start[1023] then begin
 	order=where(ind GE 16)
 	wav[0:order[0]-1]=0
 	start=where(wav GE 2.295)
 endif
 start=start[0]
 last=start[0]+1023

 ind=ind[start:last]
 wav=wav[start:last]*10000
 spec=spec[start:last]
 snr=snr[start:last]
 ;30814:31841 first area
 ;32862:33885 second area

 contf,spec,con,sbin=20,nord=4,frac=.2

 readn=8.93
 gain=2.21
 sigman=double((spec-con)*gain+readn^2)
 sigman=double(sqrt(sigman^2))
 sig=double(sqrt(sigman*snr)/2.21)

 if keyword_set(out) then begin
    writecol,'igrinspars.txt',wav,spec,con,sig,fmt='(d,d,d,d)'
 end

end