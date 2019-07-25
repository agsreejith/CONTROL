

pro cute_light_curve,file_list,lightcurve,wave_region=wave_region

  file=file_search(file_list)
  if (n_elements(file) eq 0) then begin
    logprint,'The file list is empty. Please check your input',logonly=logonly
    message,'The file list is empty. Please check your input'
  endif
 
  lc1    = dblarr(4,n_elements(file_list))
  lc2    = dblarr(4,n_elements(file_list))
  lc3    = dblarr(4,n_elements(file_list))
  lc_all = dblarr(4,n_elements(file_list))
  
  mg2=where(wave_region eq 'MgII')
  if mg2 ge 0 then lc_mg2 = dblarr(4,n_elements(file_list))
  mg1=where(wave_region eq 'MgI')
  if mg1 ge 0 then lc_mg1 = dblarr(4,n_elements(file_list))
  fe2=where(wave_region eq 'FeII')
  if fe2 ge 0 then lc_fe2 = dblarr(4,n_elements(file_list))
  
  for i=0,n_elements(file_list)-1 do begin
    wclf_spectrum_file=mrdfits(file_list[i],1,spectrum_hdr)
    if  fix(SXPAR(spectrum_hdr, 'FCALFLG',MISSING=-1)) eq -1 then begin
      logprint,'CONTROL: The FITS file specified('+file_list[i]+') is not flux calibrated.'
      logprint,'Press q to skip the creation of default light curves. Press any key to continue with light curve generation.'
      R = GET_KBRD()
      if R eq 'q' then begin
        logprint,'CONTROL: Skipping the light curve preperation.'
        goto,spectrum_loop_end
      endif
    endif
    nx=n_elements(wclf_spectrum_file.wave)
    if (tag_exist(wclf_spectrum_file,'flux')eq 1) then spectra_1dwf=wclf_spectrum_file.flux else if (tag_exist(wclf_spectrum_file,'counts') eq 1) then spectra_1dwf=wclf_spectrum_file.counts
    if (tag_exist(wclf_spectrum_file,'error')eq 1) then error_1dwf=wclf_spectrum_file.error 
    lenby3=fix(nx/3)
    lenby33=fix(nx-2*lenby3)
    st=0
    en=nx-1
    n_w=en-st+1
    x=make_array(n_w, /INDEX, /NOZERO, START=st)
    y=spectra_1dwf[st:en]
    dy=error_1dwf[st:en]
    tp=trapz_error(x,y,dy)
    lc_all[1,i]=tp[0]
    lc_all[2,i]=tp[1]
    lc_all[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    st=0
    en=lenby3-1
    n_w=en-st+1
    x=make_array(n_w, /INDEX, /NOZERO, START=st)
    y=spectra_1dwf[st:en]
    dy=error_1dwf[st:en]
    tp=trapz_error(x,y,dy)
    lc1[1,i]=tp[0]
    lc1[2,i]=tp[1]
    st=lenby3
    en=2*lenby3-1
    n_w=en-st+1
    x=make_array(n_w, /INDEX, /NOZERO, START=st)
    y=spectra_1dwf[st:en]
    dy=error_1dwf[st:en]
    tp=trapz_error(x,y,dy)
    lc2[1,i]=tp[0]
    lc2[2,i]=tp[1]
    st=2*lenby3
    en=nx-1
    n_w=en-st+1
    x=make_array(n_w, /INDEX, /NOZERO, START=st)
    y=spectra_1dwf[st:en]
    dy=error_1dwf[st:en]
    tp=trapz_error(x,y,dy)
    lc3[1,i]=tp[0]
    lc3[2,i]=tp[1]
    lc1[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    lc2[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    lc3[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    if mg2 ge 0 then begin
      ;near = Min(Abs(vector - number), index)
      ;2795,2802
      val=Min(Abs(wclf_spectrum_file.wave - 2793), st)
      val2=Min(Abs(wclf_spectrum_file.wave - 2805), en)
      n_w=en-st+1
      x=make_array(n_w, /INDEX, /NOZERO, START=st)
      y=spectra_1dwf[st:en]
      dy=error_1dwf[st:en]
      tp=trapz_error(x,y,dy)
      lc_mg2[1,i]=tp[0]
      lc_mg2[2,i]=tp[1]
      lc_mg2[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    endif
    if mg1 ge 0 then begin
      ;2852
      ;near = Min(Abs(vector - number), index)
      val=Min(Abs(wclf_spectrum_file.wave - 2850), st)
      val2=Min(Abs(wclf_spectrum_file.wave - 2854), en)
      n_w=en-st+1
      x=make_array(n_w, /INDEX, /NOZERO, START=st)
      y=spectra_1dwf[st:en]
      dy=error_1dwf[st:en]
      tp=trapz_error(x,y,dy)
      lc_mg1[1,i]=tp[0]
      lc_mg1[2,i]=tp[1]
      lc_mg1[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    endif
    if fe2 ge 0 then begin
      ;2585
      ;near = Min(Abs(vector - number), index)
      val=Min(Abs(wclf_spectrum_file.wave - 2583), st)
      val2=Min(Abs(wclf_spectrum_file.wave - 2587), en)
      n_w=en-st+1
      x=make_array(n_w, /INDEX, /NOZERO, START=st)
      y=spectra_1dwf[st:en]
      dy=error_1dwf[st:en]
      tp=trapz_error(x,y,dy)
      lc_fe2[1,i]=tp[0]
      lc_fe2[2,i]=tp[1]
      lc_fe2[0,i]=double(SXPAR( spectrum_hdr, 'OBS_TIME'))
    endif
 endfor
lightcurve={time:reform(lc_all[0,*]),full_data:reform(lc_all[1,*]),full_error:reform(lc_all[2,*]),short_data:reform(lc1[1,*]),$
            short_error:reform(lc1[2,*]),middle_data:reform(lc2[1,*]),middle_error:reform(lc2[2,*]),long_data:reform(lc3[1,*]),long_error:reform(lc3[2,*])}
if mg2 ge 0 then begin
  lightcurve=create_struct(lightcurve, 'MgII_data', reform(lc_mg2[1,*]))
  lightcurve=create_struct(lightcurve,'MgII_error', reform(lc_mg2[2,*]))
endif
if mg1 ge 0 then begin
   lightcurve=create_struct(lightcurve,'MgI_data', reform(lc_mg1[1,*]))
   lightcurve=create_struct(lightcurve,'MgI_error', reform(lc_mg1[2,*]))
endif
if fe2 ge 0 then begin
   lightcurve=create_struct(lightcurve,'FeII_data', reform(lc_fe2[1,*]))
   lightcurve=create_struct(lightcurve,'FeII_error', reform(lc_fe2[2,*]))
endif


spectrum_loop_end:
end