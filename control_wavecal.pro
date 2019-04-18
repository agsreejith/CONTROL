function control_wavecal,in_image,header,infile,wavecal_type
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    logprint,'Syntax - control_wavecal,in_image,infile,wavecal_type',logonly = logonly
    message,'Syntax - control_wavecal,in_image,infile,wavecal_type'
    err = '%CONTROL_WAVECAL: Insufficient number of parameters'
  endif
  ;Constants
  R_sun=6.957d10          ; in cm
  L_sun=3.828d33
PRINT,in_image[0]
if wavecal_type ne 'simple' then crosscor=1 else crosscor=0
wave_path=detectos(infile.wavecal_file)
if (file_test(wave_path) eq 0) then begin
  logprint,'CONTROL: Wavelength file not found. Please re-run the pipeline with actual file address',logonly=logonly
  message,'CONTROL: Wavelength file not found. Please re-run the pipeline with actual file address'
endif
readcol,wave_path,wavelength1,F='D'
wave_cal_flg=0
if crosscor eq 1 then begin
  input_flux=in_image
  target=SXPAR( header, 'TARGNAME')
  sp_path=detectos(infile.stellar_params)
  if file_test(sp_path) ne 1 then begin
    in_file=detectos(infile.model_file)
    if file_test(in_file) ne 1 then begin
      logprint,'CONTROL_WAVECAL: Input model file for corss corelation not found. Program will continue without cross corelation'
      goto,withoutcorss
    endif
    flux=dblarr(2,file_lines(in_file))
    readcol,in_file,flux[0,*],flux[1,*],F='D,D'
  endif else begin
    readcol,sp_path,star,temperature,radius,magnitude,F='A,D,D,D'
    loc=where(star eq target)
    in_file=detectos(infile.model_file)
    t_star=temperature[loc]
    V_mag=magnitude[loc]
    r_star=radius[loc]
    t_star=10000
    V_mag= 2.36
    r_star= 7.55
    t=fix(t_star)
    file=in_file+'\models\t'+string(t, Format='(I05)') +'g4.4\model.flx' ;assuming folder named by their temperature and file are named as model.flx
    if file_test(file) ne 1 then t = t+100 ;above 8000K the steps is 200K
    file=in_file+'\models\t'+string(t, Format='(I05)') +'g4.4\model.flx' ;test again to see if temperature is not in range or is in steps of 100 or 200
    if file_test(file) ne 1 then message,'Error CONTROL_WAVECAL: Invalid input of stellar temperature in stellar temperature file, Please round the temperature to nearest 100'
    length=file_lines(file)
    fdata=dblarr(3,length)
    openr,1,file
    readf,1,fdata
    close,1
    ;definitions
    flux1=dblarr(3,length)
    r_star=double(r_star*R_sun)
    flux=fdata
    flux1[1,*]=(3d18*fdata[1,*])/(fdata[0,*]*fdata[0,*]) ;convert to ergs/cm2/s/A
    flux1[2,*]=(3d18*fdata[2,*])/(fdata[0,*]*fdata[0,*]) ;convert to ergs/cm2/s/A
    flux[1,*]=flux1[1,*]*4*!pi*(r_star^2)*4*!pi ;convert to ergs/s/A second 4*!pi for steradian conversion
    flux[2,*]=flux1[2,*]*4*!pi*(r_star^2)*4*!pi  ;convert to ergs/s/A second 4*!pi for steradian conversion
    t=double(t_star)
    t4=0.0D
    t4=t^4
    r2=r_star^2
    stepahs=5.6704d-5
    L = stepahs*4*!pi*r2*t4
  endelse  
;  file_eff='extra\eff_area.txt'
;  length=file_lines(file_eff)
;  eff_area=dblarr(5,length)
;  openr,1,file_eff
;  readf,1,eff_area
;  close,1
;  aeff=interpol(eff_area[4,*],eff_area[0,*],wavelength1,/SPLINE)
    
  nx=n_elements(wavelength1)
  wave_res=dblarr(n_elements(wavelength1)/2)
  for i=0,nx-2,2 do wave_res[i/2]=mean(wavelength1[i:i+1])
  new_flux=interpol(flux[1,*],flux[0,*],wavelength1,/SPLINE)
  flux_ref = new_flux2
  flux_test = input_flux
  lmin = wavelength1[10]
  lmax = wavelength1[n_elements(wavelength1)-10]
  delta=cross_correlate(wavelength1, new_flux2, wavelength1, flux_test, lmin, lmax, CCF=ccf)
  wavelength=wavelength1-delta
  wave_cal_flg=1
  wave_shift=delta
  logprint,'CONTROL_WAVECAL: Wavelength calibration carried out using cross-correlation with a synthetic spectrum'
  sxaddpar, header, 'WCAL','Cross-correlation with synthetic','Wavelength Calibration'
  sxaddpar, header, 'WSHFT',wave_shift, 'Wavelength Cross-correlation applied'
endif
withoutcorss:
sxaddpar, header, 'WCALFLG',wave_cal_flg , 'Wavelength Calibration Flag'
spectrum={wavelength:wavelength,flux:in_image,wshift:wave_shift,header:header}
;stop
return,spectrum
end