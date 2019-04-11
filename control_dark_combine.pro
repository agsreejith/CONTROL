; NAME:
;      CONTROL_DARK_COMBINE
;
; PURPOSE:
;      Returns master bias frames from input bias files based on mean,median or mode combine as chosen by the user;
; CALLING SEQUENCE:
;      CONTROL_DARK_COMBINE,dark_list,type,mbias_file,mdark,dhdr
;
; INPUTS:
;      dark_list  = The file containing the location of dark files
;      type       = type of combining biases (mea = mean, med = median, mod = mode)
; 
; OPTIONAL INPUTS:
;      mbias_file = Location of master bias file
;      
; OUTPUT:
;      mdark      = Created master dark image
;      dhdr       = Header assosiated with the created master dark
;
; REQUIRES:
;     Hot and bad pixels have to be specified by NAN execute CONTROL_HOTBAD if needed.
;     
; PROCEDURE:
;      Master dark creator for CUTE;
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;


pro control_dark_combine,dark_list,type,mdark,mbias_str,sat_value,threshold
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    logprint,'CONTROL_DARK_COMBINE: Syntax - control_dark_combine,dark_list,type,mbias_file,mdark',logonly = logonly
    message,'CONTROL_DARK_COMBINE: Syntax - control_dark_combine,dark_list,type,mbias_file,mdark'
    err = '%control_dark_combine: Insufficient number of parameters'
    return
  endif
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_DARK_COMBINE: Type input for dark combine not defined uisng default: median'
    type='median'
  endif
  if keyword_defined(sat_value) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Saturation value set to default:7200'
    sat_value = 72000
  endif
  if keyword_defined(threshold) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Threshold for deviation set to default: 5 sigma'
    threshold = 5
  endif
  
  dummy=mrdfits(dark_list[0],0,hdr)
  nxyd=size(dummy)
  totpix=n_elements(dummy)
  limit=0.001*totpix
  ;read in master bias
  if keyword_defined(mbias_str) then begin
    mbias=mbias_str.im 
    mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
    nxyb=size(mbias)
    r=float(SXPAR( mbias_str.hdr, 'RNOISE'))
    ;to be implemented later
;    if nxyd[1] ne nxyb[1] then  
;    ycut1=SXPAR( hdr, 'YCUT1')
;    ycut2=SXPAR( hdr, 'YCUT2')
;    ylen=ycut2-ycut1+1
;    if (ylen ne ny) then begin
;
;      hbmap_nw=[*,ycut1:ycut1+ny-1]

  endif else begin
    mbias=make_array(nxyd[1],nxyd[2],/DOUBLE,value=0.0)
     mbias_flag = 0
     r=12.25
  endelse
  n=n_elements(dark_list)
  if n ne 0 then begin
    dark_ar=dblarr(nxyd[1],nxyd[2],n)
    dark_err=dblarr(nxyd[1],nxyd[2],n)
    include=intarr(n)
    for i=0, n-1 do begin
      filename=dark_list[i]
      dark=mrdfits(filename,0,hdr)
      dark_nw=rejection(dark,threshold,npix)
      if npix lt limit then begin
        dark_ar[*,*,i]=dark_nw-mbias
        dark_err[*,*,i]= dark_nw+r^2
        include[i]=i
      endif else include[i]=-1
    endfor
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    dark_arr=dblarr(nxyd[1],nxyd[2],n_frames)
    dark_er=dblarr(nxyd[1],nxyd[2],n_frames)
    dark_arr=dark_ar(*,*,incl)
    dark_er=dark_err(*,*,incl)  
    
    case type of
      'median' : begin
                 mdark_val = median(dark_arr,dimension=3,/even)
                 sigma_db=total(dark_er,3,/NAN); to modify
                 mdark_err=sqrt(sigma_db/n^2)
                 end
      'mean' : begin
               mdark_val = mean(dark_arr,dimension=3,/even,/NAN)
               sigma_db=total(dark_er,3,/NAN)
               mdark_err=sqrt(sigma_db/n^2)
               end  
      'mode' : begin
               mdark_val = mode(dark_arr,dimension=3,/even,/NAN)
               sigma_db=total(dark_er,3,/even,/NAN) ;to modify
               mdark_err=sqrt(sigma_db/n^2)
               end  
      else : message,'Invalid type input for combine: Please recheck your input'
    endcase
    ;  if type = 'med' then mdark = median(dark_arr,dimension=3,/even,/NAN)
    ;  else if type = 'mea' then mdark =mean(dark_arr,dimension=3,/even,/NAN)
    ;  else if type = 'mod' then mdark =median(dark_arr,dimension=3,/even,/NAN)
    ;  
    ;Header defnitions
    t=SXPAR(hdr,'TIME_IN')
    exptime=SXPAR(hdr,'EXPOSURE')
  
    r_noise=r
    sxaddpar, dhdr, 'Time_in_JD', t
    sxaddpar, dhdr, 'RNOISE', r_noise, 'Readout noise'
    sxaddpar, dhdr, 'SIGMA',stddev(mdark_val,/NAN), 'Standard deviation of the frame'
    sxaddpar, dhdr, 'MEAN', mean(mdark_val,/NAN), 'Mean value of the frame'
    sxaddpar, dhdr, 'MEDIAN ', median(mdark_val), 'Median value of the frame'
    sxaddpar, dhdr, 'MAX', max(mdark_val,/NAN), 'Maximum value of the frame'
    sxaddpar, dhdr, 'MIN', min(mdark_val,/NAN), 'Minimum value of the frame'
    sxaddpar, dhdr, 'EXPTIME', exptime, 'Exposure time of the frame'
    sxaddpar, dhdr, 'NFRAMES', n_frames, 'Number of frames used in dark combine'
    ;checks flags that ahve to be checked
    ;checks
    ;saturated pixels
    sat_value = 72000
    sat_loc = where(mdark_val ge sat_value)
    if total(sat_loc) eq -1 then sat_flag = 0 else sat_flag = 1
    ;deviation
    std=stddev(mdark_val)
    std_loc = where((mdark_val ge (mean(mdark_val)+threshold*std)) or (mdark_val le (mean(mdark_val)-threshold*std)))
    if total(std_loc) eq -1 then std_flag = 0 else std_flag = 1


    mdark_flag = sat_flag+std_flag+mbias_flag
    ;update flag headers
    ;sxaddpar, dhdr, 'CRFLG', cr_flag
    sxaddpar, dhdr, 'SRNFLG', sat_flag
    sxaddpar, dhdr, 'STDFLG', std_flag
    sxaddpar, dhdr, 'MBIFLG', mbias_flag
    sxaddpar, dhdr, 'MDRFLG', mdark_flag

    mdark={im:mdark_val,error:mdark_err,hdr:dhdr}

    return
  ;mwrfits,mdark,input_file+'m_bias.fits',dhdr,/create
  endif else begin
    logprint,'CONTROL_DARK_COMBINE:Dark list is empty, exiting dark combine without creating master dark'
    return
  endelse
end