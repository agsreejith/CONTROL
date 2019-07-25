; NAME:
;      CONTROL_DARK_COMBINE
;
; PURPOSE:
;      Returns master bias frames from input bias files based on mean, median or mode combine as chosen by the user;
; CALLING SEQUENCE:
;      CONTROL_DARK_COMBINE,dark_list,mdark,mbias_str,type=type,sat_value=sat_value,threshold=threshold
;
; INPUTS:
;      dark_list  = The file containing the location of dark files
;      mbias_str  = Master bias structure with the following structure im:master bias data,hdr:master bias header,dq:master bias data quality


; OPTIONAL INPUTS
;      type      = Type of combining biases (mean, median, mode). Default is median.
;      sat_value = Expected saturation point for images. Default is 72000.
;      threshold = Value specifying acceptable deviation in the frame for data quality. Default is 5. 
;                  For example a threshold of 5 will flag all pixels whose value deviates from mean of master bias by 5 sigma.     
; 
; OUTPUT:
;      mdark     = Created master dark image structure with image, header and data quality flag.
;
    
; PROCEDURE:
;      Master dark creator for CUTE;
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;      modified 05.07.2019 by A. G. Sreejith
;;#######################################################################


pro control_dark_combine,dark_list,mdark,mbias_str,type=type,sat_value=sat_value,threshold=threshold
idl_ver=float(!Version.RELEASE)
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    logprint,'CONTROL_DARK_COMBINE: Syntax - control_dark_combine,dark_list,mdark,mbias_file',logonly = logonly
    message,'CONTROL_DARK_COMBINE: Syntax - control_dark_combine,dark_list,mdark,mbias_file'
    err = '%control_dark_combine: Insufficient number of parameters'
    return
  endif
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_DARK_COMBINE: Type input for dark combine not defined uisng default: median'
    type='median'
  endif
  if keyword_defined(sat_value) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Saturation value set to default:72000'
    sat_value = 72000
endif
  if double(sat_value) eq 0 then begin
    logprint,'Saturation value cannot be zero setting it to default of 72000'
  endif

  if keyword_defined(threshold) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Threshold for deviation set to default: 5 sigma'
    threshold = 5
  endif
  if double(threshold) eq 0 then begin
    logprint,'Threshold for deviation value cannot be zero setting it to default of 5 sigma.'
  endif  

  dummy=mrdfits(dark_list[0],0,hdr,/SILENT)
  nxyd=size(dummy)
  totpix=n_elements(dummy)
  dq_arr=bytarr(nxyd[1],nxyd[2])
  limit=0.001*totpix
  ;read in master bias
  if keyword_defined(mbias_str) then begin
    if ISA(mbias_str) eq 0 then begin
      logprint,'CONTROL_DARK_COMBINE: Requires a master bias file'
      logprint,'CONTROL_DARK_COMBINE: Press any key except q to create a master bias frame with zeros or press q to quit.'
      R = GET_KBRD()
      if R eq 'q' then begin
        logprint,'CONTROL_DARK_COMBINE: Exiting as requested by the user.',logonly = logonly
        message,'CONTROL_DARK_COMBINE: Exiting as requested by the user.'
        err = '%control_dark_combine: Insufficient number of parameters'
        return
      endif else begin
        mbias=make_array(nxyd[1],nxyd[2],/DOUBLE,value=0.0)
        mbias_flag = 0
        r=12.25
        mbdq=bytarr(nxyd[1],nxyd[2])
      endelse  
    endif else begin
      mbias=mbias_str.im 
      mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
      nxyb=size(mbias)
      r=float(SXPAR( mbias_str.hdr, 'RNOISE'))
      mbdq=mbias_str.dq
      ;to be implemented later
;     if nxyd[1] ne nxyb[1] then  
;     ycut1=SXPAR( hdr, 'YCUT1')
;     ycut2=SXPAR( hdr, 'YCUT2')
;     ylen=ycut2-ycut1+1
;     if (ylen ne ny) then begin
;
;      hbmap_nw=[*,ycut1:ycut1+ny-1]
    endelse
  endif else begin
    mbias=make_array(nxyd[1],nxyd[2],/DOUBLE,value=0.0)
    mbias_flag = 0
    r=12.25
    mbdq=bytarr(nxyd[1],nxyd[2])
  endelse
  n=n_elements(dark_list)
  if n ne 0 then begin
    dark_ar=dblarr(nxyd[1],nxyd[2],n)
    dark_err=dblarr(nxyd[1],nxyd[2],n)
    include=intarr(n)
    for i=0, n-1 do begin
      filename=dark_list[i]
      dark=mrdfits(filename,0,hdr,/SILENT)
      FITS_INFO, filename,N_ext =numext,/SILENT
      if numext eq 2 then dq=mrdfits(filename,1,hdr,/SILENT) else dq=bytarr(nxyd[1],nxyd[2])
      dark_nw=rejection(dark,threshold,npix)
      if npix lt limit then begin
        dark_ar[*,*,i]=dark_nw-mbias
        dark_err[*,*,i]= dark_nw+r^2
        include[i]=i
      endif else include[i]=-1
    dq_arr+=dq
    endfor
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    dark_arr=dblarr(nxyd[1],nxyd[2],n_frames)
    dark_er=dblarr(nxyd[1],nxyd[2],n_frames)
    dark_arr=dark_ar(*,*,incl)
    dark_er=dark_err(*,*,incl)  
    if n_frames eq 0 then begin
      logprint,'CONTROL DARK COMBINE: No valid DARK file found. Terminating MASTER DARK creation.'
      return
    endif
    if n_frames le 1 then begin
      logprint,'CONTROL DARK COMBINE: Only one valid DARK file found. Do you wnat to assume it as the MASTER DARK?.'
      logprint,'Press q to skip this assumption. Press any key to continue with MASTER DARK creation with one valid DARK file.'
      R = GET_KBRD()
      if R eq 'q' then begin
        logprint,'CONTROL DARK COMBINE: Terminating MASTER DARK creation as requested by the user.'
        return
      endif
    endif
     logprint,'CONTROL DARK COMBINE: Combining '+STRTRIM(STRING(n_frames),2)+' DARK file to create MASTER DARK using '+type+' method .'
    case type of
      'median' : begin
                 if (idl_ver ge 8.3) then mdark_val = median(dark_arr,dimension=3,/even) else begin
                   totalImage = Total(dark_arr, 3,/NAN)
                   minValue = Min(dark_arr, Dimension=3,/NAN)
                   maxValue = Max(dark_arr, Dimension=3,/NAN)
                   mdark_val = totalImage - minValue - maxValue
                 endelse
                 sigma_db=total(dark_er,3,/NAN); to modify
                 mdark_err=sqrt(sigma_db/n_frames^2)
                 end
      'mean' : begin
               if (idl_ver ge 8.3) then mdark_val = mean(dark_arr,dimension=3,/NAN) else begin
                 totalImage = Total(dark_arr, 3,/NAN)
                 mdark_val = totalImage/n_frames
               endelse
               sigma_db=total(dark_er,3,/NAN)
               mdark_err=sqrt(sigma_db/n_frames^2)
               end  
      'mode' : begin
               mdark_val = mode(dark_arr,dimension=3,/NAN)
               sigma_db=total(dark_er,3,/even,/NAN) ;to modify
               mdark_err=sqrt(sigma_db/n_frames^2)
               end  
      else : begin
              errorlog,'Invalid type input for combine: Please recheck your input',logonly=1
              message,'Invalid type input for combine: Please recheck your input'
           end  
    endcase
    ;  if type = 'med' then mdark = median(dark_arr,dimension=3,/even,/NAN)
    ;  else if type = 'mea' then mdark =mean(dark_arr,dimension=3,/even,/NAN)
    ;  else if type = 'mod' then mdark =median(dark_arr,dimension=3,/even,/NAN)
    ;  
    ;Header defnitions
    t=SXPAR(hdr,'TIME_IN')
    exptime=SXPAR(hdr,'EXP_TIME')
  
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
    dq_arr=dq_arr+mbdq
    prb=where(dq_arr ge 1)
    dark_dq=bytarr(nxyd[1],nxyd[2])
    dark_dq[prb]=1
    ;checks
    ;saturated pixels
    
    sat_loc = where(mdark_val ge sat_value)
    if total(sat_loc) eq -1 then sat_flag = 0 else sat_flag = 1
     dark_dq[sat_loc]=1
    ;deviation
    std=stddev(mdark_val)
    std_loc = where((mdark_val ge (mean(mdark_val)+threshold*std)) or (mdark_val le (mean(mdark_val)-threshold*std)))
    if total(std_loc) eq -1 then std_flag = 0 else std_flag = 1
    dark_dq[std_loc]=1
    
    nan_loc = where(finite(mdark_val, /NAN) eq 1)
    if total(nan_loc) eq -1 then nan_flag = 0 else nan_flag = 1
    dark_dq[nan_loc]=1

    mdark_flag = sat_flag+std_flag+mbias_flag+nan_flag

    if mdark_flag gt 1 then mdark_flag=1
    ;update flag headers
    ;sxaddpar, dhdr, 'CRFLG', cr_flag
    sxaddpar, dhdr, 'SRNFLG', sat_flag
    sxaddpar, dhdr, 'STDFLG', std_flag
    sxaddpar, dhdr, 'MBIFLG', mbias_flag
    sxaddpar, dhdr, 'MDRFLG', mdark_flag
    logprint,'CONTROL DARK COMBINE: Saturated (values above '+STRTRIM(STRING(sat_value),2)+') pixels and pixels that deviate by '+STRTRIM(STRING(threshold),2)+' sigma from mean of  MASTER DARK have been flaged.'
    mdark={im:mdark_val,error:mdark_err,hdr:dhdr,dq:dark_dq}

    return
  ;mwrfits,mdark,input_file+'m_bias.fits',dhdr,/create
  endif else begin
    logprint,'CONTROL_DARK_COMBINE:Dark list is empty, exiting dark combine without creating master dark'
    return
  endelse
end