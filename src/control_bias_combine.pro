; NAME:
;      CONTROL_BIAS_COMBINE
;
; PURPOSE:
;      Returns master bias frames from input bias files based on mean, median or mode combine as 
;      chosen by the user
;      
; CALLING SEQUENCE:
;      CONTROL_BIAS_COMBINE,bias_list,mbias,type,sat_value,threshold
;
; INPUTS:
;      bias_list = The array containing list of bias file locations. 
; 
; OPTIONAL INPUTS
;      type      = Type of combining biases (mean, median, mode). Default is median.
;      sat_value = Expected saturation limit for the images. Default is 72000.
;      threshold = Value specifying acceptable deviation in the frame for data quality. Default is 
;                  5. For example a threshold of 5 will flag all pixels whose value deviates from 
;                  mean of master bias by 5 sigma.     
; 
; OUTPUT:
;      mbias     = Created master bias image structure with image, header and data quality flag.
;  
; PROCEDURE:
;      Master bias creator for CUTE;
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;      modified 12.03.2019 by A. G. Sreejith
;      modified 05.07.2019 by A. G. Sreejith
;##################################################################################################      

pro control_bias_combine,bias_list,mbias,type=type,sat_value=sat_value,threshold=threshold
  
  idl_ver=float(!Version.RELEASE)
  err = ''
  if N_params() LT 2 then begin             ;Need at least 2 parameters
    logprint,'CONTROL_BIAS_COMBINE: Syntax - control_bias_combine,bias_list,mbias',$
             logonly = logonly
    message,'CONTROL_BIAS_COMBINE: Syntax - control_bias_combine,bias_list,mbias
    err = '%control_bias_combine: Insufficient number of parameters'
    return
  endif
  
  ;Defining the defaults
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Type input for bias combine not defined uisng default: median'
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
  
  n=n_elements(bias_list)
  if n eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Bias list is empty,'$
             +' exiting bias combine without creating master bias'
    return
  endif else begin
    gain_val=dblarr(n)
    bias=mrdfits(bias_list[0],0,hdr,/SILENT)
    nxy=size(bias)
    nx=nxy[1]
    ny=nxy[2]
    FITS_INFO, bias_list[0],N_ext =numext,/SILENT
    if numext eq 2 then dq=mrdfits(bias_list[0],1,hdr,/SILENT) else dq=bytarr(nx,ny)
    dq_arr=bytarr(nx,ny)
    bias_ar=dblarr(nx,ny,n)
    include=intarr(n)
    totpix=n_elements(bias)
    limit=0.001*totpix     ;Condition that checks for a percentage of bad pixels.
    
    for i=0, n-1 do begin
      filename=bias_list[i] ;modify based on file structure
      bias=mrdfits(filename,0,hdr,/SILENT)
      FITS_INFO, filename,N_ext =numext,/SILENT
      if numext eq 2 then dq=mrdfits(filename,1,hdr,/SILENT) else dq=bytarr(nx,ny)
      bias_nw=rejection(bias,threshold,npix)
      if npix lt limit then begin
        bias_ar[*,*,i]=bias_nw
        gain_val[i]=SXPAR( hdr, 'CCDGAIN')
        include[i]=i
      endif else include[i]=-1
      dq_arr+=dq
    endfor
    
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    bias_arr=dblarr(nx,ny,n_frames)
    bias_arr=bias_ar(*,*,incl)
    
    if n_frames eq 0 then begin
      logprint,'CONTROL BIAS COMBINE: No valid BIAS file found. Terminating MASTER BIAS creation.'
      return
    endif
    
    if n_frames le 1 then begin
      logprint,'CONTROL BIAS COMBINE: Only one valid BIAS file found.'$
               +' Do you wnat to assume it as the MASTER BIAS?.'
      logprint,'Press q to skip this assumption.'$
               +' Press any key to continue with MASTER BIAS creation with one valid BIAS file.'
      R = GET_KBRD()
      if R eq 'q' then begin
        logprint,'CONTROL BIAS COMBINE: Terminating MASTER BIAS creation as requested by the user.'
        return
      endif  
    endif
 
    logprint,'CONTROL BIAS COMBINE: Combining '+STRTRIM(STRING(n_frames),2)+$
             ' BIAS file to create MASTER BIAS using '+type+' method .'
    case type of
      'median':begin
                  if (idl_ver ge 8.3) then mbias_val = median(bias_arr,dimension=3,/even) $
                  else begin
                    totalImage = Total(bias_arr, 3,/NAN)
                    minValue   = Min(bias_arr, Dimension=3,/NAN)
                    maxValue   = Max(bias_arr, Dimension=3,/NAN)
                    mbias_val  = totalImage - minValue - maxValue
                  endelse
                end  
      'mean' :begin
                if (idl_ver ge 8.3) then mbias_val = mean(bias_arr,dimension=3,/NAN) $
                else begin
                  totalImage = Total(bias_arr, 3,/NAN)
                  mbias_val  = totalImage/n_frames
                endelse
              end  
      'mode' : mbias_val = mode(bias_arr,dimension=3,/NAN)
    else : begin
              errorlog,'Invalid type input for combine: Please recheck your input',logonly=1
              message,'Invalid type input for combine: Please recheck your input'
           end   
    endcase
    ;Basic caclucation of readnoise. Now redundant
    ;r_noise= stddev(mbias_val,/NAN)/mean(gain_val)
    
  
    t=SXPAR(hdr,'EXP_TIME')
    r_noise=SXPAR(hdr,'RD_NOISE') 
    if datatype(r_noise,FLAG = 2) eq 0 then r_noise=12.25
    
    ;Header defnitions
    sxaddpar, bihdr, 'Time_in_JD', t
    sxaddpar, bihdr, 'RNOISE', r_noise, 'Readout noise'
    sxaddpar, bihdr, 'SIGMA',stddev(mbias_val,/NAN), 'Standard deviation of the frame'
    sxaddpar, bihdr, 'MEAN', mean(mbias_val,/NAN), 'Mean value of the frame'
    sxaddpar, bihdr, 'MEDIAN ', median(mbias_val), 'Median value of the frame'
    sxaddpar, bihdr, 'MAX', max(mbias_val,/NAN), 'Maximum value of the frame'
    sxaddpar, bihdr, 'MIN', min(mbias_val,/NAN), 'Minimum value of the frame'
    sxaddpar, bihdr, 'NFRAMES', n_frames, 'Number of frames used in bias combine'
    
    ;Necessary checks
    ;Data quality
    prb=where(dq_arr ge 1)
    bias_dq=bytarr(nx,ny)
    if total(prb) ne -1 then bias_dq[prb]=1
    
    ;Saturated pixels
    sat_loc = where(mbias_val ge sat_value)
    if total(sat_loc) eq -1 then sat_flag = 0 else sat_flag = 1
    if total(sat_loc) ne -1 then bias_dq[sat_loc]=1
    
    ;Deviation 
    std=stddev(mbias_val)
    std_loc = where((mbias_val ge (mean(mbias_val)+threshold*std)) or $
                   (mbias_val le (mean(mbias_val)-threshold*std)))
    if total(std_loc) eq -1 then std_flag = 0 else std_flag = 1
    if total(sat_loc) ne -1 then bias_dq[std_loc]=1
    nan_loc = where(finite(mbias_val, /NAN) eq 1)
    if total(nan_loc) eq -1 then nan_flag = 0 else nan_flag = 1
    if total(nan_loc) ne -1 then bias_dq[nan_loc]=1

    mbias_flag = sat_flag+std_flag+nan_flag
    if mbias_flag gt 1 then mdark_flag=1

    ;Update flag in header
    sxaddpar, bihdr, 'SRNFLG', sat_flag
    sxaddpar, bihdr, 'STDFLG', std_flag
    sxaddpar, bihdr, 'MBIFLG', mbias_flag
    
    logprint,'CONTROL BIAS COMBINE: Saturated pixels (values above '+STRTRIM(STRING(sat_value),2)+$
             ') and pixels that deviate by '+STRTRIM(STRING(threshold),2)+$
             ' sigma from mean of  MASTER BIAS have been flaged.'
    mbias={im:mbias_val,hdr:bihdr,dq:bias_dq}
  endelse
  return
end