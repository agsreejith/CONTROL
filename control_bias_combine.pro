; NAME:
;      CONTROL_BIAS_COMBINE
;
; PURPOSE:
;      Returns master bias frames from input bias files based on mean,median or mode combine as chosen by the user;
; CALLING SEQUENCE:
;      CONTROL_BIAS_COMBINE,bias_list,type,mbias,bihdr
;
; INPUTS:
;      bias_list = The file containing list of bias files. 
;      type      = type of combining biases (mea = mean, med = median, mod = mode)
; ;
; OUTPUT:
;      mbias     = Created master bias image
;      bihdr     = Header assosiated with the created master bias 
;      
; REQUIRES:
;     Hot and bad pixels have to be specified by NAN. Execute CONTROL_HOTBAD if needed.
;
; PROCEDURE:
;      Master bias creator for CUTE;
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;      modified 12.03.2019 by A. G. Sreejith
;      


pro control_bias_combine,bias_list,mbias,type,sat_value,threshold
print,'control bias'
err = ''
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    logprint,'CONTROL_BIAS_COMBINE: Syntax - control_bias_combine,bias_list,type',logonly = logonly
    message,'CONTROL_BIAS_COMBINE: Syntax - control_bias_combine,bias_list,type'
    err = '%control_bias_combine: Insufficient number of parameters'
    return
  endif
  
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Type input for bias combine not defined uisng default: median'
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
  
  n=n_elements(bias_list)
  if n eq 0 then begin
    logprint,'CONTROL_BIAS_COMBINE: Bias list is empty, exiting bias combine without creating master bias'
    return
  endif else begin
    gain_val=dblarr(n)
    bias=mrdfits(bias_list[0],0,hdr)
    nxy=size(bias)
    nx=nxy[1]
    ny=nxy[2]
    bias_ar=dblarr(nx,ny,n)
    include=intarr(n)
    
    totpix=n_elements(bias)
    limit=0.001*totpix
    for i=0, n-1 do begin
      filename=bias_list[i] ;modify based on file structure
      bias=mrdfits(filename,0,hdr)
      bias_nw=rejection(bias,threshold,npix)
      if npix lt limit then begin
        bias_ar[*,*,i]=bias_nw
        gain_val[i]=SXPAR( hdr, 'CCDGAIN')
        include[i]=i
      endif else include[i]=-1
    endfor
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    bias_arr=dblarr(nx,ny,n_frames)
    bias_arr=bias_arr(*,*,incl)
    case type of
      'median' : mbias_val = median(bias_arr,dimension=3,/even)
      'mean' : mbias_val = mean(bias_arr,dimension=3,/even,/NAN)
      'mode' : mbias_val = mode(bias_arr,dimension=3,/NAN)
    else : message,'Invalid type input for combine: Please recheck your input'
    endcase
    r_noise= stddev(mbias_val,/NAN)/mean(gain_val)
    ;if type = 'med' then mbias = median(bias_arr,dimension=3,/even,/NAN) 
    ;else if type = 'mea' then mbias =mean(bias_arr,dimension=3,/even,/NAN)
    ;else if type = 'mod' then mbias =median(bias_arr,dimension=3,/even,/NAN)
  
    t=SXPAR(hdr,'EXP_TIME')
    r_noise=SXPAR(hdr,'RD_NOISE') ;r_noise=SXPAR(hdr,'RNOISE')
    r_noise=12.25
    ;Header defnitions
    sxaddpar, bihdr, 'Time_in_JD', t
    sxaddpar, bihdr, 'RNOISE', r_noise, 'Readout noise'
    sxaddpar, bihdr, 'SIGMA',stddev(mbias_val,/NAN), 'Standard deviation of the frame'
    sxaddpar, bihdr, 'MEAN', mean(mbias_val,/NAN), 'Mean value of the frame'
    sxaddpar, bihdr, 'MEDIAN ', median(mbias_val), 'Median value of the frame'
    sxaddpar, bihdr, 'MAX', max(mbias_val,/NAN), 'Maximum value of the frame'
    sxaddpar, bihdr, 'MIN', min(mbias_val,/NAN), 'Minimum value of the frame'
    sxaddpar, bihdr, 'NFRAMES', n_frames, 'Number of frames used in bias combine'
    ;checks
    ;saturated pixels
    
    sat_loc = where(mbias_val ge sat_value)
    if total(sat_loc) eq -1 then sat_flag = 0 else sat_flag = 1
    ;deviation 
    std=stddev(mbias_val)
    std_loc = where((mbias_val ge (mean(mbias_val)+threshold*std)) or (mbias_val le (mean(mbias_val)-threshold*std)))
    if total(std_loc) eq -1 then std_flag = 0 else std_flag = 1
  
  
  
    mbias_flag = sat_flag+std_flag
    ;update flag headers
    sxaddpar, bihdr, 'SRNFLG', sat_flag
    sxaddpar, bihdr, 'STDFLG', std_flag
    sxaddpar, bihdr, 'MBIFLG', mbias_flag

    mbias={im:mbias_val,hdr:bihdr}

    return
  ;mwrfits,mbias,input_file+'m_bias.fits',bihdr,/create
    endelse
end