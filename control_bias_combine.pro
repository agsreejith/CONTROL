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
;      


pro control_bias_combine,bias_list,type,mbias
err = ''
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    print,'Syntax - control_bias_combine,bias_list,type
    err = '%control_bias_combine: Insufficient number of parameters'
    return
  endif
  n=n_elements(bias_list)
  gain_val=dblarr(n)
  for i=0, n-1 do begin
    filename=bias_list[i] ;modify based on file structure
    bias=mrdfits('filename',1,hdr)
    bias_arr[*,*,i]=bias
    gain_val[i]=SXPAR( hdr, 'CCDGAIN')
  endfor
  case type of
    'median' : mbias_val = median(bias_arr,dimension=3,/even,/NAN)
    'mean' : mbias_val = mean(bias_arr,dimension=3,/even,/NAN)
    'mode' : mbias_val = mode(bias_arr,dimension=3,/even,/NAN)
    else : print,'Invalid type input for combine: Please recheck your input'
  endcase
  r_noise= stddev(mbias_val,/NAN)/mean(gain_val)
  ;if type = 'med' then mbias = median(bias_arr,dimension=3,/even,/NAN) 
  ;else if type = 'mea' then mbias =mean(bias_arr,dimension=3,/even,/NAN)
  ;else if type = 'mod' then mbias =median(bias_arr,dimension=3,/even,/NAN)

  ;Header defnitions
  sxaddpar, bihdr, 'Time_in_JD', t
  sxaddpar, bihdr, 'RNOISE', r_noise, 'Readout noise'
  sxaddpar, bihdr, 'SIGMA',stddev(mbias_val,/NAN), 'Standard deviation of the frame'
  sxaddpar, bihdr, 'MEAN', mean(mbias_val,/NAN), 'Mean value of the frame'
  sxaddpar, bihdr, 'MEDIAN ', median(mbias_val,/NAN), 'Median value of the frame'
  sxaddpar, bihdr, 'MAX', max(mbias_val,/NAN), 'Maximum value of the frame'
  sxaddpar, bihdr, 'MIN', min(mbias_val,/NAN), 'Minimum value of the frame'
  
  ;checks
  ;saturated pixels
  sat_value = 72000
  sat_loc = where(mbias_val ge sat_value)
  if size(sat_loc) eq 0 then sat_flag = 0 else sat_flag = 1
  ;deviation 
  std=stddev(mbias_val)
  std_loc = where((mbias_val ge 5*std) or (mbias_val le 5*std))
  if size(std_loc) eq 0 then std_flag = 0 else std_flag = 1
  
  
  
  mbias_flag = sat_flag+std_flag
  ;update flag headers
  sxaddpar, bihdr, 'CRFLG', cr_flag
  sxaddpar, bihdr, 'SRNFLG', sat_flag
  sxaddpar, bihdr, 'STDFLG', std_flag
  sxaddpar, bihdr, 'MBIFLG', mbias_flag

mbias={im:mbias_val,hdr:bihdr}
  return
  ;mwrfits,mbias,input_file+'m_bias.fits',bihdr,/create

end