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


pro control_dark_combine,dark_list,type,mdark,mbias_str
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    print,'Syntax - control_dark_combine,dark_list,type,mbias_file,mdark
    err = '%control_dark_combine: Insufficient number of parameters'
    return
  endif
  
  dummy=mrdfits(dark_list[0],1,hdr)
  nxyd=size(dummy)
  ;read in master bias
  if keyword_defined(mbias) then begin
    mbias=mbias_str.im 
    mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
    nxyb=size(mbias)
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
  endelse
  n=n_elements(dark_list)
  dark_arr=dblarr(nxyd[1],nxyd[2],n)
  dark_err=dblarr(nxyd[1],nxyd[2],n)
  for i=0, n-1 do begin
    filename=dark_list[i]
    dark=mrdfits('filename',1,hdr)
    dark_arr[*,*,i]=dark-mbias
    dark_er[*,*,i]= dark+r^2
  endfor
  case type of
    'med' : begin
              mdark_val = median(dark_arr,dimension=3,/even,/NAN)
              sigma_db=total(dark_er,3,/even,/NAN); to modify
              mdark_err=sqrt(sigma_db/n^2)
            end
    'mea' : begin
              mdark_val = mean(dark_arr,dimension=3,/even,/NAN)
              sigma_db=total(dark_er,3,/even,/NAN)
              mdark_err=sqrt(sigma_db/n^2)
            end  
    'mod' : begin
              mdark_val = median(dark_arr,dimension=3,/even,/NAN)
              sigma_db=total(dark_er,3,/even,/NAN) ;to modify
              mdark_err=sqrt(sigma_db/n^2)
            end  
    else : print,'Invalid type input for combine: Please recheck your input'
  endcase
;  if type = 'med' then mdark = median(dark_arr,dimension=3,/even,/NAN)
;  else if type = 'mea' then mdark =mean(dark_arr,dimension=3,/even,/NAN)
;  else if type = 'mod' then mdark =median(dark_arr,dimension=3,/even,/NAN)
;  
  ;Header defnitions
  sxaddpar, dhdr, 'Time_in_JD', t
  sxaddpar, dhdr, 'RNOISE', r_noise, 'Readout noise'
  sxaddpar, dhdr, 'SIGMA',stddev(mdark_val,/NAN), 'Standard deviation of the frame'
  sxaddpar, dhdr, 'MEAN', mean(mdark_val,/NAN), 'Mean value of the frame'
  sxaddpar, dhdr, 'MEDIAN ', median(mdark_val,/NAN), 'Median value of the frame'
  sxaddpar, dhdr, 'MAX', max(mdark_val,/NAN), 'Maximum value of the frame'
  sxaddpar, dhdr, 'MIN', min(mdark_val,/NAN), 'Minimum value of the frame'
  sxaddpar, dhdr, 'EXPTIME', exptime, 'Exposure time of the frame'
  ;checks flags that ahve to be checked
  ;checks
  ;saturated pixels
  sat_value = 72000
  sat_loc = where(mdark_val ge sat_value)
  if size(sat_loc) eq 0 then sat_flag = 0 else sat_flag = 1
  ;deviation
  std=stddev(mdark)
  std_loc = where((mdark_val ge 5*std) or (mdark_val le 5*std))
  if size(std_loc) eq 0 then std_flag = 0 else std_flag = 1


  mdark_flag = sat_flag+std_flag+mbias_flag
  ;update flag headers
  sxaddpar, dhdr, 'CRFLG', cr_flag
  sxaddpar, dhdr, 'SRNFLG', sat_flag
  sxaddpar, dhdr, 'STDFLG', std_flag
  sxaddpar, dhdr, 'MBIFLG', mbias_flag
  sxaddpar, dhdr, 'MDRFLG', mdark_flag

   mdark={im:mdark_val,error:mdark_err,hdr:dhdr}

  return
  ;mwrfits,mdark,input_file+'m_bias.fits',dhdr,/create

end