; NAME:
;      CONTROL_FLAT_COMBINE
;
; PURPOSE:
;      Returns master flat frames from input bias files based on mean,median or mode combine as chosen by the user;
; CALLING SEQUENCE:
;      CONTROL_FLAT_COMBINE,flat_list,type,mbias_file,mdark_file,mflat,fhdr
;
; INPUTS:
;      flat_list  = The file containing list of flat files.
;      type       = Type of combining flats (mea = mean, med = median, mod = mode)
; OPTIONAL INPUTS:     
;      mbias_file = Location of master bias file
;      mdark_file = Location of master dark file
; 
; OUTPUT:
;      mflat      = Created master flat image
;      fhdr       = header file for master flat 
; REQUIRES:
;     Hot and bad pixels have to be specified by NAN execute CONTROL_HOTBAD if needed.
;
; PROCEDURE:
;      Master flat creator for CUTE;
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;      modified 14.01.2018 by A. G. Sreejith
;#######################################################################

pro control_flat_combine,flat_list,type,nmflat,mbias_str,mdark_str
  if N_params() LT 4 then begin             ;Need at least 3 parameters
    print,'Syntax - control_flat_combine,flat_list,type,mbias_file,mdark_file,mflat,fhdr
    return
  endif

  flat_dummy=mrdfits(flat_list[0],1,hdr)
  nxy=size(flat_dummy)
  ;read in master bias
  if keyword_defined(mbias_str) then begin
    mbias=mbias_str.im
     mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
  endif else begin
    mbias=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mbias_flag = 0
  endelse
  ;read in master dark
  if keyword_defined(mdark_str) then begin
    mdark=mdark_str.im
    mdark_flag = SXPAR( mdark_str.hdr, 'MDRFLG')
    mdark_err= mdark_str.error
    exp_dark=sxpar(mdark_str.hdr,'exptime')
  endif else begin
    mdark=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mdark_flag = 0
    mdark_err=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
  endelse
 
  n=n_elements(flat_list)
  flat_arr=dblarr(nxy[1],nxy[2],n)
  flat_err=dblarr(nxy[1],nxy[2],n)
  for i=0, n-1 do begin
    filename=flat_list[i]
    flat=mrdfits('filename',1,hdr)
    if (strpos(flat_type,'SKY') ge 0) then begin
      
    endif else begin
      flat_arr[*,*,i]=flat-mdark-mbias
      rs_fb=sqrt(flat+r^2)
      flat_err[*,*,i]=rs_fb +(mdark_err/exp_dark)^2
    endelse
  endfor
  case type of
    'med' : begin
              mflat = median(flat_arr,dimension=3,/even,/NAN)
              sigma_fdb=total(flat_err,3,/even,/NAN)
              mflat_err= sqrt(sigma_fdb/n^2)
              nflat=mflat_val/mean(mflat_val)
              sigma_fbdn=mflat_err/mean(mflat_val)
            end  
    'mea' : begin
              mflat_val = mean(flat_arr,dimension=3,/even,/NAN)
              sigma_fdb=total(flat_err,3,/even,/NAN)
              mflat_err= sqrt(sigma_fdb/n^2)
              nflat=mflat_val/mean(mflat_val)
              sigma_fbdn=mflat_err/mean(mflat_val)
            end
    'mod' : begin
              mflat_val = median(flat_arr,dimension=3,/even,/NAN)
              sigma_fdb=total(flat_err,3,/even,/NAN)
              mflat_err= sqrt(sigma_fdb/n^2)
              nflat=mflat_val/mean(mflat_val)
              sigma_fbdn=mflat_err/mean(mflat_val)
            end   
    else : print,'Invalid type input for combine: Please recheck your input'
  endcase
  
  error=sigma_fbdn
  nmflat=nflat
;  type eq 'med' then mflat = median(flat_arr,dimension=3,/even,/NAN)
;  else if type  eq 'mea' then mflat =mean(flat_arr,dimension=3,/even,/NAN)
;  else if type eq 'mod' then mflat =median(flat_arr,dimension=3,/even,/NAN)

 ;Header defnitions
 sxaddpar, fhdr, 'Time_in_JD', t
 sxaddpar, fhdr, 'RNOISE', r_noise, 'Readout noise'
 sxaddpar, fhdr, 'SIGMA',stddev(nmflat,/NAN), 'Standard deviation of the frame'
 sxaddpar, fhdr, 'MEAN', mean(nmflat,/NAN), 'Mean value of the frame'
 sxaddpar, fhdr, 'MEDIAN ', median(nmflat,/NAN), 'Median value of the frame'
 sxaddpar, fhdr, 'MAX', max(nmflat,/NAN), 'Maximum value of the frame'
 sxaddpar, fhdr, 'MIN', min(nmflat,/NAN), 'Minimum value of the frame'

 ;checks flags that have to be checked
 ;checks
 ;saturated pixels
 sat_value = 72000
 sat_loc = where(mflat_val ge sat_value)
 if size(sat_loc) eq 0 then sat_flag = 0 else sat_flag = 1
 ;deviation
 std=stddev(nmflat)
 std_loc = where((nmflat ge 5*std) or (nmflat le 5*std))
 if size(std_loc) eq 0 then std_flag = 0 else std_flag = 1


 mflat_flag = sat_flag+std_flag+mbias_flag+mdark_flag
 ;update flag headers
 sxaddpar, fhdr, 'CRFLG', cr_flag
 sxaddpar, fhdr, 'SRNFLG', sat_flag
 sxaddpar, fhdr, 'STDFLG', std_flag
 sxaddpar, fhdr, 'MBIFLG', mbias_flag
 sxaddpar, fhdr, 'MDRFLG', mdark_flag
 sxaddpar, fhdr, 'MFLFLG', mflat_flag

mflat={im:nmflat,error:error,hdr:fhdr}
 return

end