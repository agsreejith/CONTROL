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

pro control_flat_combine,flat_list,type,mflat,mbias_str,mdark_str,sat_value,threshold
  if N_params() LT 4 then begin             ;Need at least 3 parameters
    print,'CONTROL_FLAT_COMBINE: Syntax - control_flat_combine,flat_list,type,mbias_file,mdark_file,mflat,fhdr
    return
  endif
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_TRACE: Type input for flat combine not defined uisng default: median'
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
  
  flat_dummy=mrdfits(flat_list[0],0,hdr)
  nxy=size(flat_dummy)
  totpix=n_elements(flat_dummy)
  limit=0.001*totpix
  ;read in master bias
  if keyword_defined(mbias_str) then begin
     mbias=mbias_str.im
     mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
     r=float(SXPAR( mbias_str.hdr, 'RNOISE'))
  endif else begin
    mbias=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mbias_flag = 0
    r=12.25
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
  if n ne 0 then begin
    flat_ar=dblarr(nxy[1],nxy[2],n)
    flat_er=dblarr(nxy[1],nxy[2],n)
    include=intarr(n)
    for i=0, n-1 do begin
      filename=flat_list[i]
      flat=mrdfits(filename,0,hdr)
      flat_type=SXPAR(hdr, 'FLATYPE')
      if (strpos(flat_type,'SKY') ge 0) then begin
        ;steps needed if the flat is an observation of scanned white dwarfs
      endif else begin
        flat_nw=rejection(flat,threshold,npix) ; interpolate pixels whose values have 5 sigma variations  
        if npix lt limit then begin
          flat_ar[*,*,i]=flat_nw-mdark-mbias
          rs_fb=sqrt(flat+r^2)
          flat_er[*,*,i]=rs_fb +(mdark_err/exp_dark)^2
          include[i]=i
        endif else include[i]=-1
      endelse  
    endfor
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    flat_arr=dblarr(nxy[1],nxy[2],n_frames)
    flat_err=dblarr(nxy[1],nxy[2],n_frames)
    flat_arr=flat_ar(*,*,incl)
    flat_err=flat_er(*,*,incl)
    case type of
      'median': begin
                 mflat_val = median(flat_arr,dimension=3,/even)
                 sigma_fdb=total(flat_err,3,/NAN)
                 mflat_err= sqrt(sigma_fdb/n^2)
                 nflat=mflat_val/mean(mflat_val)
                 sigma_fbdn=mflat_err/mean(mflat_val)
               end  
      'mean' : begin
                mflat_val = mean(flat_arr,dimension=3,/even,/NAN)
                sigma_fdb=total(flat_err,3,/NAN)
                mflat_err= sqrt(sigma_fdb/n^2)
                nflat=mflat_val/mean(mflat_val)
                sigma_fbdn=mflat_err/mean(mflat_val)
               end
      'mode' : begin
                mflat_val = median(flat_arr,dimension=3,/even,/NAN)
                sigma_fdb=total(flat_err,3,/NAN)
                mflat_err= sqrt(sigma_fdb/n^2)
                nflat=mflat_val/mean(mflat_val)
                sigma_fbdn=mflat_err/mean(mflat_val)
               end   
      else : message,'Invalid type input for combine: Please recheck your input'
   endcase
  
   error=sigma_fbdn
   nmflat=nflat
   ;  type eq 'med' then mflat = median(flat_arr,dimension=3,/even,/NAN)
   ;  else if type  eq 'mea' then mflat =mean(flat_arr,dimension=3,/even,/NAN)
   ;  else if type eq 'mod' then mflat =median(flat_arr,dimension=3,/even,/NAN)
   t=SXPAR(hdr,'EXP_TIME')

   r_noise=r
   ;Header defnitions
   sxaddpar, fhdr, 'Time_in_JD', t
   sxaddpar, fhdr, 'RNOISE', r_noise, 'Readout noise'
   sxaddpar, fhdr, 'SIGMA',stddev(nmflat,/NAN), 'Standard deviation of the frame'
   sxaddpar, fhdr, 'MEAN', mean(nmflat,/NAN), 'Mean value of the frame'
   sxaddpar, fhdr, 'MEDIAN ', median(nmflat), 'Median value of the frame'
   sxaddpar, fhdr, 'MAX', max(nmflat,/NAN), 'Maximum value of the frame'
   sxaddpar, fhdr, 'MIN', min(nmflat,/NAN), 'Minimum value of the frame'
   sxaddpar, fhdr, 'NFRAMES', n_frames, 'Number of frames used in flat combine'

   ;checks flags that have to be checked
   ;checks
   ;saturated pixels
   sat_value = 72000
   sat_loc = where(mflat_val ge sat_value)
   if total(sat_loc) ne -1 then sat_flag = 0 else sat_flag = 1
   ;deviation
   std=stddev(mflat_val)
   std_loc = where((mflat_val ge (mean(mflat_val)+threshold*std)) or (mflat_val le (mean(mflat_val)-threshold*std)))
   if total(std_loc) ne -1 then std_flag = 0 else std_flag = 1
   npix_stdloc=n_elements(std_loc)

   mflat_flag = sat_flag+std_flag+mbias_flag+mdark_flag
   ;update flag headers
   ;sxaddpar, fhdr, 'CRFLG', cr_flag
   sxaddpar, fhdr, 'SRNFLG', sat_flag
   sxaddpar, fhdr, 'STDFLG', std_flag
   sxaddpar, fhdr, 'MBIFLG', mbias_flag
   sxaddpar, fhdr, 'MDRFLG', mdark_flag
   sxaddpar, fhdr, 'MFLFLG', mflat_flag

   mflat={im:nmflat,error:error,hdr:fhdr}
   return
endif else begin
  logprint,'CONTROL_FLAT_COMBINE:FLAT list is empty, exiting flat combine without creating master flat'
  return
endelse
end