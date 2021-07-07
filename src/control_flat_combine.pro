; NAME:
;      CONTROL_FLAT_COMBINE
;
; PURPOSE:
;      Returns master flat frames from input flat files based on mean, 
;      median or mode combine as chosen by the user
;      
; CALLING SEQUENCE:
;      CONTROL_FLAT_COMBINE,flat_list,mflat,mbias_str,mdark_str,type=type,sat_value=sat_value$
;                           ,threshold=threshold
;
; INPUTS:
;      flat_list  = The file containing list of flat files.
;      mbias_str  = Master bias structure with the following structure 
;                   im:master bias data,hdr:master bias header,dq:master bias data quality.
;                   If this is undefined a master bias of zeros is created  
;      mdark_str  = Master dark structure with the following structure 
;                   im:master dark data,hdr:master dark header,dq:master dark data quality.
;                   If this is undefined a master dark of zeros is created
;
; OPTIONAL INPUTS
;      type      = Type of combining flat (mean, median, mode). Default is median.
;      sat_value = Expected saturation point for images. Default is 72000.
;      threshold = Value specifying acceptable deviation in the frame for data quality. 
;                  Default is 5.
;                  For example a threshold of 5 will flag all pixels whose value deviates from 
;                  mean of master flat by 5 sigma.
;
; OUTPUT:
;      mflat     = Created master flat image structure with image, header and data quality flag.

; PROCEDURE:
;      Master flat creator for CUTE;
;
;##################################################################################################

pro control_flat_combine,flat_list,mflat,mbias_str,mdark_str,type=type,sat_value=sat_value $
                        ,threshold=threshold
common pp_version                        
idl_ver=float(!Version.RELEASE)
  if N_params() LT 4 then begin             ;Need at least 3 parameters
    print,'CONTROL_FLAT_COMBINE: Syntax - control_flat_combine,flat_list,mflat,mbias_str,mdark_str
    return
  endif
  if keyword_defined(type) eq 0 then begin
    logprint,'CONTROL_FLAT_COMBINE: Type input for flat combine not defined uisng default: median'
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
  
  flat_dummy=mrdfits(flat_list[0],0,hdr,/SILENT)
  nxy=size(flat_dummy)
  totpix=n_elements(flat_dummy)
  dq_arr=bytarr(nxy[1],nxy[2])
  limit=0.001*totpix
  ;read in master bias
  if keyword_defined(mbias_str) then begin
    if (idl_ver ge 8) then begin
      if ISA(mbias_str) eq 0 then begin
        logprint,'CONTROL_FLAT_COMBINE: Requires a master bias file'
        logprint,'CONTROL_FLAT_COMBINE: Press any key other than q to create'$
                 +' a master bias frame with zeros or press q to quit.'
        R = GET_KBRD()
        if R eq 'q' then begin
          logprint,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.',logonly = logonly
          message,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.'
          err = '%control_flat_combine: Insufficient number of parameters'
          return
        endif else begin
          logprint,'CONTROL_FLAT_COMBINE: Creating a master bias with zero values.'
          mbias=make_array(nxyd[1],nxyd[2],/DOUBLE,value=0.0)
          mbias_flag = 0
          r=12.25
          mbdq=bytarr(nxyd[1],nxyd[2])
        endelse
      endif else begin
        mbias=mbias_str.im
        mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
        r=float(SXPAR( mbias_str.hdr, 'RNOISE'))
        mbdq=mbias_str.dq
      endelse 
    endif else begin
      mbas_def=datatype(mbias_str,2)
      if mbas_def eq 0 then begin
        logprint,'CONTROL_FLAT_COMBINE: Requires a master bias file'
        logprint,'CONTROL_FLAT_COMBINE: Press any key other thasn q to create'$
                 +' a master bias frame with zeros or press q to quit.'
        R = GET_KBRD()
        if R eq 'q' then begin
          logprint,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.',logonly = logonly
          message,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.'
          err = '%control_flat_combine: Insufficient number of parameters'
          return
        endif else begin
          logprint,'CONTROL_FLAT_COMBINE: Creating a master bias with zero values.'
          mbias=make_array(nxyd[1],nxyd[2],/DOUBLE,value=0.0)
          mbias_flag = 0
          r=12.25
          mbdq=bytarr(nxyd[1],nxyd[2])
        endelse
      endif else begin
        mbias=mbias_str.im
        mbias_flag = SXPAR( mbias_str.hdr, 'MBIFLG') 
        r=float(SXPAR( mbias_str.hdr, 'RNOISE'))
        mbdq=mbias_str.dq
      endelse 
    endelse  
  endif else begin
    logprint,'CONTROL_FLAT_COMBINE: Creating a master bias with zero values.'
    mbias=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mbias_flag = 0
    r=12.25
    mbdq=bytarr(nxy[1],nxy[2])
  endelse
  ;read in master dark
  if keyword_defined(mdark_str) then begin
    if (idl_ver ge 8) then begin
      if ISA(mdark_str) eq 0 then begin
        logprint,'CONTROL_FLAT_COMBINE: Requires a master dark file'
        logprint,'CONTROL_FLAT_COMBINE: Press any key other than q to create'$
                 +' a master dark frame with zeros or press q to quit.'
        R = GET_KBRD()
        if R eq 'q' then begin
          logprint,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.',logonly = logonly
          message,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.'
          err = '%control_flat_combine: Insufficient number of parameters'
          return
        endif else begin
          logprint,'CONTROL_FLAT_COMBINE: Creating a master dark with zero values.'
          mdark=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
          mdark_flag = 0
          mdark_err=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
          mddq=bytarr(nxy[1],nxy[2])
          exp_dark=300
        endelse
      endif else begin
        mdark=mdark_str.im
        mdark_flag = SXPAR( mdark_str.hdr, 'MDRFLG')
        mdark_err= mdark_str.error
        exp_dark=sxpar(mdark_str.hdr,'exptime')
        mddq=mdark_str.dq
      endelse
    endif else begin
      mbas_def=datatype(mdark_str,2)
      if mbas_def eq 0 then begin
        logprint,'CONTROL_FLAT_COMBINE: Requires a master dark file'
        logprint,'CONTROL_FLAT_COMBINE: Press any key other than q to create'$
                 +' a master dark frame with zeros or press q to quit.'
        R = GET_KBRD()
        if R eq 'q' then begin
          logprint,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.',logonly = logonly
          message,'CONTROL_FLAT_COMBINE: Exiting as requested by the user.'
          err = '%control_flat_combine: Insufficient number of parameters'
          return
        endif else begin
          logprint,'CONTROL_FLAT_COMBINE: Creating a master dark with zero values.'
          mdark=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
          mdark_flag = 0
          mdark_err=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
          mddq=bytarr(nxy[1],nxy[2])
          exp_dark=300
        endelse
      endif else begin
        mdark=mdark_str.im
        mdark_flag = SXPAR( mdark_str.hdr, 'MDRFLG')
        mdark_err= mdark_str.error
        exp_dark=sxpar(mdark_str.hdr,'EXPTIME')
        mddq=mdark_str.dq
      endelse
    endelse  
  endif else begin
    logprint,'CONTROL_FLAT_COMBINE: Creating a master dark with zero values.'
    mdark=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mdark_flag = 0
    mdark_err=make_array(nxy[1],nxy[2],/DOUBLE,value=0.0)
    mddq=bytarr(nxy[1],nxy[2])
    exp_dark=300
  endelse
 
  n=n_elements(flat_list)
  if n ne 0 then begin
    flat_ar=dblarr(nxy[1],nxy[2],n)
    flat_er=dblarr(nxy[1],nxy[2],n)
    include=intarr(n)
    for i=0, n-1 do begin
      filename=flat_list[i]
      flat=mrdfits(filename,0,hdr,/SILENT)
      flat_type=SXPAR(hdr, 'FLATYPE')
      ccd_gain=SXPAR( hdr, 'CCDGAIN')
      if ccd_gain le 0 then ccd_gain=1
      if (strpos(flat_type,'SKY') ge 0) then begin
        ;steps needed if the flat is an observation of scanned white dwarfs
      endif else begin
        FITS_INFO, filename,N_ext =numext,/SILENT
        if numext gt 0 then dq=mrdfits(filename,1,hdr1,/SILENT) else dq=bytarr(nxy[1],nxy[2])
        flat_nw=rejection(flat,threshold,npix) ; interpolate pixels that have 5 sigma variations  
        if npix lt limit then begin
          flat_ar[*,*,i]=((flat_nw)-mbias)-mdark
          fmb=(flat_nw)-mbias
          rs_fbe = (fmb/ccd_gain))
          negative=where(rs_fbe lt 0)
          rs_fbe[negative]=r^2
          rs_fb=sqrt(rs_fbe)
          flat_er[*,*,i]=rs_fb +(mdark_err/exp_dark)^2
          include[i]=i
        endif else include[i]=-1
      endelse  
    dq_arr = dq_arr or dq
    endfor
    incl=where(include ge 0)
    n_frames=n_elements(incl)
    flat_arr=dblarr(nxy[1],nxy[2],n_frames)
    flat_err=dblarr(nxy[1],nxy[2],n_frames)
    flat_arr=flat_ar(*,*,incl)
    flat_err=flat_er(*,*,incl)
    if n_frames eq 0 then begin
      logprint,'CONTROL_FLAT_COMBINE: No valid FLAT file found. Terminating MASTER FLAT creation.'
      return
    endif
    if n_frames le 1 then begin
      logprint,'CONTROL_FLAT_COMBINE: Only one valid FLAT file found.'$
               +' Do you wnat to assume it as the MASTER FLAT?.'
      logprint,'Press q to skip this assumption.'$
               +' Press any key to continue with MASTER FLAT creation with one valid FLAT file.'
      Rd = GET_KBRD()
      if Rd eq 'q' then begin
        logprint,'CONTROL_FLAT_COMBINE: Terminating MASTER FLAT creation as requested by the user.'
        return
      endif
    endif
    logprint,'CONTROL_FLAT_COMBINE: Combining '+STRTRIM(STRING(n_frames),2)$
             +' FLAT file to create MASTER FLAT using '+type+' method .'
    
    case type of
      'median': begin
                 if (idl_ver gt 5.6) then mflat_val = median(flat_arr,dimension=3,/even) else begin
                   totalImage = Total(flat_arr, 3)
                   minValue = Min(flat_arr, Dimension=3)
                   maxValue = Max(flat_arr, Dimension=3)
                   mflat_val = totalImage - minValue - maxValue
                 endelse
                 sigma_fdb=total(flat_err,3,/NAN)
                 mflat_err= sqrt(sigma_fdb/n_frames^2)
                 nflat=mflat_val/mean(mflat_val,/NAN)
                 sigma_fbdn=mflat_err/mean(mflat_val,/NAN)
               end  
      'mean' : begin
                if (idl_ver gt 8.0) then mflat_val = mean(flat_arr,dimension=3,/NAN) else begin
                  totalImage = Total(flat_arr, 3)
                  mflat_val = totalImage/n_frames
                endelse
                sigma_fdb=total(flat_err,3,/NAN)
                mflat_err= sqrt(sigma_fdb/n_frames^2)
                nflat=mflat_val/mean(mflat_val,/NAN)
                sigma_fbdn=mflat_err/mean(mflat_val,/NAN)
               end
      'mode' : begin
                mflat_val = mode(flat_arr,dimension=3,/NAN)
                sigma_fdb=total(flat_err,3,/NAN)
                mflat_err= sqrt(sigma_fdb/n_frames^2)
                nflat=mflat_val/mean(mflat_val,/NAN)
                sigma_fbdn=mflat_err/mean(mflat_val,/NAN)
               end   
      else : begin
              errorlog,'Invalid type input for combine: Please recheck your input',logonly=1
              message,'Invalid type input for combine: Please recheck your input'
           end  
   endcase
  
   error=sigma_fbdn
   nmflat=nflat
   ;  type eq 'med' then mflat = median(flat_arr,dimension=3,/even,/NAN)
   ;  else if type  eq 'mea' then mflat =mean(flat_arr,dimension=3,/even,/NAN)
   ;  else if type eq 'mod' then mflat =median(flat_arr,dimension=3,/even,/NAN)
   r_noise=SXPAR(hdr,'RNOISE')
   exptime=SXPAR(hdr,'EXPTIME')

   r_noise=r
   ;Header defnitions
   sxaddpar, fhdr, 'TELESCOP', SXPAR(hdr,'TELESCOP'),'Telescope name'
   sxaddpar, fhdr, 'ROOTNAME', SXPAR(hdr,'ROOTNAME'),'Root directory'
   sxaddpar, fhdr, 'FILENAME', SXPAR(hdr,'FILENAME'),'Filename'
   sxaddpar, fhdr, 'PRGRM_ID', SXPAR(hdr,'PRGRM_ID'),'Program ID'
   sxaddpar, fhdr, 'TARGT_ID', SXPAR(hdr,'TARGT_ID'),'Target ID'
   sxaddpar, fhdr, 'EXP_ID  ', SXPAR(hdr,'EXP_ID'),'Exposure ID'
   sxaddpar, fhdr, 'OBS_ID  ', SXPAR(hdr,'OBS_ID'),'Observation ID'
   sxaddpar, fhdr, 'FILETYPE', 'MFLAT','Type of observation'
   sxaddpar, fhdr, 'RNOISE  ', r_noise, 'Readout noise'
   sxaddpar, fhdr, 'SIGMFLAT', stddev(mflat_val,/NAN), 'Standard deviation of the frame'
   sxaddpar, fhdr, 'MEANFLAT', mean(mflat_val,/NAN), 'Mean value of the frame'
   sxaddpar, fhdr, 'MDNFLAT ', median(mflat_val), 'Median value of the frame'
   sxaddpar, fhdr, 'MAXFLAT ', max(mflat_val,/NAN), 'Maximum value of the frame'
   sxaddpar, fhdr, 'MINFLAT ', min(mflat_val,/NAN), 'Minimum value of the frame'
   sxaddpar, fhdr, 'EXPTIME', exptime, 'Exposure time of the single frame'
   sxaddpar, fhdr, 'NFRAMES', n_frames, 'Number of frames used in flat combine'
   sxaddpar, fhdr, 'FLATTYP ',type,'Type of flat combine employed'
   sxaddpar, fhdr, 'FLATSAT ',sat_value,'Saturation limit used in flat frames'
   sxaddpar, fhdr, 'FLATSIG ',threshold,'Deviation limit for good flat frames'
   sxaddpar, fhdr, 'CCDGAIN ',SXPAR(hdr,'CCDGAIN'),'CCD gain'
   sxaddpar, fhdr, 'CCDTEMP ',SXPAR(hdr,'CCDTEMP'),'CCD temperature'
   sxaddpar, fhdr, 'TECBTEM ',SXPAR(hdr,'TECBTEM'),'TEC backside temperature'
   sxaddpar, fhdr, 'RADTEMP ',SXPAR(hdr,'RADTEMP'),'Radiator temperature'
   sxaddpar, fhdr, 'SHTRSTS ',SXPAR(hdr,'SHTRSTS'),'Shutter status'
   sxaddpar, fhdr, 'PIPENUM',version,'Pipeline version number used to reduce the data'
   sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V2.0 .',hdr,/COMMENT
   ;checks flags that have to be checked
   dq_arr=dq_arr or mbdq
   dq_arr=dq_arr or mddq
   prb=where(dq_arr ge 1)
   flat_dq=bytarr(nxy[1],nxy[2])
   if total(prb) ne -1 then flat_dq[prb]=dq_arr[prb] or 128b
   ;checks
   ;saturated pixels
   
   sat_loc = where(mflat_val ge sat_value)
   if total(sat_loc) ne -1 then sat_flag = 0 else sat_flag = 1
   if total(sat_loc) ne -1 then flat_dq[sat_loc]=flat_dq[sat_loc] or 8b
   ;deviation
   std=stddev(mflat_val)
   std_loc = where((mflat_val ge (mean(mflat_val)+threshold*std)) or $
                  (mflat_val le (mean(mflat_val)-threshold*std)))
   if total(std_loc) ne -1 then std_flag = 0 else std_flag = 1
   npix_stdloc=n_elements(std_loc)
   if total(std_loc) ne -1 then flat_dq[std_loc]=flat_dq[std_loc] or 8b
   nan_loc = where(finite(mflat_val, /NAN) eq 1)
   if total(nan_loc) eq -1 then nan_flag = 0 else nan_flag = 1
   if total(nan_loc) ne -1 then flat_dq[nan_loc]=flat_dq[nan_loc] or 8b

   mflat_flag = sat_flag+std_flag+mbias_flag+nan_flag+mdark_flag
   if mflat_flag gt 1 then mflat_flag=1
 
   ;update flag headers
   ;sxaddpar, fhdr, 'CRFLG', cr_flag
   sxaddpar, fhdr, 'SRNFLG', sat_flag
   sxaddpar, fhdr, 'STDFLG', std_flag
   sxaddpar, fhdr, 'MBIFLG', mbias_flag
   sxaddpar, fhdr, 'MDRFLG', mdark_flag
   sxaddpar, fhdr, 'MFLFLG', mflat_flag
   logprint,'CONTROL FLAT COMBINE: Saturated (values above '+STRTRIM(STRING(sat_value),2)$
            +') pixels and pixels that deviate by '+STRTRIM(STRING(threshold),2)$
            +' sigma from mean of  MASTER FLAT have been flaged.'
   mflat={im:nmflat,error:error,hdr:fhdr,dq:flat_dq}
 
   return
endif else begin
  logprint,'CONTROL_FLAT_COMBINE:FLAT list is empty,'$
           +' exiting flat combine without creating master flat'
  return
endelse
end