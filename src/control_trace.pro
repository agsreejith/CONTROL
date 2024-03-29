; NAME:
;      CONTROL_TRACE
;
; PURPOSE:
;      Returns extracted spectrum and assosiated products.
;      
; CALLING SEQUENCE:
;      spectrum = CONTROL_TRACE(in_image,infile,trace_type,filename)
;
; INPUTS:
;      in_image = The 2D array from which spectrum is to extracted.
;      infile   = Parameter file for CONTROL
;      filename = Filename coresponding to the in_image 
;
; OPTIONAL INPUTS
;      trace_type = Type of spectrum extraction(simple,fixed,variable,function). Default is simple.
;      
; OUTPUT:
;      spectrum     = Structure containing extracted 1D spectrum and 1D background with assosiated
;                     error, data quality flags and header.
;
; PROCEDURE:
;      Spectrum trace and extractor for CUTE;
;
;##################################################################################################

function arrayor,array,key
;function to OR of elements across array
;INPUTS:
;     array = input array
;     key   = determines the direction for OR operation, 1 across row, 2 across column, 3 across row and column
;OUTPUT:
;     returns the OR ed array/value     
;#########################################
  axy=size(array)
  ax=axy[1]
  ay=axy[2]
  if key eq 1 then begin
    or_val=dblarr(ay)
    for aj=0,ay-1 do begin
      for ai =0, ax-1 do begin
        or_val[aj]=or_val[aj] or array[ai,aj]
      endfor
    endfor
  endif
  if key eq 2 then begin
    or_val=dblarr(ax)
    for ai=0,ax-1 do begin
      for aj =0, ay-1 do begin
        or_val[ai]=or_val[ai] or array[ai,aj]
      endfor
    endfor
  endif
  if key eq 2 then begin
    or_val=0
    for ai=0,ax-1 do begin
      for aj =0, ay-1 do begin
        or_val=or_val or array[ai,aj]
      endfor
    endfor
  endif
  return, or_val
end

;Simplest extraction a box with fixed width, centroid and slope
function extract_box,im,sigma,dq,centroid,slope,width,back_trace
  err = ''
  if N_params() LT 5 then begin             ;Need at least 4 parameters
    message,'Syntax - extract_box,im,sigma,slope,width,centroid
    err = '%extract_box: Insufficient number of parameters'

  endif
  nxy=size(im)
  if nxy[0] ne 2 then begin
    message,'Input file error, the image file is not 2-D'
    err = '%extract_box: Input file error'

  endif
  nxy2=size(sigma)
  if nxy2[0] ne 2 then begin
    message,'Input file error, the error file is not 2-D'
    err = '%extract_box: Input file error'

  endif
  nx=nxy[1]
  ny=nxy[2]
  multiple=back_trace.type
  if (datatype(back_trace.shift_upper) eq 'INT') then ushift=fix(back_trace.shift_upper)
  if (datatype(back_trace.shift_lower) eq 'INT') then lshift=fix(back_trace.shift_lower)
  if (datatype(back_trace.shift_up_width) eq 'INT') then begin
    uwidth=fix(back_trace.shift_up_width)
    if (total(uwidth) eq n_elements(uwidth)) then begin
      single_upix=1
      ubg_length=double(n_elements(uwidth))
    endif else begin
      single_upix=0
      ubg_length=double(2*total(uwidth)+n_elements(uwidth))
    endelse
  endif
  if (datatype(back_trace.shift_low_width) eq 'INT') then begin
    lwidth=fix(back_trace.shift_low_width)
    if (total(lwidth) eq n_elements(lwidth)) then begin
      single_lpix=1
      lbg_length=double(n_elements(lwidth))
    endif else begin
      single_lpix=0
      lbg_length=double(2*total(lwidth)+n_elements(lwidth))
    endelse
  endif  
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  dq_val=bytarr(nx)
  dq_bval=bytarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  pixel=indgen(nx)
  if keyword_defined(centroid) then centroid=centroid else begin
    logprint,'CONTROL_TRACE: centroid for box extraction not specified. Calculating on the go'$
            ,logonly=logonly
    colsum=total(im,1)
    maximum=max(colsum,loc)
    centroid=make_array(nx,value=loc)
  endelse

  if width le 0 then begin
    logprint,'Width specified is less than or equal to zero.'$
            +' CONTROL will reset it to default value of 10.'
    width =10
  endif
  stp= fix(centroid+slope*pixel+width)
  str=fix(centroid+slope*pixel-width)
  centroidn=fix(centroid+slope*pixel)
  centroid=centroidn
  out_frame_up=where(stp gt ny-1)
  out_frame_down=where(str lt 0)
  if total(out_frame_up) ge 0 or total(out_frame_down) ge 0 then begin
    errorlog,'CONTROL_TRACE: extraction width specified is outside the CCD area.'$
      +' Please re-enter a new width for extraction',logonly=logonly
    message,' CONTROL_TRACE: extraction width specified is outside the CCD area.'$
      +'Please re-enter a new width for extraction'
  endif else begin
    for j=0,nx-1 do begin
      spectrum_val[j]=total(im[j,str[j]:stp[j]], /NAN)
      noise[j]=total(sigma[j,str[j]:stp[j]]^2, /NAN)
      dq_val[j]= dq_val[j] or arrayor(dq[j,str[j]:stp[j]],2)
      ;dq_val[j]=total(dq[j,str[j]:stp[j]], /NAN)
      ;for background use the input from lookup table as well
      if multiple eq 0 then begin
        if ushift ne 0 then begin
          if ((stp[j]+ushift) le ny-1) then begin
            back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift], /NAN)
            noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2, /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,str[j]+ushift:stp[j]+ushift], 2)
            ;dq_bval[j]+=total(dq[j,str[j]+ushift:stp[j]+ushift], /NAN)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new upper offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                    +' Please re-enter a new upper offset for backgrounds'
          endelse
          background[j]=back_u[j]
          backg_error[j]=sqrt(noise_u[j])
        endif
        if lshift ne 0 then begin
          if ((str[j]-lshift) ge 0) then begin
            back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift], /NAN)
            noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2, /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,str[j]-lshift:stp[j]-lshift],2)
            ;dq_bval[j]+=total(dq[j,str[j]-lshift:stp[j]-lshift], /NAN)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new lower offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new lower offset for backgrounds'
          endelse
          background[j]=back_l[j]
          backg_error[j]=sqrt(noise_l[j])
        endif
      endif else begin
        for k=0,n_elements(ushift)-1 do begin
          if single_upix eq 0 then begin
            ubgstr=centroidn[j]+ushift[k]-uwidth[k]
            ubgstop=centroidn[j]+ushift[k]+uwidth[k]
          endif else begin
            ubgstr=centroidn[j]+ushift[k]
            ubgstop=centroidn[j]+ushift[k]
          endelse  
          if ((ubgstop) le ny-2) then begin
            back_u[j]=back_u[j]+total(im[j,ubgstr:ubgstop], /NAN)
            noise_u[j]=noise_u[j]+total(sigma[j,ubgstr:ubgstop]^2, /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,ubgstr:ubgstop],2)
            ;dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new upper offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new upper offset for backgrounds'
          endelse
        endfor
        for k=0,n_elements(lshift)-1 do begin
          if single_lpix eq 0 then begin
            lbgstr=centroidn[j]+lshift[k]-lwidth[k]
            lbgstop=centroidn[j]+lshift[k]+lwidth[k]
          endif else begin
            lbgstr=centroidn[j]+lshift[k]
            lbgstop=centroidn[j]+lshift[k]+1
          endelse
          if ((lbgstr) ge 1) then begin
            back_l[j]=back_l[j]+total(im[j,lbgstr:lbgstop], /NAN)
            noise_l[j]=noise_l[j]+total(sigma[j,lbgstr:lbgstop]^2, /NAN)
            dq_bval[j]+=arrayor(dq[j,lbgstr:lbgstop],2)
            ;dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new lower offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new lower offset for backgrounds'
          endelse
        endfor
        bg_length=double(ubg_length+lbg_length)
        spec_length=double(stp[j]-str[j])
        bgscale=float(bg_length/spec_length)
        background[j]=(back_u[j]+back_l[j])/bgscale
        backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
      endelse
    endfor
  endelse
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=dq_val(prb) or 128b
  undefine,prb
  prbb=where(dq_bval ge 1)
   if total(prbb) ne -1 then dq_bval(prbb)= dq_bval(prbb) or 128b
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,$
            lower:str,upper:stp,dq:dq_val,dq_bcg:dq_bval}
  return,spectrum
end

;extraction with pixel dependednd centroid and fixed width
function extract_sum,im,sigma,dq,centroid,slope,upper,lower,back_trace
  err = ''
  if N_params() LT 5 then begin             ;Need at least 4 parameters
    message,'Syntax - extract_sum,im,sigma,centroid,upper,lower
    err = '%extract_sum: Insufficient number of parameters'

  endif
  nxy=size(im)
  if nxy[0] ne 2 then begin
    message,'Input file error, the image file is not 2-D'
    err = '%extract_sum: Input file error'

  endif
  nxy2=size(sigma)
  if nxy2[0] ne 2 then begin
    print,'Input file error, the error file is not 2-D'
    err = '%extract_sum: Input file error'

  endif
  nx=nxy[1]
  ny=nxy[2]
  multiple=back_trace.type
  if (datatype(back_trace.shift_upper) eq 'INT') then ushift=fix(back_trace.shift_upper)
  if (datatype(back_trace.shift_lower) eq 'INT') then lshift=fix(back_trace.shift_lower)
  if (datatype(back_trace.shift_up_width) eq 'INT') then begin
    uwidth=fix(back_trace.shift_up_width)
    if (total(uwidth) eq n_elements(uwidth)) then begin
      single_upix=1
      ubg_length=double(n_elements(uwidth))
    endif else begin
      single_upix=0
      ubg_length=double(2*total(uwidth)+n_elements(uwidth))
    endelse
  endif
  if (datatype(back_trace.shift_low_width) eq 'INT') then begin
    lwidth=fix(back_trace.shift_low_width)
    if (total(lwidth) eq n_elements(width)) then begin
      single_lpix=1
      lbg_length=double(n_elements(width))
    endif else begin
      single_lpix=0
      lbg_length=double(2*total(width)+n_elements(lwidth))
    endelse
  endif  
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  dq_val=bytarr(nx)
  dq_bval=bytarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  pixel=indgen(nx)
  if n_elements(centroid) ne nx then begin
    message,'Input dimension error, centroid has different dimension than image array '
    err = '%extract_sum: Input diamension error'
  endif
  if upper le 0 then begin
    logprint,'Upper width specified is less than or equal to zero.'$
            +' CONTROL will reset it to default value of 10.'
    upper =10
  endif
  if lower le 0 then begin
    logprint,'Lower width specified is less than or equal to zero.'$
            +' CONTROL will reset it to default value of 10.'
    lower =10
  endif
  stp= fix(centroid+slope*pixel+upper)
  str=fix(centroid+slope*pixel-lower)
  centroid=fix(centroid+slope*pixel)
  out_frame_up=where(stp gt ny-1)
  out_frame_down=where(str lt 0)
  if total(out_frame_up) ge 0 or total(out_frame_down) ge 0 then begin
    errorlog,'CONTROL_TRACE: extraction width specified is outside the CCD area.'$
      +' Please re-enter a new width for extraction',logonly=logonly
    message,' CONTROL_TRACE: extraction width specified is outside the CCD area.'$
      +'Please re-enter a new width for extraction'
  endif else begin
    for j=0,nx-1 do begin
      spectrum_val[j]=total(im[j,str[j]:stp[j]], /NAN)
      noise[j]=total(sigma[j,str[j]:stp[j]]^2, /NAN)
      ;dq_val[j]=total(dq[j,str[j]:stp[j]], /NAN)
      dq_val[j]=dq_val[j] or arrayor(dq[j,str[j]:stp[j]],2)
      ;for background use the input from lookup table as well
      if multiple eq 0 then begin
        if ushift ne 0 then begin
          if ((stp[j]+ushift) le ny-1) then begin
            back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift], /NAN)
            noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2, /NAN)
            ;dq_bval[j]+=total(dq[j,str[j]+ushift:stp[j]+ushift], /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,str[j]+ushift:stp[j]+ushift], 2)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area'$
                    ,logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new upper offset for backgrounds'
          endelse
          background[j]=back_u[j]
          backg_error[j]=sqrt(noise_u[j])
        endif
        if lshift ne 0 then begin  
            if ((str[j]-lshift) ge 0) then begin
              back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift], /NAN)
              noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2, /NAN)
              ;dq_bval[j]+=total(dq[j,str[j]-lshift:stp[j]-lshift], /NAN)
              dq_bval[j]=dq_bval[j] or arrayor(dq[j,str[j]-lshift:stp[j]-lshift], 2)
            endif else begin
              errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                      +' Please re-enter a new lower offset for backgrounds',logonly=logonly
              message,' Background offset specified is outside the CCD area.'$
                     +' Please re-enter a new lower offset for backgrounds'
            endelse
            background[j]=(back_l[j])
            backg_error[j]=sqrt(noise_l[j])
        endif
      endif else begin
        for k=0,n_elements(ushift)-1 do begin
          if single_upix eq 0 then begin
            ubgstr=centroid[j]+ushift[k]-uwidth[k]
            ubgstop=centroid[j]+ushift[k]+uwidth[k]
          endif else begin
            ubgstr=centroid[j]+ushift[k]
            ubgstop=centroid[j]+ushift[k]+1
          endelse
          if ((ubgstop) le ny-2) then begin
            back_u[j]=total(im[j,ubgstr:ubgstop], /NAN)
            noise_u[j]=total(sigma[j,ubgstr:ubgstop], /NAN)
            ;dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,ubgstr:ubgstop],2)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new upper offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new upper offset for backgrounds'
          endelse
        endfor
        for k=0,n_elements(lshift)-1 do begin
          if single_lpix eq 0 then begin
            lbgstr=centroid[j]+lshift[k]-lwidth[k]
            lbgstop=centroid[j]+lshift[k]+lwidth[k]
          endif else begin
            lbgstr=centroid[j]+lshift[k]
            lbgstop=centroid[j]+lshift[k]+1
          endelse
          if ((lbgstr) ge 1) then begin
            back_u[j]=total(im[j,lbgstr:lbgstop], /NAN)
            noise_u[j]=total(sigma[j,lbgstr:lbgstop], /NAN)
            dq_bval[j]=dq_bval[j] or arrayor(dq[j,lbgstr:lbgstop],2)
            ;dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          endif else begin
            errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                    +' Please re-enter a new upper offset for backgrounds',logonly=logonly
            message,' Background offset specified is outside the CCD area.'$
                   +' Please re-enter a new lower offset for backgrounds'
          endelse
        endfor
        bg_length=double(ubg_length+lbg_length)
        spec_length=double(stp[j]-str[j])
        bgscale=float(bg_length/spec_length)
        background[j]=(back_u[j]+back_l[j])/bgscale
        backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
      endelse
    endfor
  endelse
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=dq_val(prb) or 128b
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=dq_bval(prbb) or 128b
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,$
            lower:str,upper:stp,dq:dq_val,dq_bcg:dq_bval}
  return,spectrum
end

;extraction with pixel dependend maxvalue and pixel dependend width

function extract_varsum,im,sigma,dq,centroid,threshold,back_trace
  err = ''
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    message,'Syntax - extract_varsum,im,sigma,centroid,threshold
    err = '%extract_varsum: Insufficient number of parameters'

  endif
  nxy=size(im)
  if nxy[0] ne 2 then begin
    message,'Input file error, the image file is not 2-D'
    err = '%extract_varsum: Input file error'

  endif
  nxy2=size(sigma)
  if nxy2[0] ne 2 then begin
    message,'Input file error, the error file is not 2-D'
    err = '%extract_varsum: Input file error'

  endif

  nx=nxy[1]
  ny=nxy[2]
  multiple=back_trace.type

  if (datatype(back_trace.shift_upper) eq 'INT') then ushift=fix(back_trace.shift_upper)
  if (datatype(back_trace.shift_lower) eq 'INT') then lshift=fix(back_trace.shift_lower)
  if (datatype(back_trace.shift_up_width) eq 'INT') then begin
    uwidth=fix(back_trace.shift_up_width)
    if (total(uwidth) eq n_elements(uwidth)) then begin
      single_upix=1
      ubg_length=double(n_elements(uwidth))
    endif else begin
      single_upix=0
      ubg_length=double(2*total(uwidth)+n_elements(uwidth))
    endelse
  endif
  if (datatype(back_trace.shift_low_width) eq 'INT') then begin
    lwidth=fix(back_trace.shift_low_width)
    if (total(lwidth) eq n_elements(width)) then begin
      single_lpix=1
      lbg_length=double(n_elements(width))
    endif else begin
      single_lpix=0
      lbg_length=double(2*total(width)+n_elements(lwidth))
    endelse
  endif  
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  dq_val=bytarr(nx)
  dq_bval=bytarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  upper=intarr(nx)
  lower=intarr(nx)
  centroid=intarr(nx)
  m=max(im,loc,dimension=2)
  xloc=loc mod nx
  yloc=loc/nx
  centroid=yloc
  for j=1,nx-2 do begin ;valid for simulated data change to: j=0,nx-1 for actual data
    k=centroid[j]
    if k ge ny then k=ny-1 
    while( im[j,k] ge threshold*im[j,centroid[j]]) do k++
    upper[j]=k
    k=centroid[j]
    while( im[j,k] ge threshold*im[j,centroid[j]]) do k--
    lower[j]=k
  endfor
  upper[0]=upper[1] ;valid for simulated data
  lower[0]=lower[1] ;valid for simulated data
  upper[nx-1]=upper[nx-2] ;valid for simulated data
  lower[nx-1]=lower[nx-2] ;valid for simulated data
  
  k=0
  m_upper=dblarr(nx-8)
  m_lower=dblarr(nx-8)
  for j=0, nx-9 do begin
    m_upper[k]=max(upper[j:j+7])
    m_lower[k]=min(lower[j:j+7])
    k++
  endfor
  dev=2
  uu_thresh=median(m_upper)+dev*stddev(m_upper)
  ul_thresh=median(m_upper)-dev*stddev(m_upper)
  lu_thresh=median(m_lower)+dev*stddev(m_lower)
  ll_thresh=median(m_lower)-dev*stddev(m_lower)
  pix=indgen(n_elements(m_upper))

  ;now to find outliers in the max values
  for j=0, n_elements(m_upper)-1 do begin
    if (ul_thresh ge m_upper[j] le uu_thresh) then begin
      if j ge 4 and j le n_elements(m_upper)-5 then  m_upper[j]= median(m_upper[j-4:j+4]) $
        else if j lt 4 then m_upper[j]= median(m_upper[j:j+8]) $
        else m_upper[j]= median(m_upper[j-8:j])
    endif
    if (ll_thresh ge m_lower[j] le lu_thresh) then begin
      if j ge 4 and j le n_elements(m_lower)-5 then  m_lower[j]= median(m_lower[j-4:j+4]) $
        else if  j lt 4 then m_lower[j]= median(m_lower[j:j+8]) $
        else m_lower[j]= median(m_lower[j-8:j])
    endif
  endfor
  pixel=indgen(nx)
  degree=1
  mn_upper=dblarr(n_elements(m_upper))
  mn_lower=dblarr(n_elements(m_lower))
  fit=poly_fit(pix,m_upper,2,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma_fit,$
               STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_upper+=fit[i]*pix^i

  fit=poly_fit(pix,m_lower,2,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma_fit,$
               STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_lower+=fit[i]*pix^i


  n_upper=fix(interpol(mn_upper,pix,pixel))+1
  n_lower=fix(interpol(mn_lower,pix,pixel))
  
  out_frame_up=where(n_upper gt ny-1)
  out_frame_down=where(n_lower lt 0)
  if total(out_frame_up) ge 0 then n_upper[out_frame_up] = ny-1 
  if total(out_frame_down) ge 0 then n_lower[out_frame_down] = 0
  for j=0, nx-1 do begin
     spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]], /NAN)
     noise[j]=total(sigma[j,n_lower[j]:n_upper[j]]^2, /NAN)
     dq_val[j]=dq_val[j] or arrayor(dq[j,n_lower[j]:n_upper[j]], 2)
     ;dq_val[j]=total(dq[j,n_lower[j]:n_upper[j]], /NAN
     ;for background use the input from lookup table as well
     if multiple eq 0 then begin
       if ushift ne 0 then begin 
         if ((n_upper[j]+ushift) le ny-1) then begin
           back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
           noise_u[j]=total(sigma[j,n_lower[j]+ushift:n_upper[j]+ushift]^2, /NAN)
           ;dq_bval[j]+=total(dq[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
           dq_bval[j]=dq_bval[j] or arrayor(dq[j,n_lower[j]+ushift:n_upper[j]+ushift],2)
         endif else begin
           errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                   +' Please re-enter a new upper offset for backgrounds' ,logonly=logonly
           message,' Background offset specified is outside the CCD area.'$
                  +' Please re-enter a new upper offset for backgrounds'
         endelse
         background[j]=(back_u[j])
         backg_error[j]=(sqrt(noise_u[j]))
       endif
       if lshift ne 0 then begin 
         if ((n_lower[j]-lshift) ge 0) then begin
           back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
           noise_l[j]=total(sigma[j,n_lower[j]-lshift:n_upper[j]-lshift]^2, /NAN)
           dq_bval[j]=dq_bval[j] or arrayor(dq[j,n_lower[j]-lshift:n_upper[j]-lshift],2)
           ;dq_bval[j]+=total(dq[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
         endif else begin
           errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                   +' Please re-enter a new lower offset for backgrounds',logonly=logonly
           message,'Background offset specified is outside the CCD area.'$
                  +' Please re-enter a new lower offset for backgrounds'
         endelse
         background[j]=(back_l[j])
         backg_error[j]=(sqrt(noise_l[j]))
       endif
     endif else begin
       for k=0,n_elements(ushift)-1 do begin
         if single_upix eq 0 then begin
           ubgstr=centroid[j]+ushift[k]-uwidth[k]
           ubgstop=centroid[j]+ushift[k]+uwidth[k]
         endif else begin
           ubgstr=centroid[j]+ushift[k]
           ubgstop=centroid[j]+ushift[k]+1
         endelse
         if (ubgstop le ny-2) then begin
           back_u[j]=total(im[j,ubgstr:ubgstop], /NAN)
           noise_u[j]=total(sigma[j,ubgstr:ubgstop], /NAN)
           dq_bval[j]=dq_bval[j] or arrayor(dq[j,ubgstr:ubgstop],2)
           ;dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
         endif else begin
           errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                   +' Please re-enter a new upper offset for backgrounds',logonly=logonly
           message,'Background offset specified is outside the CCD area.'$
                  +' Please re-enter a new upper offset for backgrounds'
         endelse
       endfor
       for k=0,n_elements(lshift)-1 do begin
         if single_lpix eq 0 then begin
           lbgstr=centroid[j]+lshift[k]-lwidth[k]
           lbgstop=centroid[j]+lshift[k]+lwidth[k]
         endif else begin
           lbgstr=centroid[j]+lshift[k]
           lbgstop=centroid[j]+lshift[k]+1
         endelse
         if (lbgstr ge 1) then begin
           back_u[j]=total(im[j,lbgstr:lbgstop], /NAN)
           noise_u[j]=total(sigma[j,lbgstr:lbgstop], /NAN)
           dq_bval[j]=dq_bval[j] or arrayor(dq[j,lbgstr:lbgstop],2)
           ;dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
         endif else begin
           errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                   +' Please re-enter a new lower offset for backgrounds',logonly=logonly
           message,'Backround offset specified is outside the CCD area.'$
                  +' Please re-enter a new lower offset for backgrounds'
         endelse
       endfor
       bg_length=double(ubg_length+lbg_length)
       spec_length=double(n_upper[j]-n_lower[j])
       bgscale=float(bg_length/spec_length)
       background[j]=(back_u[j]+back_l[j])/bgscale
       backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
     endelse
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=dq_val(prb) or 128b
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=dq_bval(prbb) or 128b
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,$
            lower:n_lower,upper:n_upper,dq:dq_val,dq_bcg:dq_bval}

  return,spectrum
end

;extraction with pixel dependend centroid and a cross dispersion function
function extract_func,im,sigmaim,dq,centroid,threshold,back_trace
  err = ''
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    logprint,'Syntax - extract_func,im,sigma,centroid,threshold,',logonly = logonly
    message,'Syntax - extract_func,im,sigma,centroid,threshold'
    err = '%extract_func: Insufficient number of parameters'
  endif
  nxy=size(im)
  if nxy[0] ne 2 then begin
    print,'Input file error, the image file is not 2-D'
    err = '%extract_func: Input file error'
  endif
  nxy2=size(sigmaim)
  if nxy2[0] ne 2 then begin
    print,'Input file error, the error file is not 2-D'
    err = '%extract_func: Input file error'
  endif
  nx=nxy[1]
  ny=nxy[2]
  multiple=back_trace.type
  if (datatype(back_trace.shift_upper) eq 'INT') then ushift=fix(back_trace.shift_upper)
  if (datatype(back_trace.shift_lower) eq 'INT') then lshift=fix(back_trace.shift_lower)
  if (datatype(back_trace.shift_up_width) eq 'INT') then begin
    uwidth=fix(back_trace.shift_up_width)
    if (total(uwidth) eq n_elements(uwidth)) then begin
      single_upix=1
      ubg_length=double(n_elements(uwidth))
    endif else begin
      single_upix=0
      ubg_length=double(2*total(uwidth)+n_elements(uwidth))
    endelse
  endif
  if (datatype(back_trace.shift_low_width) eq 'INT') then begin
    lwidth=fix(back_trace.shift_low_width)
    if (total(lwidth) eq n_elements(width)) then begin
      single_lpix=1
      lbg_length=double(n_elements(width))
    endif else begin
      single_lpix=0
      lbg_length=double(2*total(width)+n_elements(lwidth))
    endelse
  endif  
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  dq_val=bytarr(nx)
  dq_bval=bytarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  upper=intarr(nx)
  lower=intarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  if n_elements(centroid) ne nx then begin
    logprint,'Input dimension error, centroid has different dimension than image array'$
            ,logonly = logonly
    message,'Input dimension error, centroid has different dimension than image array'
    err = '%extract_func: Input diamension error'

  endif
  cent_centroid=fix(mean(centroid))
  new_im=im[*,cent_centroid-20:cent_centroid+20]
  x=indgen(41,start=-20)
  for j=1,nx-2 do begin ;valid for simulated data change to: j=0,nx-1 for actual data
    y=new_im[j,*]
    params=mpfitpeak(x,y,a)
    peak_loc=where(im[j,*] eq a[0])
    peak_loc=a[1]
    peak_loc_pixel=(where(x eq fix(peak_loc)))
    k=peak_loc_pixel
    while( y[k] ge threshold*a[0] and k le peak_loc_pixel+20) do k++
    upper[j]=cent_centroid+x[peak_loc_pixel]+(k-peak_loc_pixel)
    k=peak_loc_pixel
    while( y[k] ge threshold*a[0] and k ge peak_loc_pixel-20) do k--
    lower[j]=cent_centroid-x[peak_loc_pixel]-(peak_loc_pixel-k)
  endfor
  centroid[j]=peak_loc
  k=0
  upper[0]=upper[1] ;valid for simulated data
  lower[0]=lower[1] ;valid for simulated data
  upper[nx-1]=upper[nx-2] ;valid for simulated data
  lower[nx-1]=lower[nx-2] ;valid for simulated data
  m_upper=dblarr(nx-8)
  m_lower=dblarr(nx-8)

  for j=0, nx-9 do begin
    m_upper[k]=max(upper[j:j+7])
    m_lower[k]=min(lower[j:j+7])
    k++
  endfor
  dev=2
  uu_thresh=median(m_upper)+dev*stddev(m_upper)
  ul_thresh=median(m_upper)-dev*stddev(m_upper)
  lu_thresh=median(m_lower)+dev*stddev(m_lower)
  ll_thresh=median(m_lower)-dev*stddev(m_lower)
  pix=indgen(n_elements(m_upper))

  ;now to find outliers in the max values
  for j=0, n_elements(m_upper)-1 do begin
    if (ul_thresh ge m_upper[j] le uu_thresh) then begin
      if j ge 4 and j le n_elements(m_upper)-5 then  m_upper[j]= mean(m_upper[j-4:j+4]) $
      else if j lt 4 then m_upper[j]= mean(m_upper[j:j+8]) else m_upper[j]= mean(m_upper[j-8:j])
    endif
    if (ll_thresh ge m_lower[j] le lu_thresh) then begin
      if j ge 4 and j le n_elements(m_lower)-5 then  m_lower[j]= mean(m_lower[j-4:j+4]) $
        else if  j lt 4 then m_lower[j]= mean(m_lower[j:j+8]) else m_lower[j]= mean(m_lower[j-8:j])
    endif
  endfor

  pixel=indgen(nx)
  degree=1
  mn_upper=dblarr(n_elements(m_upper))
  mn_lower=dblarr(n_elements(m_lower))
  fit=poly_fit(pix,m_upper,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma,$
               STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_upper+=fit[i]*pix^i

  fit=poly_fit(pix,m_lower,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma,$
               STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_lower+=fit[i]*pix^i


  n_upper=fix(interpol(mn_upper,pix,pixel))+1
  n_lower=fix(interpol(mn_lower,pix,pixel))
  n_upper=fix(interpol(m_upper,pix,pixel))
  n_lower=fix(interpol(m_lower,pix,pixel))
  diffu=intarr(n_elements(n_upper))
  diffl=intarr(n_elements(n_upper))
  out_frame_up=where(n_upper gt ny-1)
  out_frame_down=where(n_lower lt 0)
  if total(out_frame_up) ge 0 then n_upper[out_frame_up] = ny-1
  if total(out_frame_down) ge 0 then n_lower[out_frame_down] = 0
  for j=0, nx-1 do begin
    spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]], /NAN)
    noise[j]=total(sigmaim[j,n_lower[j]:n_upper[j]]^2, /NAN)
    dq_val[j]=dq_val[j] or arrayor(dq[j,n_lower[j]:n_upper[j]],2)
    ;dq_val[j]=total(dq[j,n_lower[j]:n_upper[j]])
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ushift ne o0 then begin
        if ((n_upper[j]+ushift) le ny) then begin
          back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
          noise_u[j]=total(sigmaim[j,n_lower[j]+ushift:n_upper[j]+ushift]^2, /NAN)
          dq_bval[j]=dq_bval[j] or arrayor(dq[j,n_lower[j]+ushift:n_upper[j]+ushift],2)
          ;dq_bval[j]+=total(dq[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                  +' Please re-enter a new upper offset for backgrounds' ,logonly=logonly
          message,' Background offset specified is outside the CCD area.'$
                 +' Please re-enter a new upper offset for backgrounds'
        endelse
        background[j]=(back_u[j])
        backg_error[j]=(sqrt(noise_u[j]))
      endif 
      if lshift ne 0 then begin 
        if ((n_lower[j]-lshift) ge 0) then begin
          back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
          noise_l[j]=total(sigmaim[j,n_lower[j]-lshift:n_upper[j]-lshift]^2, /NAN)
          dq_bval[j]=dq_bval[j] or arrayor(dq[j,n_lower[j]-lshift:n_upper[j]-lshift],2)
          ;dq_bval[j]+=total(dq[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                  +' Please re-enter a new lower offset for backgrounds',logonly=logonly
          message,'Backround offset specified is outside the CCD area.'$
                 +' Please re-enter a new lower offset for backgrounds'
        endelse
        background[j]=(back_l[j])
        backg_error[j]=(sqrt(noise_l[j]))
      endif
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if single_upix eq 0 then begin
          ubgstr=centroid[j]+ushift[k]-uwidth[k]
          ubgstop=centroid[j]+ushift[k]+uwidth[k]
        endif else begin
          ubgstr=centroid[j]+ushift[k]
          ubgstop=centroid[j]+ushift[k]+1
        endelse
        if (ubgstop le ny-2) then begin
          back_u[j]=total(im[j,ubgstr:ubgstop], /NAN)
          noise_u[j]=total(sigmaim[j,ubgstr:ubgstop], /NAN)
          dq_bval[j]=dq_bval[j] or arrayor(dq[j,ubgstr:ubgstop],2)
          ;dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                  +' Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,'Background offset specified is outside the CCD area.'$
                 +' Please re-enter a new upper offset for backgrounds'
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if single_lpix eq 0 then begin
          lbgstr=centroid[j]+lshift[k]-lwidth[k]
          lbgstop=centroid[j]+lshift[k]+lwidth[k]
        endif else begin
          lbgstr=centroid[j]+lshift[k]
          lbgstop=centroid[j]+lshift[k]+1
        endelse
        if (lbgstr gt 0) then begin
          back_u[j]=total(im[j,lbgstr:lbgstop], /NAN)
          noise_u[j]=total(sigmaim[j,lbgstr:lbgstop], /NAN)
          dq_bval[j]= dq_bval[j] or arrayor(dq[j,lbgstr:lbgstop],2)
          ;dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area.'$
                  +' Please re-enter a new lower offset for backgrounds',logonly=logonly
          message,'Background offset specified is outside the CCD area.'$
                 +' Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=double(ubg_length+lbg_length)
      spec_length=double(n_upper[j]-n_lower[j])
      bgscale=float(bg_length/spec_length)
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=dq_val(prb) or 128b
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=dq_bval(prbb) or 128b
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,$
            lower:n_lower,upper:n_upper,dq:dq_val,dq_bcg:dq_bval}
  return,spectrum
end

function control_trace,in_image,infile,trace_type,filename
  idl_ver=float(!Version.RELEASE)
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    print,'Syntax - control_trace,in_image,infile,trace_type,,filename
    err = '%CONTROL_TRACE: Insufficient number of parameters'

  endif
  if keyword_defined(trace_type) eq 0 then begin
    logprint,'CONTROL_TRACE: Type input for trace not defined using default: simple'
    trace_type='simple'
  endif

  if(tag_exist(infile,'data_path') eq 1) then data_path=infile.data_path
  if(tag_exist(infile,'temp_path') eq 1) then inter_path=infile.temp_path
  if n_elements(data_path) eq 0 then begin
    dpflg=1
    CD, Current=data_path
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': data_path=data_path+'\' ;WINDOWS
      'UNIX': data_path=data_path+'/'; UNIX.
    ENDCASE
  endif else data_path=detectos(data_path)
  if n_elements(data_path) eq 0 then begin
    CD, Current=data_path
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': data_path=data_path+'\' ;WINDOWS
      'UNIX': data_path=data_path+'/'; UNIX.
    ENDCASE
  endif else data_path=detectos(data_path)

  if n_elements(inter_path) eq 0 then begin
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': inter_path=data_path+'temp\' ;WINDOWS
      'UNIX': inter_path=data_path+'temp/'; UNIX.
    ENDCASE
    logprint,'CONTROL: No temporary file path found.'$
            +' Temporary file directory created in data directory'
  endif else inter_path=detectos(inter_path)

  im=in_image.data
  error=in_image.error
  hdr=in_image.header
  dq=in_image.dq
  ycut1=SXPAR( hdr, 'YCUT1')
  ycut2=SXPAR( hdr, 'YCUT2')

  constu1=1
  constu2=0
  constu3=1
  constl1=1
  constl2=0
  constl3=1

  if(tag_exist(infile,'cent_poly_degree') eq 1) then degree=infile.cent_poly_degree else begin
    degree=1  ;setting degree of polynomial for centroid default is 1
    logprint,'CONTROL_TRACE: Polynomial degree for centroid assumed to be one.'
  endelse
  nxy=size(im)
  nx=nxy[1]
  ny=nxy[2]
  k=0

  if(tag_exist(infile,'centroid') eq 1) then begin
    cent_value=double(infile.centroid)
    cent_value=cent_value-ycut1

    if (cent_value le 0  or double(infile.centroid) ge ycut2) then begin
      logprint,'CONTROL_TRACE: Centroid value specified in the parameter file is outside the frame.'
      return,0
    endif
    centroid=make_array(nx,value=cent_value)
    pixel=indgen(nx)
    y=indgen(ny)
    border=15
    low_range=centroid[0]-border
    high_range=centroid[0]+border
    
    if (cent_value le border or cent_value ge ny-border+1) then begin
      logprint,'CONTROL_TRACE: Centroid value specified in the parameter file is'$
              +' close to the edge of the frame.'
      logprint,'Press q to skip the trace and extraction for current spectrum.'$
              +' Press any key to continue trace and extraction for the current spectrum.'
      R = GET_KBRD()
      if R eq 'q' then begin
        logprint,'CONTROL: Skipping the trace and extraction for current spectrum as'$
                +' requested by the user.'
        return,0
      endif else begin
        if cent_value le border then low_range = 0 else if cent_value ge ny-border+1 then $
           high_range = ny-1
      endelse
    endif
    if tag_exist(infile,'slope') eq 1 then slope=float(infile.slope)
    if ny le 100 then implot=im else implot=im[*,centroid[0]-50:centroid[0]+50]
    if ny le 100 then yranges = [0,ny] else yranges= [centroid[0]-50,centroid[0]+50]    
    centroidp=centroid+(slope*pixel)
    window,xsize=1800,ysize=600
    cgimage,implot,/AXES,/Save,yrange=yranges,xrange=[pixel[0]$
           ,pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2.5,stretch=clip
    cgoplot,pixel,centroidp,color='red';,yrange=[245,265]
    cgLegend, Title=['Centroid'], LineStyle=[0], SymSize=2, Color=['red']$
            , Location=[0.8,0.88],thick=1.5,Length=0.05$
            , VSpace=2.0, /Box,/Background, BG_Color='white'
    write_png,inter_path+filename+'_centroid.png',TVRD(/TRUE)
  endif else begin
    logprint,'CONTROL_TRACE: centroid not specified. Calculating on the go'
    pixel=indgen(nx)
    rows=indgen(ny)
    mim=median(im,dimension=1)
    fit=linfit(rows,mim)
    fits=fit[0]+fit[1]*rows
    diff=mim-fits
    mdiff=max(diff[0:-2],locc)
    if ny le 100 then im_trim=im else im_trim=im[*,locc-50:locc+50]
    t_cent=dblarr(nx-8+1)
    for i=0, nx-8 do begin
      c_sum=median(im_trim[i:i+7,*],dimension=1)
      m=max(c_sum,loc)
      t_cent[k]=loc
      k++
    endfor
    trace_cen=dblarr(n_elements(t_cent))
    pix=indgen(n_elements(t_cent))
    fit=poly_fit(pix,t_cent,degree,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror $
                ,SIGMA=sigma,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
    for i=0,degree do trace_cen+=fit[i]*pix^i
    centroid_full=interpol(trace_cen,pix,pixel)
    centroid=fix(centroid_full)
    if ny le 100 then centroid = centroid else centroid=centroid+locc-50
    y=indgen(ny)
    if ny le 100 then implot=im else implot=im[*,centroid[0]-50:centroid[0]+50]
    if ny le 100 then yranges = [0,ny] else yranges= [centroid[0]-50,centroid[0]+50]
    window,xsize=1800,ysize=600
    cgimage,implot,/AXES,/Save,yrange=yranges,xrange=[pixel[0],pixel[nx-1]],$
            xtitle='X pixels',ytitle='Y pixels',charsize=2.5,stretch=clip
    cgoplot,pixel,centroid+0.5,color='red';,yrange=[245,265]
    cgLegend, Title=['Centroid'], LineStyle=[0], SymSize=2, Color=['red']$
            , Location=[0.8,0.88] ,thick=1.5,Length=0.05$
            , VSpace=2.0, /Box,/Background, BG_Color='white'
    write_png,inter_path+filename+'centroid.png',TVRD(/TRUE)
  endelse

  ;background stuff
  if tag_exist(infile,'background_trace') eq 0 then begin
    errorlog,'CONTROL_TRACE: Background trace information is required.'$
            +' Please re-run the simulator with the same',logonly=1
    message,'CONTROL_TRACE: Background trace information is required.'$
           +' Please re-run the simulator with the same'
  endif
  if (file_test(infile.background_trace)eq 0) then begin
    logprint,'CONTROL_TRACE: Background trace input is pixel deviations from centroid'
    str_bcg=strsplit(infile.background_trace,',',/EXTRACT)
    if n_elements(str_bcg) eq 1 then begin
      shift_val=fix(str_bcg)
      if shift_val ge 0 then begin 
        ushift = shift_val
        lshift = 0
      endif else begin
        ushift = 0
        lshift = -1*shift_val
      endelse    
      if ushift eq 0 and lshift eq 0 then begin
        errorlog,'CONTROL_TRACE:Shift value for background trace is 0 or not specified'
        message,'CONTROL_TRACE:Shift value for background trace is 0 or not specified'
      endif
      multiple=0
    endif else begin
      if (idl_ver ge 8) then begin
        bcg_lt=strsplit(str_bcg,':',/EXTRACT)
        str_bcg_arr=bcg_lt.toarray()
        shift_val=fix(str_bcg_arr[*,0])
        bcg_arr_dim=size(str_bcg_arr)
        if bcg_arr_dim[0] ge 2 then shift_width=fix(str_bcg_arr[*,1]) $
        else shift_width = make_array(bcg_arr_dim[1],Value=1,/INTEGER) 
      endif else begin
        shift_val=intarr(n_elements(str_bcg))
        shift_width=intarr(n_elements(str_bcg))
        for jj=0,n_elements(str_bcg)-1 do begin
          bcg_lt=strsplit(str_bcg[jj],':',/EXTRACT)
          shift_val[jj]=bcg_lt[0]
          if n_elements(bcg_lt) gt 1 then shift_width[jj]=bcg_lt[1] else shift_width[jj]=1
        endfor
      endelse    
      ush_loc=where(shift_val ge 0,/NULL)
      lsh_loc=where(shift_val lt 0,/NULL)
      if ush_loc eq !NULL then undefine,ushift,uwidth else begin 
        ushift=shift_val[ush_loc]
        uwidth=shift_width[ush_loc]
      endelse  
      if lsh_loc eq !NULL then undefine,lshift,lwidth else begin 
        lshift=shift_val[lsh_loc]
        lwidth=shift_width[lsh_loc]
      endelse  
      multiple=n_elements(str_bcg)
    endelse
  endif else begin
    bt_file=detectos(infile.background_trace)
    readcol,bt_file,star,ushift,lshift,u_offset,l_offset,F='A,I,I,D,D'
    ;fnd the target
    target=SXPAR( hdr, 'TARGNAME')
    s=STRCOMPRESS(target, /REMOVE_ALL)
    bg_loc=WHERE(STRMATCH(star, s, /FOLD_CASE) EQ 1)

    if (n_elements(bg_loc) eq 1) then multiple=0 else multiple=bg_loc
    if ycut1 le 171 then begin
      ushift=fix(constu1*u_offset+ushift)
      lshift=fix(constl1*l_offset+lshift)
    endif else if ycut1 le 343 then begin
      ushift=fix(constu2*u_offset+ushift)
      lshift=fix(constl2*l_offset+lshift)
    endif else begin
      ushift=fix(constu3*u_offset+ushift)
      lshift=fix(constl3*l_offset+lshift)
    endelse
  endelse
  if (idl_ver ge 8) then begin
    if typename(lshift) eq 'UNDEFINED' then lshift = 'undefined'
    if typename(ushift) eq 'UNDEFINED' then ushift = 'undefined'
    if typename(uwidth) eq 'UNDEFINED' then uwidth = 'undefined'
    if typename(lwidth) eq 'UNDEFINED' then lwidth = 'undefined'
  endif else begin
    if datatype(lshift,2) eq 0 then lshift = 'undefined'
    if datatype(ushift,2) eq 0 then ushift = 'undefined'
  endelse
  back_trace={type:multiple,shift_upper:ushift,shift_lower:lshift,shift_up_width:uwidth, $
              shift_low_width:lwidth}
  ;default extraction of box (similar to COS extraction)
  ;impliment slope change with extraction cut
  if tag_exist(infile,'slope') eq 1 then slope=float(infile.slope)
  if tag_exist(infile,'width') eq 1 then width=fix(infile.width)
  if tag_exist(infile,'upper') eq 1 then upper=fix(infile.upper)
  if tag_exist(infile,'lower') eq 1 then lower=fix(infile.lower)
  if tag_exist(infile,'threshold') eq 1 then threshold=fix(infile.threshold)


  logprint,'CONTROL_TRACE: Centroid value used is '+strtrim(string(centroid[0]),2)+$
           ' at the edge of the spectrum.'
  if datatype(slope,2) ne 0 then logprint,'Slope value used is '+strtrim(string(slope),2)+'.'
  if datatype(width,2) ne 0 then if width gt 0 then $
    logprint,'CONTROL_TRACE: Width value used for extraction type simple is '$
             +strtrim(string(width),2)+'.'
  if datatype(upper,2) ne 0 then $
    logprint,'CONTROL_TRACE: Upper width value used for extraction type fixed is '$
            +strtrim(string(fix(upper)),2)+'.'
  if datatype(lower,2) ne 0 then $
    logprint,'CONTROL_TRACE: Lower width value used for extraction type fixed is '$
            +strtrim(string(fix(lower)),2)+'.'
  if datatype(threshold,2) ne 0 then $
    logprint,'CONTROL_TRACE: Threshold value used for extraction type variable/function is '$
            +strtrim(string(double(threshold)),2)+'.'

  if datatype(slope,2) eq 0 then begin
    logprint,'Slope value not found. Using default value of -9.76e-4.'
    slope=-9.76e-4
  endif
  if datatype(width,2) eq 0 then begin
    logprint,'Width value not found. Using default value of 10.'
    width=10
  endif
  if(tag_exist(infile,'centroid') eq 0) then begin
    logprint,'Slope value not required as centroid calculated on the go'$
            +' and have slope information. Setting it to 0.'
    slope=0
  endif  
  ;def_spectrum = extract_box(im,error,dq,centroid,slope,width,back_trace)
  centroid_def=centroid
  ;extraction
  case trace_type of
    'simple'  : begin
      extraction = extract_box(im,error,dq,centroid,slope,width,back_trace)
      sxaddpar, hdr,'EXTRTYP', 'simple',   'Type of extraction method'
      sxaddpar, hdr,'EXTRCNT', extraction.centroid[0],'Centroid of extraction spectrum'
      sxaddpar, hdr,'EXTPAR1', slope,      'Extraction parameter 1'
      sxaddpar, hdr,'EXTPAR2', width,      'Extraction parameter 2'
      sxaddpar, hdr,'EXTPAR3', 'NA',       'Extraction parameter 3'
    end
    'fixed'   : begin
      if tag_exist(infile,'upper') eq 1 then upper=fix(infile.upper) else begin
        logprint,'Upper width value not found. Using default value of 10.'
        upper=10
      endelse
      if tag_exist(infile,'lower') eq 1 then lower=fix(infile.lower) else begin
        logprint,'Lower width value not found. Using default value of 10.'
        lower=10
      endelse
      extraction = extract_sum(im,error,dq,centroid,slope,upper,lower,back_trace)
      sxaddpar, hdr,'EXTRTYP', 'fixed',    'Type of extraction method'
      sxaddpar, hdr,'EXTRCNT', extraction.centroid[0],'Centroid of extraction spectrum'
      sxaddpar, hdr,'EXTPAR1', slope,      'Extraction parameter 1'
      sxaddpar, hdr,'EXTPAR2', upper,      'Extraction parameter 2'
      sxaddpar, hdr,'EXTPAR3', lower,      'Extraction parameter 3'
    end
    'variable': begin
      if tag_exist(infile,'threshold') eq 1 then threshold=double(infile.threshold) else begin
        logprint,'Threshold value not found. Using default value of 0.01'
        threshold=0.01
      endelse
      extraction = extract_varsum(im,error,dq,centroid,threshold,back_trace)
      sxaddpar, hdr,'EXTRTYP', 'variable', 'Type of extraction method'
      sxaddpar, hdr,'EXTRCNT', extraction.centroid[0],'Centroid of extraction spectrum'
      sxaddpar, hdr,'EXTPAR1', threshold,  'Extraction parameter 1'
      sxaddpar, hdr,'EXTPAR2', mean(extraction.lower),'Extraction parameter 2'
      sxaddpar, hdr,'EXTPAR3', mean(extraction.upper),'Extraction parameter 3'
    end
    'function': begin
      if tag_exist(infile,'threshold') eq 1 then threshold=double(infile.threshold) else begin
        logprint,'Threshold value not found. Using default value of 0.01'
        threshold=0.01
      endelse
      extraction = extract_func(im,error,dq,centroid,threshold,back_trace)
      sxaddpar, hdr,'EXTRTYP', 'function', 'Type of extraction method'
      sxaddpar, hdr,'EXTRCNT', extraction.centroid[0],'Centroid of extraction spectrum'
      sxaddpar, hdr,'EXTPAR1', threshold,  'Extraction parameter 1'
      sxaddpar, hdr,'EXTPAR2', mean(extraction.lower),'Extraction parameter 2'
      sxaddpar, hdr,'EXTPAR3', mean(extraction.upper),'Extraction parameter 3'
    end
    else :      begin
      logprint,'CONTROL_TRACE: Invalid type input for trace: Please recheck your input'$
              ,logonly=logonly
      message,' Invalid type input for trace: Please recheck your input'
    end
  endcase
  if ny le 100 then implot=im else implot=im[*,centroid[0]-50:centroid[0]+50]
  if ny le 100 then yranges = [0,ny] else yranges= [centroid[0]-50,centroid[0]+50]
  
  window,xsize=2000,ysize=500
  cgimage,implot,/AXES,/Save,yrange=yranges,xrange=[pixel[0]$
         ,pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2,CTINDEX=0,$
         stretch=2 
  cgoplot,pixel,extraction.lower,color='red',symsize=2;,yrange=[245,265]
  cgoplot,pixel,extraction.upper+1,color='red',symsize=2;,yrange=[245,265]
  if back_trace.type eq 0 then begin
    if back_trace.shift_upper eq 0 then begin
      cgoplot,pixel,extraction.lower-back_trace.shift_lower[0],color='green',symsize=2;,yrange=[245,265]
      cgoplot,pixel,extraction.upper+1-back_trace.shift_lower[0],color='green',symsize=2;,yrange=[245,265]
    endif
  endif
  cgLegend, Title=['Extraction region', 'Background region'], LineStyle=[0,0], SymSize=2, Color=['red','green']$
          , Location=[0.77,0.89],thick=1.5,Length=0.05$
          , VSpace=2.0, /Box,/Background, BG_Color='white'

  write_png,inter_path+filename+'_trace.png',TVRD(/TRUE)
  
  colaps=median(im,dimension=1)
  colaps=median(im[0:100,*],dimension=1)
  ypix=indgen(n_elements(colaps))
  window,xsize=1400,ysize=700
  cgplot, ypix,(colaps),xrange=[min(ypix),max(ypix)],xtitle='Y pixels',ytitle='counts'$ ;[10!u4!n]
        , Title='Target and background  extraction region.',symsize=2,charsize=2.5$
        , charthick=2.5,xthick=2.5,ythick=2.5
  extr_low=min(extraction.lower[0:100])
  extr_hi=max(extraction.upper)
  ;extr_low = extraction.lower[100]
  extr_hi = extraction.upper[100]
  cgoplot,[extr_low,extr_low],[!Y.CRange[0],!Y.CRange[1]], color='red'
  cgoplot,[extr_hi,extr_hi],[!Y.CRange[0],!Y.CRange[1]], color='red'
  if multiple eq 0 then begin
    if lshift ne 0 then begin
      cgoplot,[extr_low-lshift,extr_low-lshift],[!Y.CRange[0],!Y.CRange[1]], color='green'
      cgoplot,[extr_hi-lshift,extr_hi-lshift],[!Y.CRange[0],!Y.CRange[1]], color=' green'
    endif
    if ushift ne 0 then begin
      cgoplot,[extr_low+ushift,extr_low+ushift],[!Y.CRange[0],!Y.CRange[1]], color='green'
      cgoplot,[extr_hi+ushift,extr_hi+ushift],[!Y.CRange[0],!Y.CRange[1]], color='green'
    endif  
  endif else begin
    for i=0, n_elements(ushift)-1 do begin
      bg_plot_loc=ushift[i]+extraction.centroid[-1]
      cgoplot,[bg_plot_loc,bg_plot_loc],[!Y.CRange[0],!Y.CRange[1]], color='green'
    endfor
    for i=0, n_elements(lshift)-1 do begin
      bg_plot_loc=extraction.centroid[0]+lshift[i]
      cgoplot,[bg_plot_loc,bg_plot_loc],[!Y.CRange[0],!Y.CRange[1]], color='green'
    endfor
  endelse
  
  write_png,inter_path+filename+'_trace2.png',TVRD(/TRUE)
  ext_data=extraction.data
  ext_error=extraction.error
  ext_bck=extraction.bck
  ext_bck_error=extraction.bck_error
  ext_dq=extraction.dq
  ext_bg_dq=extraction.dq_bcg
  ;comparison to make sure extraction is proper, comparison with default
  ;plt=plot(ext_data)
  ;comparison of pixel vs extraction region
  sd=dblarr(52)
  fl=dblarr(52)
  pi=dblarr(52)
  j=0
;  for i=1, 20 do begin
;    cen_ext=centroid_def
;    extractions = extract_box(im,error,dq,cen_ext,slope,i,back_trace)
;    flux=extractions.data
;    unc=extractions.error
;    if (idl_ver ge 8) then begin
;      continum_flux=flux[-150:-10]
;      continum_error=unc[-150:-10]
;    endif else begin
;      continum_flux=flux[n_elements(flux)-150:n_elements(flux)-10]
;      continum_error=unc[n_elements(flux)-150:n_elements(flux)-10]
;    endelse
;      
;    x=indgen(n_elements(continum_flux))
;    tp=trapz_error(x,continum_flux,continum_error)
;    sd[j]=stddev(continum_flux)
;    fl[j]=tp[0]
;    pi[j]=i
;    j++
;  endfor
;  window,xsize=1400,ysize=700
;  cgplot, pi,(fl/1E4),psym=15,Color='red7',xtitle='pixels from centroid' ,symsize=2$
;        , title='Variation of flux & uncertanity vs extraction window (1898-2038 pix)', $
;        charsize=2.5,charthick=2.5,xthick=2.5,ythick=2.5, $
;         ytitle='Flux of extracted spectrum [10!u4!n]'$
;          ,YRange=[double(min(fl/1E4))-0.01, double(max(fl/1E4))+1],yticks=6
;        
;  cgAxis, YAxis=1,ytitle=' Sigma of extracted spectrum', YRange=[double(min(sd))-10$
;        , double(max(sd))+10],charsize=2.5,charthick=2.5,xthick=2.5,ythick=2.5,yticks=6,/Save
;  cgoplot,pi,sd,psym=16,Color='blue',symsize=2
;  write_png,inter_path+filename+'_pixel_vs_flux.png',TVRD(/TRUE)
  sxaddpar, hdr,'EXTFLF',1,'Spectrum extraction flag'
  
  spectrum={data:ext_data,error:ext_error,header:hdr,background:ext_bck,bck_error:ext_bck_error,$
            dq:ext_dq,dq_bg:ext_bg_dq}
            
  return,spectrum
end
