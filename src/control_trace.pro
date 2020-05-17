;data quality carried out only for spectrum. Modification of background and data quality remaining


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
    logprint,'CONTROL_TRACE: centroid for box extraction not specified. Calculating on the go',logonly=logonly
    colsum=total(im,1)
    maximum=max(colsum[0],loc)
    centroid=make_array(nx,value=loc)
  endelse

 if width le 0 then begin
  logprint,'Width specified is less than or equal to zero. CONTROL will reset it to default value of 10.'
  width =10
 endif
  stp= fix(centroid+slope*pixel+width)
  str=fix(centroid+slope*pixel-width)
  out_frame_up=where(stp gt ny-1)
  out_frame_down=where(str lt 0)
  for j=0,nx-1 do begin
    spectrum_val[j]=total(im[j,str[j]:stp[j]], /NAN)
    noise[j]=total(sigma[j,str[j]:stp[j]]^2, /NAN)
    dq_val[j]=total(dq[j,str[j]:stp[j]], /NAN)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((stp[j]+ushift) le ny-1) then begin
        back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift], /NAN)
        noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2, /NAN)
        dq_bval[j]+=total(dq[j,str[j]+ushift:stp[j]+ushift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
      endelse
      if ((str[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift], /NAN)
        noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2, /NAN)
        dq_bval[j]+=total(dq[j,str[j]-lshift:stp[j]-lshift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if ((centroid[j]+ushift[k]) le ny-2) then begin
          back_u[j]=total(im[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          noise_u[j]=total(sigma[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if ((centroid[j]-lshift[k]) ge 1) then begin
          back_l[j]=total(im[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          noise_l[j]=total(sigma[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=double(n_elements(ushift)+n_elements(lshift))
      spec_length=double(stp[j]-str[j])
      bgscale=float(bg_length/spec_length)
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
    
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=1
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=1
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:str,upper:stp,dq:dq_val,dq_bcg:dq_bval}
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
    logprint,'Upper width specified is less than or equal to zero. CONTROL will reset it to default value of 10.'
    upper =10
  endif
  if lower le 0 then begin
    logprint,'Lower width specified is less than or equal to zero. CONTROL will reset it to default value of 10.'
    lower =10
  endif
  stp= fix(centroid+slope*pixel+upper)
  str=fix(centroid+slope*pixel-lower)
  out_frame_up=where(stp gt ny-1)
  out_frame_down=where(str lt 0)
  for j=0,nx-1 do begin
    spectrum_val[j]=total(im[j,str[j]:stp[j]], /NAN)
    noise[j]=total(sigma[j,str[j]:stp[j]]^2, /NAN)
    dq_val[j]=total(dq[j,str[j]:stp[j]], /NAN)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((str[j]+ushift) le ny-1) then begin
        back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift], /NAN)
        noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2, /NAN)
        dq_bval[j]+=total(dq[j,str[j]+ushift:stp[j]+ushift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
      endelse
      if ((stp[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift], /NAN)
        noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2, /NAN)
        dq_bval[j]+=total(dq[j,str[j]-lshift:stp[j]-lshift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse  
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if ((centroid[j]+ushift[k]) le ny-2) then begin
          back_u[j]=total(im[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          noise_u[j]=total(sigma[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if ((centroid[j]-lshift[k]) ge 1) then begin
          back_u[j]=total(im[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          noise_u[j]=total(sigma[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
           dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=double(n_elements(ushift)+n_elements(lshift))
      spec_length=double(stp[j]-str[j])
      bgscale=float(bg_length/spec_length)
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=1
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=1
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:str,upper:stp,dq:dq_val,dq_bcg:dq_bval}
  return,spectrum
end

;extraction with pixel dependend centroid and pixel dependend width

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
  if n_elements(centroid) ne nx then begin
    print,'Input dimension error, centroid has different dimension than image array '
    err = '%extract_varsum: Input diamension error'
  endif
  for j=1,nx-2 do begin ;valid for simulated data change to: j=0,nx-1 for actual data
    k=centroid[j]
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
      if j ge 4 and j le n_elements(m_upper)-5 then  m_upper[j]= median(m_upper[j-4:j+4]) else if j lt 4 then m_upper[j]= median(m_upper[j:j+8]) else m_upper[j]= median(m_upper[j-8:j])
    endif
    if (ll_thresh ge m_lower[j] le lu_thresh) then begin
      if j ge 4 and j le n_elements(m_lower)-5 then  m_lower[j]= median(m_lower[j-4:j+4]) else if  j lt 4 then m_lower[j]= median(m_lower[j:j+8]) else m_lower[j]= median(m_lower[j-8:j])
    endif
  endfor
  pixel=indgen(nx)
  degree=1
  mn_upper=dblarr(n_elements(m_upper))
  mn_lower=dblarr(n_elements(m_lower))
  fit=poly_fit(pix,m_upper,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma_fit,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_upper+=fit[i]*pix^i
  
  fit=poly_fit(pix,m_lower,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma_fit,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_lower+=fit[i]*pix^i


  n_upper=fix(interpol(mn_upper,pix,pixel))+1
  n_lower=fix(interpol(mn_lower,pix,pixel))
 
  for j=0, nx-1 do begin
    spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]], /NAN)
    noise[j]=total(sigma[j,n_lower[j]:n_upper[j]]^2, /NAN)
    dq_val[j]=total(dq[j,n_lower[j]:n_upper[j]], /NAN)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((n_upper[j]+ushift) le ny-1) then begin
        back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
        noise_u[j]=total(sigma[j,n_lower[j]+ushift:n_upper[j]+ushift]^2, /NAN)
        dq_bval[j]+=total(dq[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' ,logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
      endelse
      if ((n_lower[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
        noise_l[j]=total(sigma[j,n_lower[j]-lshift:n_upper[j]-lshift]^2, /NAN)
        dq_bval[j]+=total(dq[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
        message,'Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if (centroid[j]+ushift[k] le ny-2) then begin
          back_u[j]=total(im[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          noise_u[j]=total(sigma[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,'Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if (centroid[j]-lshift[k] ge 1) then begin
          back_u[j]=total(im[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          noise_u[j]=total(sigma[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
          message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=double(n_elements(ushift)+n_elements(lshift))
      spec_length=double(n_upper[j]-n_lower[j])
      bgscale=float(bg_length/spec_length)
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=1
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=1
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:n_lower,upper:n_upper,dq:dq_val,dq_bcg:dq_bval}

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
    logprint,'Input dimension error, centroid has different dimension than image array',logonly = logonly
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
      if j ge 4 and j le n_elements(m_upper)-5 then  m_upper[j]= mean(m_upper[j-4:j+4]) else if j lt 4 then m_upper[j]= mean(m_upper[j:j+8]) else m_upper[j]= mean(m_upper[j-8:j])
    endif
    if (ll_thresh ge m_lower[j] le lu_thresh) then begin
      if j ge 4 and j le n_elements(m_lower)-5 then  m_lower[j]= mean(m_lower[j-4:j+4]) else if  j lt 4 then m_lower[j]= mean(m_lower[j:j+8]) else m_lower[j]= mean(m_lower[j-8:j])
    endif
  endfor
  
  pixel=indgen(nx)
  degree=1
  mn_upper=dblarr(n_elements(m_upper))
  mn_lower=dblarr(n_elements(m_lower))
  fit=poly_fit(pix,m_upper,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_upper+=fit[i]*pix^i
  
  fit=poly_fit(pix,m_lower,1,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do mn_lower+=fit[i]*pix^i


  n_upper=fix(interpol(mn_upper,pix,pixel))+1
  n_lower=fix(interpol(mn_lower,pix,pixel))
  n_upper=fix(interpol(m_upper,pix,pixel))
  n_lower=fix(interpol(m_lower,pix,pixel))
  diffu=intarr(n_elements(n_upper))
  diffl=intarr(n_elements(n_upper))

  for j=0, nx-1 do begin
    spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]], /NAN)
    noise[j]=total(sigmaim[j,n_lower[j]:n_upper[j]]^2, /NAN)
    dq_val[j]=total(dq[j,n_lower[j]:n_upper[j]])
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((n_upper[j]+ushift) le ny) then begin
        back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
        noise_u[j]=total(sigmaim[j,n_lower[j]+ushift:n_upper[j]+ushift]^2, /NAN)
        dq_bval[j]+=total(dq[j,n_lower[j]+ushift:n_upper[j]+ushift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' ,logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
      endelse
      if ((n_lower[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
        noise_l[j]=total(sigmaim[j,n_lower[j]-lshift:n_upper[j]-lshift]^2, /NAN)
        dq_bval[j]+=total(dq[j,n_lower[j]-lshift:n_upper[j]-lshift], /NAN)
      endif else begin
        errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
        message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if (centroid[j]+ushift[k] le ny-2) then begin
          back_u[j]=total(im[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          noise_u[j]=total(sigmaim[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]+ushift[k]:centroid[j]+ushift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds',logonly=logonly
          message,'Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if (centroid[j]-lshift[k] gt 0) then begin
          back_u[j]=total(im[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          noise_u[j]=total(sigmaim[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
          dq_bval[j]+=total(dq[j,centroid[j]-lshift[k]:centroid[j]-lshift[k]+1], /NAN)
        endif else begin
          errorlog,'CONTROL_TRACE: background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds',logonly=logonly
          message,'Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=double(n_elements(ushift)+n_elements(lshift))
      spec_length=double(n_upper[j]-n_lower[j])
      bgscale=float(bg_length/spec_length)
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  prb=where(dq_val ge 1)
  if total(prb) ne -1 then dq_val(prb)=1
  undefine,prb
  prbb=where(dq_bval ge 1)
  if total(prbb) ne -1 then dq_bval(prbb)=1
  undefine,prbb
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:n_lower,upper:n_upper,dq:dq_val,dq_bcg:dq_bval}
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
    logprint,'CONTROL: No temporary file path found. Temporary file directory created in data directory'
  endif else inter_path=detectos(inter_path)
   
im=in_image.data
error=in_image.error
hdr=in_image.header
dq=in_image.dq
;im=mrdfits('D:\simulator_output\kelt9\orbit00\raw_image00000.fits',0,hdr)
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
    logprint,'CONTROL_TRACE: Centroid value specified in the parameter file is close to the edge of the frame.'
    logprint,'Press q to skip the trace and extraction for current spectrum. Press any key to continue trace and extraction for the current spectrum.'
    R = GET_KBRD()
    if R eq 'q' then begin
      logprint,'CONTROL: Skipping the trace and extraction for current spectrum as requested by the user.'
      return,0
    endif else begin
      if cent_value le border then low_range = 0 else if cent_value ge ny-border+1 then high_range = ny-1
    endelse
  endif   
  implot=im[*,low_range:high_range]
  window,xsize=1800,ysize=600
  cgimage,implot,/AXES,/Save,yrange=[centroid[0]-15,centroid[0]+15],xrange=[pixel[0],pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2.5
  cgoplot,pixel,centroid+0.5,color='red';,yrange=[245,265]
  cgLegend, Title=['Centroid'], LineStyle=[0], SymSize=2, Color=['red'], Location=[pixel[nx-300],centroid[0]+14],/DATA ,thick=1.5,Length=0.05, VSpace=2.0, /Box,/Background, BG_Color='white'
  write_png,inter_path+filename+'_centroid.png',TVRD(/TRUE)
endif else begin
  logprint,'CONTROL_TRACE: centroid not specified. Calculating on the go'
  t_cent=dblarr(nx-8+1)
  for i=0, nx-8 do begin
    c_sum=total(im[i:i+7,*],1)
    m=max(c_sum,loc)
    t_cent[k]=loc
    k++
  endfor
  trace_cen=dblarr(n_elements(t_cent))
  pix=indgen(n_elements(t_cent))
  pixel=indgen(nx)
  fit=poly_fit(pix,t_cent,degree,chisq=chisq,/DOUBLE, MEASURE_ERRORS=merror,SIGMA=sigma,STATUS=status, YBAND=yband,YERROR=yerror,YFIT=yfit)
  for i=0,degree do trace_cen+=fit[i]*pix^i
  centroid_full=interpol(trace_cen,pix,pixel)
  centroid=fix(centroid_full)
  y=indgen(ny)
  implot=im[*,centroid[0]-15:centroid[0]+15]
  window,xsize=1800,ysize=600
  cgimage,implot,/AXES,/Save,yrange=[centroid[0]-15,centroid[0]+15],xrange=[pixel[0],pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2.5
  cgoplot,pixel,centroid+0.5,color='red';,yrange=[245,265]
  cgLegend, Title=['Centroid'], LineStyle=[0], SymSize=2, Color=['red'], Location=[pixel[nx-300],centroid[0]+14],/DATA ,thick=1.5,Length=0.05, VSpace=2.0, /Box,/Background, BG_Color='white'
  write_png,inter_path+filename+'centroid.png',TVRD(/TRUE)
endelse

;background stuff
if tag_exist(infile,'background_trace') eq 0 then begin
  errorlog,'CONTROL_TRACE: Background trace information is required. Please re-run the simulator with the same',logonly=1
  message,'CONTROL_TRACE: Background trace information is required. Please re-run the simulator with the same'
endif  
if (file_test(infile.background_trace)eq 0) then begin
  logprint,'CONTROL_TRACE: Background trace input is pixel deviation from centroid'
  str_bcg=strsplit(infile.background_trace,',',/EXTRACT)
  if n_elements(str_bcg) eq 1 then begin
    ushift = fix(infile.background_trace)
    lshift = fix(infile.background_trace) 
    if ushift eq 0 then begin
      logprint,'Shift value for background trace is 0 or not specified'  
    endif  
    multiple=0
  endif else begin
    shift_val=fix(str_bcg)
    ush_loc=where(shift_val ge 0,/NULL)
    lsh_loc=where(shift_val lt 0,/NULL)
    if ush_loc eq !NULL then undefine,ushift else ushift=shift_val[ush_loc]
    if lsh_loc eq !NULL then undefine,lshift else lshift=-1*shift_val[lsh_loc]
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
endif else begin
  if datatype(lshift,2) eq 0 then lshift = 'undefined'
  if datatype(ushift,2) eq 0 then ushift = 'undefined'
endelse
back_trace={type:multiple,shift_upper:ushift,shift_lower:lshift}
;default extraction of box (similar to COS extraction)
;impliment slope change with extraction cut
if tag_exist(infile,'slope') eq 1 then slope=float(infile.slope)
if tag_exist(infile,'width') eq 1 then width=fix(infile.width)
if tag_exist(infile,'upper') eq 1 then upper=fix(infile.upper)
if tag_exist(infile,'lower') eq 1 then lower=fix(infile.lower)
if tag_exist(infile,'threshold') eq 1 then threshold=fix(infile.threshold)


logprint,'Centroid value used is '+strtrim(string(centroid[0]),2)+' at the edge of the spectrum.'
if datatype(slope,2) ne 0 then logprint,'Slope value used is '+strtrim(string(slope),2)+'.'
if datatype(width,2) ne 0 then if width gt 0 then logprint,'Width value used for extraction type simple is '+strtrim(string(width),2)+'.'
if datatype(fix(infile.upper),2) ne 0 then logprint,'Upper width value used for extraction type fixed is '+strtrim(string(fix(infile.upper)),2)+'.'
if datatype(fix(infile.lower),2) ne 0 then logprint,'Lower width value used for extraction type fixed is '+strtrim(string(fix(infile.lower)),2)+'.'
if datatype(double(infile.threshold),2) ne 0 then logprint,'Threshold value used for extraction type variable/function is '+strtrim(string(double(infile.threshold)),2)+'.'

if datatype(slope,2) eq 0 then begin
  logprint,'Slope value not found. Using default value of -9.76e-4.'
  slope=-9.76e-4
endif
if datatype(width,2) eq 0 then begin
  logprint,'Width value not found. Using default value of 10.'
  width=10
endif
def_spectrum = extract_box(im,error,dq,centroid,slope,width,back_trace)

;extraction
case trace_type of
    'simple'  : begin
                  extraction = extract_box(im,error,dq,centroid,slope,width,back_trace)
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
                end  
    'variable': begin
                  if tag_exist(infile,'threshold') eq 1 then threshold=double(infile.threshold) else begin
                    logprint,'Threshold value not found. Using default value of 0.01'
                    threshold=0.01
                  endelse
                  extraction = extract_varsum(im,error,dq,centroid,threshold,back_trace)
                end  
    'function': begin
                  if tag_exist(infile,'threshold') eq 1 then threshold=double(infile.threshold) else begin
                    logprint,'Threshold value not found. Using default value of 0.01'
                    threshold=0.01
                  endelse
                  extraction = extract_func(im,error,dq,centroid,threshold,back_trace)
                end
    else :      begin
                  logprint,'CONTROL_TRACE: Invalid type input for trace: Please recheck your input',logonly=logonly
                  message,' Invalid type input for trace: Please recheck your input'
                end  
endcase
implot=im[*,centroid[0]-15:centroid[0]+15]
window,xsize=1800,ysize=600
cgimage,implot,/AXES,/Save,yrange=[centroid[0]-15,centroid[0]+15],xrange=[pixel[0],pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2.5
cgoplot,pixel,extraction.lower,color='red';,yrange=[245,265]
cgoplot,pixel,extraction.upper+1,color='red';,yrange=[245,265]
cgLegend, Title=['Trace region'], LineStyle=[0], SymSize=2, Color=['red'], Location=[pixel[nx-400],centroid[0]+14],/DATA ,thick=1.5,Length=0.05, VSpace=2.0, /Box,/Background, BG_Color='white'
write_png,inter_path+filename+'_trace.png',TVRD(/TRUE)  
ext_data=extraction.data
ext_error=extraction.error
ext_bck=extraction.bck
ext_bck_error=extraction.bck_error
ext_dq=extraction.dq
ext_bg_dq=extraction.dq_bcg
;comparison to make sure extraction is proper, comparison with default 
;plt=plot(ext_data)

spectrum={data:ext_data,error:ext_error,header:hdr,background:ext_bck,bck_error:ext_bck_error,dq:ext_dq,dq_bg:ext_bg_dq} 
return,spectrum  
end