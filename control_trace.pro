
;Simplest extraction a box with fixed width, centroid and slope
function extract_box,im,sigma,centroid,slope,width,back_trace
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
  ushift=fix(back_trace.shift_upper)
  lshift=fix(back_trace.shift_lower)
  
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
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
  stp= fix(centroid+slope*pixel+width)
  str=fix(centroid+slope*pixel-width)
 
  for j=0,nx-1 do begin
    spectrum_val[j]=total(im[j,str[j]:stp[j]])
    noise[j]=total(sigma[j,str[j]:stp[j]]^2)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((stp+ushift) le ny-1) then begin
        back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift])
        noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
      endelse
      if ((str-lshift) ge 0) then begin
        back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift])
        noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if ((centroid+ushift[k]) le ny-2) then begin
          back_u[j]=total(im[j,centroid+ushift[k]:centroid+ushift[k]+1])
          noise_u[j]=total(sigma[j,centroid+ushift[k]:centroid+ushift[k]+1])
        endif else begin
          logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if ((centroid-lshift[k]) ge 1) then begin
          back_l[j]=total(im[j,centroid-lshift[k]:centroid-lshift[k]-1])
          noise_l[j]=total(sigma[j,centroid-lshift[k]:centroid-lshift[k]-1])
        endif else begin
          logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=n_elements(ushift)+n_elements(lshift)
      spec_length=stp-str
      bgscale=bg_length/spec_length
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
    
  endfor
  error=sqrt(noise)

  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:str,upper:stp}
  return,spectrum
end

;extraction with pixel dependednd centroid and fixed width
function extract_sum,im,sigma,centroid,upper,lower,back_trace
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
  ushift=fix(back_trace.shift_upper)
  lshift=fix(back_trace.shift_lower)
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  pixel=indgen(nx)
  if n_elements(centroid) ne nx then begin
    message,'Input diamension error, centroid has different diamension than image array '
    err = '%extract_sum: Input diamension error'
  endif
  stp= fix(centroid+slope*pixel+upper)
  str=fix(centroid+slope*pixel-lower)

  for j=0,nx-1 do begin
    spectrum_val[j]=total(im[j,str[j]:stp[j]])
    noise[j]=total(sigma[j,str[j]:stp[j]]^2)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((str+ushift) le ny-1) then begin
        back_u[j]=total(im[j,str[j]+ushift:stp[j]+ushift])
        noise_u[j]=total(sigma[j,str[j]+ushift:stp[j]+ushift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
      endelse
      if ((stp-lshift) ge 0) then begin
        back_l[j]=total(im[j,str[j]-lshift:stp[j]-lshift])
        noise_l[j]=total(sigma[j,str[j]-lshift:stp[j]-lshift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse  
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if ((centroid+ushift[k]) le ny-2) then begin
          back_u[j]=total(im[j,centroid+ushift[k]:centroid+ushift[k]+1])
          noise_u[j]=total(sigma[j,centroid+ushift[k]:centroid+ushift[k]+1])
        endif else begin
          logprint,'CONTROL_TRACE: background offset specified is outside the CCD area',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds'
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if ((centroid-lshift[k]) ge 1) then begin
          back_u[j]=total(im[j,centroid-lshift[k]:centroid-lshift[k]-1])
          noise_u[j]=total(sigma[j,centroid-lshift[k]:centroid-lshift[k]-1])
        endif else begin
          logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
          message,' Background offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=n_elements(ushift)+n_elements(lshift)
      spec_length=stp-str
      bgscale=bg_length/spec_length
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:centroid-lower,upper:centroid+upper}
  return,spectrum
end

;extraction with pixel dependend centroid and pixel dependend width

function extract_varsum,im,sigma,centroid,threshold,back_trace
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
  ushift=fix(back_trace.shift_upper)
  lshift=fix(back_trace.shift_lower)
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  upper=intarr(nx)
  lower=intarr(nx)
  if n_elements(centroid) ne nx then begin
    print,'Input diamension error, centroid has different diamension than image array '
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
    spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]])
    noise[j]=total(sigma[j,n_lower[j]:n_upper[j]]^2)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((n_upper[j]+ushift) le ny-1) then begin
        back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift])
        noise_u[j]=total(sigma[j,n_lower[j]+ushift:n_upper[j]+ushift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
      endelse
      if ((n_lower[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift])
        noise_l[j]=total(sigma[j,n_lower[j]-lshift:n_upper[j]-lshift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
        message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if (centroid+ushift[k] le ny-2) then begin
          back_u[j]=total(im[j,centroid+ushift[k]:centroid+ushift[k]+1])
          noise_u[j]=total(sigma[j,centroid+ushift[k]:centroid+ushift[k]+1])
        endif else begin
          logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
          message,'Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if (centroid-lshift[k] ge 1) then begin
          back_u[j]=total(im[j,centroid-lshift[k]:centroid-lshift[k]-1])
          noise_u[j]=total(sigma[j,centroid-lshift[k]:centroid-lshift[k]-1])
        endif else begin
          logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
          message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=n_elements(ushift)+n_elements(lshift)
      spec_length=n_upper[j]-n_lower[j]
      bgscale=bg_length/spec_length
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:n_lower,upper:n_upper}

  return,spectrum
end

;extraction with pixel dependend centroid and a cross dispersion function
function extract_func,im,sigma,centroid,threshold,back_trace
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
  nxy2=size(sigma)
  if nxy2[0] ne 2 then begin
    print,'Input file error, the error file is not 2-D'
    err = '%extract_func: Input file error'
  endif
  nx=nxy[1]
  ny=nxy[2]
  multiple=back_trace.type
  ushift=fix(back_trace.shift_upper)
  lshift=fix(back_trace.shift_lower)
  spectrum_val=dblarr(nx)
  noise=dblarr(nx)
  back_u=dblarr(nx)
  noise_u=dblarr(nx)
  back_l=dblarr(nx)
  noise_l=dblarr(nx)
  upper=intarr(nx)
  lower=intarr(nx)
  background=dblarr(nx)
  backg_error=dblarr(nx)
  if n_elements(centroid) ne nx then begin
    logprint,'Input diamension error, centroid has different diamension than image array',logonly = logonly
    message,'Input diamension error, centroid has different diamension than image array'
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
    spectrum_val[j]=total(im[j,n_lower[j]:n_upper[j]])
    noise[j]=total(sigma[j,n_lower[j]:n_upper[j]]^2)
    ;for background use the input from lookup table as well
    if multiple eq 0 then begin
      if ((n_upper[j]+ushift) le ny) then begin
        back_u[j]=total(im[j,n_lower[j]+ushift:n_upper[j]+ushift])
        noise_u[j]=total(sigma[j,n_lower[j]+ushift:n_upper[j]+ushift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
        message,' Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
      endelse
      if ((n_lower[j]-lshift) ge 0) then begin
        back_l[j]=total(im[j,n_lower[j]-lshift:n_upper[j]-lshift])
        noise_l[j]=total(sigma[j,n_lower[j]-lshift:n_upper[j]-lshift]^2)
      endif else begin
        logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
        message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
      endelse
      background[j]=(back_u[j]+back_l[j])/2
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/4)
    endif else begin
      for k=0,n_elements(ushift)-1 do begin
        if (centroid+ushift[k] le ny-2) then begin
          back_u[j]=total(im[j,centroid+ushift[k]:centroid+ushift[k]+1])
          noise_u[j]=total(sigma[j,centroid+ushift[k]:centroid+ushift[k]+1])
        endif else begin
          logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
          message,'Background offset specified is outside the CCD area. Please re-enter a new upper offset for backgrounds' 
        endelse
      endfor
      for k=0,n_elements(lshift)-1 do begin
        if (centroid-lshift[k] gt 0) then begin
          back_u[j]=total(im[j,centroid-lshift[k]:centroid-lshift[k]-1])
          noise_u[j]=total(sigma[j,centroid-lshift[k]:centroid-lshift[k]-1])
        endif else begin
          logprint,'CONTROL_TRACE: backhround offset specified is outside the CCD area',logonly=logonly
          message,'Backround offset specified is outside the CCD area. Please re-enter a new lower offset for backgrounds'
        endelse
      endfor
      bg_length=n_elements(ushift)+n_elements(lshift)
      spec_length=n_upper[j]-n_lower[j]
      bgscale=bg_length/spec_length
      background[j]=(back_u[j]+back_l[j])/bgscale
      backg_error[j]=(sqrt(noise_u[j]+noise_l[j])/bgscale^2)
    endelse
  endfor
  error=sqrt(noise)
  spectrum={data:spectrum_val,error:error,bck:background,bck_error:backg_error,centroid:centroid,lower:n_lower,upper:n_upper}
  return,spectrum
end

function control_trace,in_image,infile,trace_type,filename
  if N_params() LT 4 then begin             ;Need at least 4 parameters
    print,'Syntax - control_trace,in_image,infile,trace_type,,filename
    err = '%CONTROL_TRACE: Insufficient number of parameters'
    
  endif
  if keyword_defined(filename) eq 0 then begin
    logprint,'CONTROL_TRACE: Type input for trace not defined uisng default: variable'
    trace_type='variable'
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
;im=mrdfits('D:\simulator_output\kelt9\orbit00\raw_image00000.fits',0,hdr)
ycut1=SXPAR( hdr, 'YCUT1')
ycut2=SXPAR( hdr, 'YCUT2')

constu1=1
constu2=0
constu3=1
constl1=1
constl2=0
constl3=1

if(tag_exist(infile,'cent_poly_degree') eq 1) then degree=infile.cent_poly_degree else degree=1  ;setting degree of polynomial for centroid default is 1
nxy=size(im)
nx=nxy[1]
ny=nxy[2]
k=0

if(tag_exist(infile,'centroid') eq 1) then begin
  cent_value=double(infile.centroid)
  cent_value=cent_value+ycut1 
  centroid=make_array(nx,value=cent_value) 
  pixel=indgen(nx)
  y=indgen(ny)
  implot=im[*,centroid[0]-15:centroid[0]+15]
  window,xsize=1800,ysize=600
  cgimage,implot,/AXES,/Save,yrange=[centroid[0]-15,centroid[0]+15],xrange=[pixel[0],pixel[nx-1]],xtitle='X pixels',ytitle='Y pixels',charsize=2.5
  cgoplot,pixel,centroid+0.5,color='red';,yrange=[245,265]
  cgLegend, Title=['Centroid'], LineStyle=[0], SymSize=2, Color=['red'], Location=[pixel[nx-300],centroid[0]+14],/DATA ,thick=1.5,Length=0.05, VSpace=2.0, /Box,/Background, BG_Color='white'
  write_png,inter_path+filename+'centroid.png',TVRD(/TRUE)
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
if valid_num(infile.background_trace) eq 1 then begin
  ushift = infile.background_trace
  lshift = infile.background_trace 
  multiple=0
endif else begin  
  bt_file=detectos(infile.background_trace)
  readcol,bt_file,star,ushift,lshift,u_offset,l_offset,F='A,I,I,D,D'
  ;fnd the target
  target=SXPAR( hdr, 'TARGNAME')
  s=STRCOMPRESS(target, /REMOVE_ALL)
  bg_loc=WHERE(STRMATCH(star, s, /FOLD_CASE) EQ 1)

  if (n_elements(bg_loc) eq 1) then multiple=0 else multiple=bg_loc
  if ycut1 le 171 then begin
    ushift=constu1*u_offset+ushift
    lshift=constl1*l_offset+lshift
  endif else if ycut1 le 343 then begin
    ushift=constu2*u_offset+ushift
    lshift=constl2*l_offset+lshift
  endif else begin
    ushift=constu3*u_offset+ushift
    lshift=constl3*l_offset+lshift
  endelse
endelse  
back_trace={type:multiple,shift_upper:ushift,shift_lower:lshift}
;default extraction of box (similar to COS extraction)
;impliment slope change with extraction cut
slope=float(infile.slope)
width=fix(infile.width)

def_spectrum = extract_box(im,error,centroid,slope,width,back_trace)
;extraction
case trace_type of
    'simple'  : begin
                  extraction = extract_box(im,error,centroid,slope,width,back_trace)
                end  
    'fixed'   : begin
                  upper=fix(infile.upper)
                  lower=fix(infile.lower)
                  extraction = extract_sum(im,error,centroid,upper,lower,back_trace)
                end  
    'variable': begin
                  threshold=double(infile.threshold)
                  extraction = extract_varsum(im,error,centroid,threshold,back_trace)
                end  
    'function': begin
                  threshold=double(infile.threshold)
                  extraction = extract_func(im,error,centroid,threshold,back_trace)
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
write_png,inter_path+filename+'trace.png',TVRD(/TRUE)  
ext_data=extraction.data
ext_error=extraction.error
ext_bck=extraction.bck
ext_bck_error=extraction.bck_error
;comparison to make sure extraction is proper, comparison with default 

spectrum={data:ext_data,error:ext_error,header:hdr,background:ext_bck,bck_error:ext_bck_error} 
return,spectrum  
end