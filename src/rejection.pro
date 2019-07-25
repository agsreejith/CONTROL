function rejection,input,sigma,npix
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    print,'Syntax - REJECTION,input,sigma,npix
    err = '%REJECTION: Insufficient number of parameters'
    return,-1
  endif
  elem= (size(input))[0]
  if (elem eq 2) then begin
    nx = (size(input))[1]
    ny = (size(input))[2]
  endif
  data=dblarr(nx,ny)  
    npix=0
  thresh_up=mean(input,/DOUBLE)+5*stddev(input,/DOUBLE)
  thresh_lw=mean(input,/DOUBLE)-5*stddev(input,/DOUBLE)
  for i = 0, ny-1 do begin
    bv = where((input[0:nx-1, i] gt thresh_up) or (input[0:nx-1, i] lt thresh_lw))
    gv = where(input[0:nx-1, i] le thresh_up) and (input[0:nx-1, i] ge thresh_lw)
    datav = input[0:nx-1, i]
    if (bv ne [-1]) then datav[bv] = interpol(datav[gv], gv, bv)
    npix+=n_elements(bv)
    data[*,i]=datav
  endfor
  return,data
end