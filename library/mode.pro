function mode, array,dimension=dimension,NAN = nan

  ; Calculates the MODE (value with the maximum frequency distribution)
  ; of an array. Works ONLY with integer data.

  On_Error, 2
  
  ; Check for arguments.
  if n_elements(array) eq 0 then message, 'Must pass an array argument.'
  
  nxy=size(array)
  if keyword_defined(dimension) eq 0 then  begin 
    mode_arr=mode_val(array)
    return,mode_arr
  endif
  if dimension gt nxy[0] then message, 'dimension miss match with array.'
  ;if datatype(array) gt t then message,'Invalid datatype for array'
  place_holder=nxy[0]
  n=nxy[1:place_holder-1]
  ndia=n[where(n ne dimension)]

  if total(ndia) ne -1 then mode_arr=make_array(ndia,TYPE=datatype(array,2))

  nele=nxy(dimension)
  if dimension gt 3 then message, 'Larger than third dimensional arrays not supported at the moment.'
 
  if place_holder eq 3 then begin
    for i=0,ndia[0]-1 do begin
      for j=0,ndia[1]-1 do begin
        if dimension eq 1 then begin
          input_array=array[*,i,j]
          mode_arr[i,j]=mode_val(input_array)
        endif else if dimension eq 2 then begin
          input_array=array[i,*,j]
          mode_arr[i,j]=mode_val(input_array)
        endif else begin
          input_array=array[i,j,*]
          mode_arr[i,j]=mode_val(input_array)

        endelse
      endfor
    endfor
  endif else if place_holder eq 2 then begin
    for i=0,ndia[0] do begin
      if dimension eq 1 then begin
        input_array=array[*,i]
        mode_arr[i]=mode_val(input_array)
      endif else if dimension eq 2 then begin
        input_array=array[i,*]
        mode_arr[i]=mode_val(input_array)
      endif
    endfor
  endif else if place_holder eq 1 then begin
    mode_arr=mode_val(array)
  
  endif

  return, mode_arr
end