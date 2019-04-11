function mode, array,diamension=diamension,NAN = nan

  ; Calculates the MODE (value with the maximum frequency distribution)
  ; of an array. Works ONLY with integer data.

  On_Error, 2

  ; Check for arguments.
  if n_elements(array) eq 0 then message, 'Must pass an array argument.'
  
  nxy=size(array)
  
  if diamension gt nxy[0] then message, 'Diamension miss match with array.'
  if datatype(array) gt t then message,'Invalid datatype for array'
  
  n=nxy[1:nxy[0]]
  ndia=where(n ne diamension)
  if total(ndia) ne -1 then mode_arr=make_array(ndia,TYPE=datatype(array))
  nele=nxy(diamension)
  if diamensions gt 3 then message, 'Larger than third diamensional arrays not supported at the moment.'
  
  if nxy[0] eq 3 then begin
    for i=0,ndia[0] do begin
      for j=0,ndia[1] do begin
        if diamension eq 1 then begin
          input_array=array[*,i,j]
          mode_arr[i,j]=mode_val(input_array)
        endif else if diamension eq 2 then begin
          input_array=array[i,*,j]
          mode_arr[i,j]=mode_val(input_array)
        endif else begin
          input_array=array[i,j,*]
          mode_arr[i,j]=mode_val(input_array)
        endelse
      endfor
    endfor
  endif else if nxy[0] eq 2 then begin
    for i=0,ndia[0] do begin
      if diamension eq 1 then begin
        input_array=array[*,i]
        mode_arr[i]=mode_val(input_array)
      endif else if diamension eq 2 then begin
        input_array=array[i,*]
        mode_arr[i]=mode_val(input_array)
      endif
    endfor
  endif else if nxy[0] eq 1 then mode_arr=mode_val(array)
       
  return, mode_arr
end