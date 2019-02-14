; NAME:
;      CONTROL_HOTBAD
;
; PURPOSE:
;      Carried out hot and bad pixel correction for CCD frames based on a hot/bad pixel map. This module also finds new hot or bad pixels if any.
; 
; CALLING SEQUENCE:
;      out= CONTROL_FLAT_COMBINE,input,map,/append
;
; INPUTS:
;      input  = The frame on which hot/bad pixels are to be corrected.
;      hbmap    = The location of map of cosmetic defects on the CCD as a FITS file with the same sizes as the input CCD frame byte-type pixels (8-bit) and values of 1B for good pixels and 0B for bad onces. 
;
; OPTIONAL OUTPUT:
;     type = 
;   
; OUTPUT:
;      out    = hot/bad pixel corrected frame.
; REQUIRES:
;     map of hot/bad pixel
;
; PROCEDURE:
;      Master flat creator for CUTE;
; MODIFICATION HISTORY:
;      created 9.01.2019 by A. G. Sreejith
;      modified 15.01.2019 by A. G. Sreejith
;#######################################################################

pro CONTROL_HOTBAD,input_file,hbmap,out,type

  err = ''
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    print,'Syntax - CONTROL_HOTBAD,input,hbmap,append=append
    err = '%control_hotbad: Insufficient number of parameters'
    return
  endif
  if not keyword_defined(type) then type ='interpolate' 
  input=mrdfits(input_file,0,hdr)
  ; checking the inputs
  elem= (size(input))[0]
  
  if (elem eq 2) then begin 
    nx = (size(input))[1]
    ny = (size(input))[2]
    if not keyword_defined(hbmap) then hbmap = bytarr(nx, ny) + 1
  endif else begin
    print,'CONTROL_HOTBAD input file error'
    return
  endelse
  if nx ne size(hbmap)[1] then begin
    print,'CONTROL: Input file error, Input file:'+input file' and hot and bad pixel maps have different columns'
    sxaddpar, hdr, 'COLMNERR', 1. ;Hot and bad pixel correction flag
    return
  endif
  out_data = dblarr(nx,ny) 
  ycut1=SXPAR( hdr, 'YCUT1')
  ycut2=SXPAR( hdr, 'YCUT2')
  ylen=ycut2-ycut1+1
  if (ylen ne ny) then hbmap_nw=[*,ycut1:ycut1+ny-1] else hbmap_nw=[*,ycut1:ycut2]
  for i = 0, ny-1 do begin
    bv = where(hbmap_nw[0:nx-1, i] eq 0)
    gv = where(hbmap_nw[0:nx-1, i] eq 1)
    datav = input[0:nx-1, i]
    ;check for bad rows 
    
    
    case type of
      'interpolate': if (bv ne [-1]) then datav[bv] = interpol(datav[gv], gv, bv)
      'average'    : if ((bv ne [-1]) and (rflg ne 1)) then datav[bv] = (datav[bv-1]+datav[bv+1])/2 else datav[bv] = hb_avg
      else : print,'Invalid type input for combine: Please recheck your input'
    endcase
    out_data[*,i]=datav
  endfor
 out={data:out_data,header:hdr}
 return
end