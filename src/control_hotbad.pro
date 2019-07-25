; NAME:
;      CONTROL_HOTBAD
;
; PURPOSE:
;      Carried out hot and bad pixel correction for CCD frames based on a hot/bad pixel map.
; 
; CALLING SEQUENCE:
;      out= CONTROL_HOTBAD,input,map
;
; INPUTS:
;      input  = The frame on which hot/bad pixels are to be corrected.
;      hbmap    =  map of cosmetic defects on the CCD as a 2-D array with the same sizes as the input CCD frame byte-type pixels (8-bit) and values of 0B for good pixels and 1B for bad onces.
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
;      modified 21.02.2019 by A. G. Sreejith
;#######################################################################

pro CONTROL_HOTBAD,input_file,hbmap,out,type

  err = ''
  if N_params() LT 3 then begin             ;Need at least 4 parameters
    logprint,'Syntax - CONTROL_HOTBAD,input_file,hbmap,out,type',logonly = logonly
    message,'Syntax - CONTROL_HOTBAD,input_file,hbmap,out,type'
    err = '%control_hotbad: Insufficient number of parameters'
    return
  endif
  if not keyword_defined(type) then type ='interpolate' 
  input=mrdfits(input_file,0,hdr,/SILENT)
  ; checking the inputs
  elem= (size(input))[0]
  
  if (elem eq 2) then begin 
    nx = (size(input))[1]
    ny = (size(input))[2]
    mnx= (size(hbmap))[1]
    mny=  (size(hbmap))[1]
    if not keyword_defined(hbmap) then hbmap = bytarr(nx, ny) + 1
  endif else begin
    logprint,'CONTROL_HOTBAD input file dimension error',logonly = logonly
    message,'CONTROL_HOTBAD input file dimension error'
    
  endelse
  if nx ne mnx then begin
    logprint,'CONTROL: Input file error, Input file:'+input_file+' and hot and bad pixel maps have different columns',logonly = logonly
    message,'CONTROL: Input file error, Input file:'+input_file+' and hot and bad pixel maps have different columns'
    ;sxaddpar, hdr, 'COLMNERR', 1. ;Hot and bad pixel correction flag
  endif
  col_data = dblarr(nx,ny)
  row_data = dblarr(nx,ny)
     
  ycut1=SXPAR( hdr, 'YCUT1')
  ycut2=SXPAR( hdr, 'YCUT2')
  ylen=ycut2-ycut1+1
  if (ylen ne ny) then hbmap_nw=hbmap[*,ycut1:ycut1+ny-1] else hbmap_nw=hbmap[*,ycut1:ycut2]
  ;execulte bad pixel removal row wise takes care of bad columns
  for i = 0, ny-1 do begin
    bv = where(hbmap_nw[0:nx-1, i] eq 1)
    gv = where(hbmap_nw[0:nx-1, i] eq 0)
    datav = input[0:nx-1, i]
    if gv[0] ne -1 then begin
      case type of
        'interpolate': if (bv ne [-1]) then datav[bv] = interpol(datav[gv], gv, bv)
        'average'    : if (bv ne [-1]) then datav[bv] = (datav[bv-1]+datav[bv+1])/2 
        else : logprint,'Invalid type input for combine: Please recheck your input'
      endcase
    endif  
      col_data[*,i]=datav
    endfor

  ;execulte bad pixel removal column wise takes care of bad rows
  for i = 0, nx-1 do begin
    bv = where(hbmap_nw[i,0:ny-1] eq 1)
    gv = where(hbmap_nw[i,0:ny-1] eq 0)
    datav = col_data[i,0:ny-1]
    if gv[0] ne -1 then begin
      case type of
        'interpolate': if (bv ne [-1]) then datav[bv] = interpol(datav[gv], gv, bv)
        'average'    : if (bv ne [-1]) then datav[bv] = (datav[bv-1]+datav[bv+1])/2 
        else : logprint,'Invalid type input for combine: Please recheck your input'
      endcase
    endif  
    row_data[i,*]=datav
  endfor

  out={data:row_data,header:hdr}
 return
end