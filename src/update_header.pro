function update_header,hdr
  hdr_nw=hdr
  sxdelpar,hdr_nw,'XTENSION'
  sxdelpar,hdr_nw,'SIMPLE' 
  sxdelpar,hdr_nw,'BITPIX'
  sxdelpar,hdr_nw,'NAXIS' 
  sxdelpar,hdr_nw,'NAXIS1'
  sxdelpar,hdr_nw,'NAXIS2' 
  sxdelpar,hdr_nw,'EXTEND'
  sxdelpar, hdr_nw,'TFIELDS'
  for i=0, 10 do begin
    sxdelpar, hdr_nw,'TTYPE'+STRTRIM(STRING(i),2)
    sxdelpar, hdr_nw,'TUNIT'+STRTRIM(STRING(i),2)
    sxdelpar, hdr_nw,'TFORM'+STRTRIM(STRING(i),2)
    sxdelpar, hdr_nw,'TDISP'+STRTRIM(STRING(i),2)
    sxdelpar, hdr_nw,'TNULL'+STRTRIM(STRING(i),2)
    sxdelpar, hdr_nw,'TDIM'+STRTRIM(STRING(i),2)
  endfor
  sxdelpar,hdr_nw,'PCOUNT'
  sxdelpar,hdr_nw,'GCOUNT'                                      
  return,hdr_nw
end
