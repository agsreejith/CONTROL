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
; sxaddpar,hdr_nw,'OBS_TIME',sxpar(hdr,'OBS_TIME'),'Time of observation'
; sxaddpar,hdr_nw,'STAR_TEM',STAR_TEM,'Stellar temperature'
; sxaddpar,hdr_nw,'STAR_RAD=        2.36000000000 /Stellar radius                                   
; sxaddpar,hdr_nw,'STAR_MAG=        7.55000000000 /Stellar magnitude
; sxaddpar,hdr_nw,'RA_TRG  =        307.859791700 /Right ascension     
; sxaddpar,hdr_nw,'DEC_TRG =        39.9388330000 /Declination
; sxaddpar,hdr_nw,'TYPE_OF_= 'Default: From STIS backgrounds' /Type of background                   SLIT_POS=       0.000000000000 /Position of star on slit wrt slit center
;  SPEC_RES=             0.800000 /Spectral resolution                              RT_NOISE=        12.2500000000 /Read out noise of observation
;  EXP_TIME=                  300 /Exposure time                                    MID_TIME=      2457095.6857200 /Mid transit time of observation
;  FILETYPE= 'OBJECT  '           /Type of observation                              TARGNAME= 'TEST    '           /Target name of observation
;  YCUT1   =                    0 /Bottom of science box extraction                 YCUT2   =                  514 /Top of science box extraction
;  WAVESHFT=       0.000000000000 /Wavelength shift in Angstroms                    IDLVER  =              8.50000 /Version of IDL used to produce the frame
;  ORBPRD  =              5400.00 /Orbit period of satellite                        FRAORB  =       0.000000000000 /Time as a fraction of orbital period
;  BCFLG   =                    1 /BIAS CORRECTION FLAG                            
   sxdelpar,hdr_nw,'PCOUNT'
   sxdelpar,hdr_nw,'GCOUNT'                                           ;CRFLG   =                    1 /COSMIC RAY CORRECTION FLAG
  ;DCFLG   =                    1 /DARK CORRECTION FLAG                             FCFLG   =                    1 /FLAT CORRECTION FLAG

return,hdr_nw
end