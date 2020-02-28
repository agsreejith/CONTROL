; NAME:
;      CONTROL_HBPFIND
;
; PURPOSE:
;      Determines hot and bad pixels for CCD frames based on a hot/bad pixel map and update it.
;
; CALLING SEQUENCE:
;      CONTROL_HBPFIND,input,hbmap
;
; INPUTS:
;      input  = The list of frame on which hot/bad pixels are to be found.
;      hbmap  = The location of map of cosmetic defects on the CCD as a FITS 
;               file with the same sizes as the input CCD frame byte-type pixels 
;               (8-bit) and values of 1B for good pixels and 0B for bad onces.
;               This map is updated at the end of the procedure.
;
; REQUIRES:
;     map of hot/bad pixel
;
; PROCEDURE:
;      Procedure to update hot and bad pixel map for CUTE 
;   
;      
;#################################################################################################

pro control_hbpfind,file_list,mask_file
file_list=file_list
mask=mrdfits(mask_file,0,hdr)
length=n_elements(file_list)
control_hotbad,file_list[0],mask,out
out_im=out.data
std=stddev(out_im)
std_loc = where((out_im ge 5*std) or (out_im le 5*std))
for i=0,length do begin
  control_hotbad,file_list[i],mask,out
  out_im=out.data
  std=stddev(out_im)
  loc = where((out_im ge 5*std) or (out_im le 5*std))
  match2,std_loc,loc,similar
  std_loc=similar 
endfor
new_loc=where(std_loc ne -1)
mask[new_loc]=0b
sxaddpar, hdr, 'UPDTIME', SYSTIME(/JULIAN, /UTC) 
writefits,mask_file,mask,hdr ;hotbad file
end
