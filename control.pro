; NAME:
;      CONTROL
;      
; PURPOSE:
;      
; CALLING SEQUENCE:
;      
;
; INPUTS:
;      
; OUTPUT:
;      
; REQUIRES:
;     
; PROCEDURE:
;      
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;

pro control
print,' _______  _______  __    _  _______  ______    _______  ___'     
print,'|       ||       ||  |  | ||       ||    _ |  |       ||   |    
print,'|       ||   _   ||   |_| ||_     _||   | ||  |   _   ||   |    
print,'|       ||  | |  ||       |  |   |  |   |_||_ |  | |  ||   |    
print,'|      _||  |_|  ||  _    |  |   |  |    __  ||  |_|  ||   |___ 
print,'|     |_ |       || | |   |  |   |  |   |  | ||       ||       |
print,'|_______||_______||_|  |__|  |___|  |___|  |_||_______||_______|
print,'---------------------------------------------------------------'
print,'    CONTROL v0.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE      '
print,'---------------------------------------------------------------' 
print,'                                                                 '                                                     
;######################################################################
;check section
  ;data_path='D:\simulator_output\hd209\'
  hbmask=1
  mask=bytarr(2048,515)
  mask[10,25]=1
;######################################################################
  ;path selection
  if n_elements(data_path) eq 0 then begin 
        CD, Current=data_path 
        CASE StrUpCase(!Version.OS_Family) OF
          'WINDOWS': data_path=data_path+'\' ;WINDOWS
          'UNIX': data_path=data_path+'/'; UNIX. 
        ENDCASE  
  endif
  
  if n_elements(inter_path) eq 0 then begin
    CASE StrUpCase(!Version.OS_Family) OF
          'WINDOWS': inter_path=data_path+'temp\ ' ;WINDOWS
          'UNIX': inter_path=data_path+'temp/ '; UNIX. 
        ENDCASE 
    print,'CONTROL: No intermediate file path found. Temporary file directory created in data directory'
  endif
  
  if n_elements(out_path) eq 0 then begin
        CASE StrUpCase(!Version.OS_Family) OF
          'WINDOWS': out_path=data_path+'output\ ' ;WINDOWS
          'UNIX': out_path=data_path+'output/ '; UNIX. 
        ENDCASE 
    print,'CONTROL: No output file path found. Output file directory created in data directory'
  endif
  ;creates paths if they donr exist
  if (file_test(data_path,/DIRECTORY)eq 0)) then FILE_MKDIR,data_path
  if (file_test(inter_path,/DIRECTORY)eq 0)) then FILE_MKDIR,inter_path
  if (file_test(out_path,/DIRECTORY)eq 0)) then FILE_MKDIR,out_path
  
  
  file_list=file_search(data_path+'*.fits') ;Get file list
  ;unique files?
  if n_elements(file_list eq 0) then begin
    print,'CONTROL: File acess error, No fits file found!!'
    print,'CONTROL: Verify data path.'
    return
  endif
  
  file_names=file_basename(file_list,'.fits') ;Get file names
  
  ;carry out bad and hot pixel corrections
  if(hbmask) then begin 
    for i=0,n_elements(file_list)-1 do begin 
      control_hotbad,file_list[i],mask,out
      hdr=out.header
      out_im=out.data
      sxaddpar, hdr, 'HBFLG', 1. ;Hot and bad pixel correction flag
      writefits,inter_path+file_names[i]+'_hb.fits',out_im,hdr ;hotbad file 
      file_type[i]=hierarch(hdr, 'HIERARCH OF CUTE');Find object description
    endfor
    file_list_hb=file_search(inter_path+'*_hb.fits')
  endif else begin
    file_type[i]=hierarch(h, 'HIERARCH OF CUTE')
    file_list_hb=file_search(data_path+'*.fits')  
  endelse
  
  ;Automatic file classification
  
  ;Compile a list of BIAS frames
  biaslist=file_list_hb[where(strpos(file_type,'BIAS') ge 0)]
  ;bias_file = prefix+inst_mode+'.bias.fits'         ;Master bias filename

  ;Compile a list of DARK frames
  darklist=file_list_hb[where(strpos(file_type,'DARK') ge 0)]
  ;dark_file = prefix+inst_mode+'.dark.fits'         ;Master bias filename

  ;Compile a list of FLAT frames
  flatlist=file_list_hb[where(strpos(file_type,'FLAT') ge 0)]
  ;flat_file=prefix+inst_mode+'.flat.fits'           ;Master flat filename
  ;norm_flat_file=prefix+inst_mode+'.flat.norm.fits' ;Normalized master flat filename

  ;Assume that all remaining frames are science data
  spectra=file_list_hb[where(file_type eq 'OBJECT')]
  spectra_name=file_basename(spectra,'.fits')
  
  ;bias section
  control_bias_combine,biaslist,b_type,mbias
   mbias_file=inter_path+'mbias.fits' 
  writefits,mbias_file,mbias.im,mbias.hdr
 
  ;dark section
  control_dark_combine,darklist,d_type,mdark,mbias
  mdark_file=inter_path+'mdark.fits'
  writefits,mdark_file,mdark.im,mdark.hdr
  writefits,mdark_file,mdark.error,mdark.hdr,/APPEND
  sigma_d=mdark.error
  ;flat section
  control_flat_combine,flat_list,f_type,mflat,mbias,mdark
  mflat_file=inter_path+'mflat.fits'
  writefits,mflat_file,mflat.im,mflat.hdr
  writefits,mflat_file,mflat.error,mflat.hdr,/APPEND
  nflat=mflat.im
  sigma_fbdn=mflat.error  
  ;bias correction carried out earlier as la cosmic requires bias corrected images
  for i=0,n_elements(spectra)-1 do begin
    raw_im=mrdfits(spectra[i],1,hdr)
    raw_imb=raw_im-mbias.im
    sigma_imb=sqrt(raw_im+r^2)
    ;add error as well
    sxaddpar, hdr, 'BCFLG', 1 ;bias correction flag
    writefits,inter_path+spectra_name+'_b.fits',raw_imb,hdr
    writefits,inter_path+spectra_name+'_b.fits',sigma_imb,hdr,/APPEND
  endfor
  
  rawb_list= inter_path+spectra_name+'_b.fits'
  rawbc_list= inter_path+spectra_name+'_bc.fits' 
  ;cosmic ray correction
  la_cosmic,rawb_list,outlist=rawbc_list; add other parameters after testing
  ;dark and flat correction
  for i=0,n_elements(spectra)-1 do begin
    raw_imb=mrdfits(rawbc_list[i],1,hdr)
    sigma_imb=mrdfits(rawbc_list[i],2,hdr) ; error after bias correction
    raw_imbd=raw_imb-mdark.im
    sigma_imbd=sqrt((sigma_imb)^2 + (sigma_d)^2) ;error after dark correction
    sxaddpar, hdr, 'DCFLG', 1 ;dark correction flag
    if save_temp_files eq 1 then writefits,inter_path+spectra_name+'_bcd.fits',raw_imbd,hdr
    if save_temp_files eq 1 then writefits,inter_path+spectra_name+'_bcd.fits',sigma_imbd,hdr,/APPEND
    raw_imbdf=raw_imbd/nflat
    sigma_imbdf=sqrt((sigma_imbd/nflat)^2 + ((raw_imbd/nflat^2)* sigma_fbdn)^2) ;error
    sxaddpar, hdr, 'FCFLG', 1 ;dark correction flag
    if save_temp_files eq 1 then writefits,inter_path+spectra_name+'_bcdf.fits',raw_imbdf,hdr
    if save_temp_files eq 1 then writefits,inter_path+spectra_name+'_bcdf.fits',sigma_imbdf,hdr,/APPEND
    
    
    ;spectrum identification and extraction
  
    ;background correction steps
  
    ;extraction steps
  
    ;wavelength calibration
  
  ;flux calibration
 endfor 
  
    
    
    
    
  

  


end