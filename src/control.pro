; NAME:
;      CONTROL
;      
; PURPOSE:
;      
; CALLING SEQUENCE: 
;     CONTROL,parameter file
;      
;
; INPUTS: 
;   Parameter file: Parameter file provided with the distribution providing details defining the pipeline operation.
;      
; OUTPUT:
;   Produce science quality dataproducts for CUTE mission (1D flux and wavelength calibrated spectra, 2D calibrated spectra, Transit light curves).   
;     
; PROCEDURE:
;      
; MODIFICATION HISTORY:
;      created 23.12.2018 by A. G. Sreejith
;data quality implemented till trace

pro control,file,help=help
!quiet = 1
close,/all
if keyword_set(help) then begin
  print,'                               CUTE AUTONOMOUS DATA REDUCTION PIPELINE'
  print,''
  print,'************************************************************************************************************************'
  print, '  Options avaliable in configuration file'
  print,' This software is intented to be fully automated, aimed at producing science-quality output with a single command line with zero user interference for CUTE data.'
  print,' It can be easily used for any single order spectral data in any wavelength without any modification.'
  print,''
  print,'01. Path locations: Data path, intermediate file path and output file path.'
  print,'    If data path is not provided it is assumed to be the current directory.'
  print,'    If intermediate file path and output file path are not provided they are created in the data path.'
  print,'02. The user has the option to select individual steps during reduction. The options are as follows'
  print,'    hb_correction,create_mbias,create_mdark,create_mflat,cr_bias,cr_dark,cr_flat,cr_cosmic,extract,bg_sub,wcalib,fluxcalib,light_curve,retrieval,level3,all'
  print,'    hb_correction: Hot and bad pixel correction'
  print,'    create_mbias: Create master bias.'
  print,'    create_mdark: Create master dark.'
  print,'    create_mflat: Create master flat.'
  print,'    cr_bias: Correct for bias in science frames.'
  print,'    cr_dark: Correct for dark in science frames.'
  print,'    cr_flat: Correct for flat in science frames.'
  print,'    cr_cosmic: Correct for cosmic rays in science frames.'
  print,'    extract: Define the trace and extract the spectrum.'
  print,'    bg_sub: Subtract background from the spectrum.'
  print,'    wcalib: Do wavelength calibration.'
  print,'    fluxcalib: Do Flux calibration.'
  print,'    light_curve: Create three light curves form the data.'
  print,'    retrieval: Process the flux and wavelength calibrated spectra to obtain transmission spectra.'
  print,'    level3: Carry out all processes required for generating Level 3 data of CUTE, ie., flux and wavelength calibrated 1-D spectra.'
  print,'    all: Execute all the steps. This is also the default option if steps keyword is absent.'
  print,'03. Location of hot and bad pixel map. The location of the map of cosmetic defects on the CCD as a FITS file with the same size as the input CCD frame'
  print,'   byte-type pixels (8-bit) and values of 0B for good pixels and 1B for bad ones.'
  print,'04. Location of master calibration files if available.'
  print,'   This includes locations of master bias, master dark and/or master flat files (as fits files).'
  print,'05. Set save_temp_files equal to 1 if you need all the intermediate files to be saved to the intermediate file directory. Set it to 0 if you do not need them. Default is 1.'
  print,'06. Options for master file creation.'
  print,'   Statistical method to create master files. Options are: median, mean or mode'
  print,'   Threshold for rejection of pixels: sigma_deviation (how many sigmas from mean value of the frame).'
  print,'   Saturation_level sets the saturation limit of CCD. If not set default of 72000 is used.'
  print,'07. Cosmic Ray correction options: Level of cr clipping for LA cosmic algorithm. Default value set in the code is 8.'
  print,'08. Trace parameters'
  print,'   Degree for centroid polynomial: Default is 1.'
  print,'   Trace type used for extraction. The options are: simple,fixed,variable,function. For details of these options refer to the manual.'
  print,'   Additional parameters required for different trace options:'
  print,'     centroid  : Centroid will be calculated if not provided'
  print,'     slope     : Required for option: simple and fixed, which defines the slope of the spectrum'
  print,'     width     : Required for option: simple, which defines the width of the spectrum'
  print,'     upper     : Required for option: fixed, which defines the upper width of the spectrum'
  print,'     lower     : Required for option: fixed, which defines the lower width of the spectrum'
  print,'     threshold : Required for option: variable or function, which defines the threshold from the maximum/peak of the spectrum.'
  print,'   File location for information of where to extract background. The user could also provide it as fixed number which corresponds to the shift in pixels from centroid.'
  print,'09. Wavelength calibration variables'
  print,'   Type of wavelength calibration required: Options are simple and crscor'
  print,'   Location of wavelength file. A text file with length equal to number of pixels in cross dispersion direction representing wavelength to pixel mapping.'
  print,'   A location of look up table for stellar parameters. Required only for the option crscor so as to compare with a model data.'
  print,'   If this look up table is not provided then the model file should be a two column model file of wavelength vs flux.'
  print,'   Location of synthetic spectrum (model_file) for cross-correlation. CONTROL assumes that folders are named by their temperature and files are named as model.flx or'
  print,'   a two column model file of wavelength vs flux, if stellar parameters look up table is not provided.'
  print,'10. Flux calibration variables'
  print,'   Location of flux calibration file which provides wavelength vs response relation.'
  print,'  To remove problematic files from processing, create the file badfiles.txt in the data folder'
  print,'*************************************************************************************************************************************************************************'
  
  print,'Operation steps of pipeline are as follows'
  print,' The program works in a series of steps following standard CCD reduction techniques the user can also select individual modules using the step function described above.'
  print,' A reduction log is created to help the users with processes carried out and mentioning the different parameters used if any.'
  print,' It also creates an errorlog, listen all errors that occured during the process.'
  print,' The steps involved in the program execution are as follows:'
  print,''
  print,'   1.  Prepare: Check for the different input parameters are set variables accordingly. If data files are present classify them according to file type.'
  print,'   2.  Hot and bad pixel correction: Correct for hot and bad pixels in the frames.'
  print,'   3.  Create master bias: Create a master bias if master bias file is not found from a set of bias files.'
  print,'   4.  Create master dark: Create a master dark if master dark file is not found from a set of dark files.'
  print,'   5.  Create master flat: Create a master flat if master flat file is not found from a set of flat files.'
  print,'   6.  Correct for bias: Correction of the effect of bias in science frames.'
  print,'   7.  Correct for dark: Correction of the effect of dark in science frames.'
  print,'   8.  Correct for flat: Correction of the effect of flat in science frames.'
  print,'   9.  Correct for cosmic rays: Correction for the effect of cosmic rays using LA cosmic algorithm.'
  print,'   10. Extract spectrum: Define spectrum trace and extract the spectrum.'
  print,'   11. Subtract Background: Subtract background from spectrum.'
  print,'   12. Wavelength calibration: Do wavelength calibration'
  print,'   13. Flux calibration: Do flux calibration'
  print,'   14. Default Light curve: Create 3 light curves (short, middle and long wavelength)'
  print,'   15. Transmission spectrum: Retrieve transmission spectrum.'
  print,'***************************************************************************************************************************************************************************'
  print, '  For additional help and more options, go to the following website: https://github.com/agsreejith/CONTROL'
  return
endif

if keyword_defined(file) eq 0 then begin
  CD, Current=cur_dir
  file=cur_dir+'control_parameters.txt'
  print,'Configuration file assumed to be in current directory'
endif
if file_test(file) eq 0 then begin
  print,'CONTROL: Missing configuration file. Please re-run CONTROL with a configuration file'
  print,'CONTROL: Users should have received a configuration file with this software distribution if not please refer to the configuration file format below.'
  print,'CONTROL: The configuration file should have the following parameters'
  print,'#########################################################################################################################'
  print,'#'
  print,'#         CONFIGURATION FILE FOR CONTROL v1.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE'
  print,'#'
  print,'#************************************************************************************************************************'
  print,'print,#Please note: Comments and empty lines have to start with #.'
  print,'#Use # also for parameters not used.'
  print,'#'
  print,'#************************************************************************************************************************'
  print,'#'
  print,'data_path       = '
  print,'#temp_path      ='
  print,'#out_path       ='
  print,'steps           = '
  print,'hb_map          = '
  print,'#master_bias_file   ='
  print,'#master_dark_file   ='
  print,'#master_flat_file   ='
  print,'save_temp_files     = '
  print,'bias_combine_type   = '
  print,'dark_combine_type   = '
  print,'flat_combine_type   = '
  print,'saturation_limit  = '
  print,'sigma_deviation     = '
  print,'#cosmic_ray_clip    ='
  print,'#cent_poly_degree   ='
  print,'trace_type      = '
  print,'centroid      = '
  print,'slope         = '
  print,'width         = '
  print,'upper         = '
  print,'background_trace  = '
  print,'wavecal_mode    = '
  print,'wavecal_file    = '
  print,'stellar_params    = '
  print,'model_file      = '
  print,'flux_cal      = '
  print,'#************************************************************************************************************************'
  print,'Run CONTROL with the option help to get the details about configuration file keywords or please refer to the software manual.'
  message,'CONTROL: Exiting as configuration file is not avaliable'
endif
infile=gm_read_textstructure(file)
dpflg=0
t=systime(/UTC)
tsec=systime(/seconds)
if (tag_exist(infile,'steps') eq 0) then in_step='all' else in_step=infile.steps
steps=strsplit(in_step,',',/EXTRACT)
hb=1
crmb=1
crmd=1
crmf=1
crb=1
crd=1
crf=1
cr_cor=1
extr=1
bgs=1
wcl=1
fcl=1
sci_out=1
cr_s=1

if where(steps eq 'all') eq -1 then begin
   ;hb_correction,create_mbias,create_mdark,create_mflat,cr_bias,cr_dark,cr_flat,cr_cosmic,extract,wcalib,fluxcalib,light_curve,cr_systematics
  if where(steps eq 'hb_correction') eq -1 then hb=0
  if where(steps eq 'create_mbias') eq -1 then crmb=0  
  if where(steps eq 'create_mdark') eq -1 then crmd=0
  if where(steps eq 'create_mflat') eq -1 then crmf=0
  if where(steps eq 'cr_bias') eq -1 then crb=0
  if where(steps eq 'cr_dark') eq -1 then crd=0
  if where(steps eq 'cr_flat') eq -1 then crf=0
  if where(steps eq 'cr_cosmic') eq -1 then cr_cor=0 else begin
    cr_cor=1
    crb=1
  endelse 
  if where(steps eq 'extract') eq -1 then extr=0
  if where(steps eq 'bg_sub') eq -1 then  bgs=0
  if where(steps eq 'wcalib') eq -1  then wcl=0
  if where(steps eq 'fluxcalib') eq -1 then fcl=0
  if where(steps eq 'light_curve') eq -1 then sci_out=0
  if where(steps eq 'retrieval') eq -1 then cr_s=0
  if where(steps eq 'level3') ne -1 then cr_s=0
endif

idl_ver=float(!Version.RELEASE)
;reading parameters
if(tag_exist(infile,'data_path') eq 1) then data_path=infile.data_path
if(tag_exist(infile,'temp_path') eq 1) then inter_path=infile.temp_path
if(tag_exist(infile,'out_path') eq 1) then out_path=infile.out_path
if(tag_exist(infile,'hbpix_corection_type') eq 1) then hb_type=infile.hbpix_corection_type else hb_type='interpolate'  
if(tag_exist(infile,'bias_combine_type') eq 1) then b_type=infile.bias_combine_type
if(tag_exist(infile,'dark_combine_type') eq 1) then d_type=infile.dark_combine_type
if(tag_exist(infile,'flat_combine_type') eq 1) then f_type=infile.flat_combine_type
if(tag_exist(infile,'cosmic_ray_clip') eq 1) then crclip=infile.cosmic_ray_clip else crclip = 8
if(tag_exist(infile,'cosmic_before_dark') eq 1) then cscr_bfd =fix(infile.cosmic_before_dark) else cscr_bfd = 0  
if(tag_exist(infile,'save_temp_files') eq 1) then save_temp_files=infile.save_temp_files else save_temp_files = 1
;if(tag_exist(infile,'wave_calibration') eq 1) then wcl=infile.wave_calibration else wcl=1
;if(tag_exist(infile,'flux_calibration') eq 1) then fcl=infile.flux_calibration else fcl=1
if(tag_exist(infile,'hb_map') eq 1) then hbmask=infile.hb_map 
if(tag_exist(infile,'saturation_limit') eq 1) then sat_value=infile.saturation_limit
if(tag_exist(infile,'sigma_deviation') eq 1) then threshold=infile.sigma_deviation
science= sci_out
;if(tag_exist(infile,'cr_correction') eq 1) then cr_cor=infile.cr_correction else cr_cor=1
if n_elements(data_path) eq 0 then begin
  dpflg=1
  CD, Current=data_path
  CASE StrUpCase(!Version.OS_Family) OF
    'WINDOWS': data_path=data_path+'\' ;WINDOWS
    'UNIX': data_path=data_path+'/'; UNIX.
  ENDCASE
endif else begin 
  data_path=detectos(data_path)
  if (file_test(data_path,/DIRECTORY)eq 0) then begin
    logprint,'CONTROL: No data file path found. Creating data file path based on configuration file'
    FILE_MKDIR,data_path
  endif
endelse
if n_elements(out_path) eq 0 then begin
  CASE StrUpCase(!Version.OS_Family) OF
    'WINDOWS': out_path=data_path+'output\' ;WINDOWS
    'UNIX': out_path=data_path+'output/'; UNIX.
  ENDCASE
  logprint,'CONTROL: No output file path found. Output file directory created in data directory'
endif else out_path=detectos(out_path)

if (file_test(out_path,/DIRECTORY)eq 0) then begin
  errorlog,'CONTROL: No output file path found. Creating output file path based on configuration file'
  FILE_MKDIR,out_path
endif

logprint,' _______  _______  __    _  _______  ______    _______  ___',logfile=out_path+'control_log'+string(tsec,Format='(I12)')+'.txt'     
logprint,'|       ||       ||  |  | ||       ||    _ |  |       ||   |'    
logprint,'|       ||   _   ||   |_| ||_     _||   | ||  |   _   ||   |'    
logprint,'|       ||  | |  ||       |  |   |  |   |_||_ |  | |  ||   |'    
logprint,'|      _||  |_|  ||  _    |  |   |  |    __  ||  |_|  ||   |___' 
logprint,'|     |_ |       || | |   |  |   |  |   |  | ||       ||       |'
logprint,'|_______||_______||_|  |__|  |___|  |___|  |_||_______||_______|'
logprint,'---------------------------------------------------------------'
logprint,'    CONTROL v1.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE      '
logprint,'---------------------------------------------------------------'
logprint,' 
logprint,'Control execution log dated:',logonly=1
logprint,t,logonly=1  

errorlog,' _______  _______  __    _  _______  ______    _______  ___',logfile=out_path+'control_error_log'+string(tsec,Format='(I12)')+'.txt',logonly=1
errorlog,'|       ||       ||  |  | ||       ||    _ |  |       ||   |',logonly=1
errorlog,'|       ||   _   ||   |_| ||_     _||   | ||  |   _   ||   |',logonly=1
errorlog,'|       ||  | |  ||       |  |   |  |   |_||_ |  | |  ||   |',logonly=1
errorlog,'|      _||  |_|  ||  _    |  |   |  |    __  ||  |_|  ||   |___',logonly=1
errorlog,'|     |_ |       || | |   |  |   |  |   |  | ||       ||       |',logonly=1
errorlog,'|_______||_______||_|  |__|  |___|  |___|  |_||_______||_______|',logonly=1
errorlog,'---------------------------------------------------------------',logonly=1
errorlog,'    CONTROL v1.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE      ',logonly=1
errorlog,'---------------------------------------------------------------',logonly=1
errorlog,'',logonly=1
errorlog,'Control error log dated:',logonly=1
errorlog,t,logonly=1
                                                                                                                 
if dpflg eq 1 then errorlog,'CONTROL: No data file path found. Data files assumed to be in current directory''
  ;path selection
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
    errorlog,'CONTROL: No temporary file path found. Temporary file directory created in data directory'
  endif else inter_path=detectos(inter_path)
  
  
  if n_elements(out_path) eq 0 then begin
        CASE StrUpCase(!Version.OS_Family) OF
          'WINDOWS': out_path=data_path+'output\' ;WINDOWS
          'UNIX': out_path=data_path+'output/'; UNIX. 
        ENDCASE 
    errorlog,'CONTROL: No output file path found. Output file directory created in data directory'
  endif else out_path=detectos(out_path)
  
  ;creates paths if they dont exist
  if (file_test(data_path,/DIRECTORY)eq 0) then begin
    errorlog,'CONTROL: No data file path found. Creating data file path based on configuration file'
    FILE_MKDIR,data_path
  endif  
  if (file_test(inter_path,/DIRECTORY)eq 0) then begin 
    errorlog,'CONTROL: No temporary file path found. Creating temporary file path based on configuration file'
    FILE_MKDIR,inter_path
  endif
  if (file_test(out_path,/DIRECTORY)eq 0) then begin
    errorlog,'CONTROL: No output file path found. Creating output file path based on configuration file'
    FILE_MKDIR,out_path
  endif
  if save_temp_files lt 0 or save_temp_files gt 1 then begin
    logprint,'Allowed values for save temporary files are 1 or 0. For any other values default of 1 will be assumed.'
    save_temp_files=1
  endif  
  
  all_file_list=file_search(data_path+'*.fits') ;Get file list
  ;unique files?
  if (n_elements(all_file_list) eq 0) then begin
    errorlog,'CONTROL: File acess error, No fits file found!!',logonly = logonly
    errorlog,'CONTROL: Verify data path.'
    message,'File acess error, No fits file found for reduction!!'
    return
  endif
  badfilefile =data_path+'badfiles.txt'  
  if file_test(badfilefile) then readcol, badfilefile, badfiles, format='a' else badfiles=string(0)
  all_file_names=file_basename(all_file_list,'.fits')
   all_file_names= all_file_names+'.fits'
  match2,all_file_names,badfiles,suba,subb
  file_list_loc=where(suba eq -1)
  file_list=all_file_list[file_list_loc]
  file_names=file_basename(file_list,'.fits') ;Get file names
  file_type=strarr(n_elements(file_list))
  if n_elements(file_list) eq 0 then begin
    errorlog,'CONTROL: File acess error, No good fits file found!!',logonly = 1
    errorlog,'CONTROL: Verify data path.'
    message,'File acess error, No good fits file found for reduction!!'
  endif 

  ;carry out bad and hot pixel corrections
  if(hb) then begin
    if n_elements(hbmask) eq 0 then begin
      errorlog,'CONTROL: No mask file found for hot and bad pixel correction',logonly = 1
      message,'CONTROL: No mask file found for hot and bad pixel correction'
    endif  
    hbmask=detectos(hbmask) 
    hbmask=string(hbmask)
    if file_test(hbmask) then begin
      mask=mrdfits(hbmask,0,mask_hdr,/SILENT)
      logprint,'Applying hot and bad pixel correction based on input mask file.'
      logprint,'Hot and bad pixel corrected using '+hb_type+' method'
      for i=0,n_elements(file_list)-1 do begin
        control_hotbad,file_list[i],mask,out,hb_type
        hdr=out.header
        out_im=out.data
        sxaddpar, hdr, 'HBFLG', 1. ;Hot and bad pixel correction flag
        writefits,inter_path+file_names[i]+'_hb.fits',out_im,hdr ;hotbad file 
        file_type[i]=SXPAR(hdr, 'FILETYPE',MISSING=-1);Find object description
      endfor
      file_list_hb=file_search(inter_path+'*_hb.fits')
    endif else begin
      logprint,'No input mask file found for hot and bad pixel correction, proceeding without hot and bad pixel correction.'
      for i=0,n_elements(file_list)-1 do begin
        hdr=HEADFITS(file_list[i],exten=0) 
        file_type[i]=SXPAR(hdr, 'FILETYPE',MISSING=-1)
        file_list_hb=file_search(data_path+'*.fits')
      endfor
    endelse    
  endif else begin
    logprint,'Proceeding without hot and bad pixel correction as requested by the user.'
    for i=0,n_elements(file_list)-1 do begin
      hdr=HEADFITS(file_list[i],exten=0) 
      file_type[i]=SXPAR(hdr, 'FILETYPE',MISSING=-1)
      file_list_hb=file_search(data_path+'*.fits')
    endfor  
  endelse
  
  
  ;Automatic file classification

  ;Compile a list of BIAS frames
  biaslist_ind=where(strpos(file_type,'BIAS') ge 0)
  if biaslist_ind[0] ne -1 then biaslist=file_list_hb[biaslist_ind] else begin
    errorlog,'CONTROL: No BIAS files found.'
    ;message,'CONTROL: No BIAS files found. CONTROL will exit now.'
  endelse  
  ;bias_file = prefix+inst_mode+'.bias.fits'         ;Master bias filename

  ;Compile a list of DARK frames
  darklist_ind=where(strpos(file_type,'DARK') ge 0)
  if darklist_ind[0] ne -1 then darklist=file_list_hb[darklist_ind] else begin
    errorlog,'CONTROL: No DARK files found.'
    ;message,'CONTROL: No DARK files found. CONTROL will exit now.'
  endelse
  ;dark_file = prefix+inst_mode+'.dark.fits'         ;Master bias filename

  ;Compile a list of FLAT frames
  flatlist_ind=where(strpos(file_type,'FLAT') ge 0)
  if flatlist_ind[0] ne -1 then flatlist=file_list_hb[flatlist_ind] else begin
    errorlog,'CONTROL: No FLAT files found.'
    ;message,'CONTROL: No DARK files found. CONTROL will exit now.'
  endelse
  ;flat_file=prefix+inst_mode+'.flat.fits'           ;Master flat filename
  ;norm_flat_file=prefix+inst_mode+'.flat.norm.fits' ;Normalized master flat filename

  ;Assume that all remaining frames are science data
  spectra_ind=where(strpos(file_type, 'OBJECT') ge 0)
  if spectra_ind[0] ne -1 then begin
    spectra=file_list_hb[spectra_ind] 
    spectra_name=file_basename(spectra,'.fits')
  endif else begin
    errorlog,'CONTROL: No SCIENCE files found.'
    ;message,'CONTROL: No DARK files found. CONTROL will exit now.'
  endelse
  ;calculate readout noise and gain 
  if n_elements(biaslist) eq 0 then logprint,'CONTROL: Cannot calculate readnoise and gain without bias files. These values will be read from headers of master files' else begin
     b1=mrdfits(biaslist[0],0,hdr,/SILENT)
     b2=mrdfits(biaslist[1],0,hdr,/SILENT)   
     if n_elements(flatlist) eq 0 then begin
      logprint,'CONTROL: Cannot calculate gain without flat files. These values will be read from headers of master files'
      if tag_exist(infile,'master_flat_file') eq 0 then logprint,'CONTROL: No master flat file found to read gain values reading them from bias files'
      g_calc=SXPAR( hdr, 'CCDGAIN')
     endif else begin
       f1=mrdfits(flatlist[0],0,hdr,/SILENT)
       f2=mrdfits(flatlist[1],0,hdr,/SILENT) 
       g_calc=gain(f1,f2,b1,b2)
       g_calc=1
     endelse  
     R=rnoise(b1,b2,g_calc)
  endelse      
  for i=0,n_elements(biaslist)-1 do begin
    h = headfits(biaslist[i])
    sxaddpar,h,'RNOISE',R
    sxaddpar,h,'CCDGAIN',g_calc
    modfits,biaslist[i],0,h
  endfor

  ;bias section
  if (crmb) then begin
    logprint,'CONTROL will create master bias now'
    if tag_exist(infile,'master_bias_file') eq 0 then begin
      if n_elements(biaslist) eq 0 then begin
        errorlog,'CONTROL: File access error, No master bias or bias file found',logonly = logonly
        logprint,'CONTROL: Verify data path or check bias file header. CONTROL checks for FILETYPE keyword in header to clasify files'
        message,'CONTROL: File access error, No master bias or bias file found!'
      endif
      control_bias_combine,biaslist,mbias,type=b_type,sat_value=sat_value,threshold=threshold
      mbias_file=inter_path+'mbias.fits' 
      mwrfits,mbias.im,mbias_file,mbias.hdr,/create, /SILENT
      mwrfits,mbias.dq,mbias_file,mbias.hdr
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
      mb_dq=mbias.dq
    endif else begin
      mbfile_path=detectos(infile.master_bias_file)
      mbfile=file_search(mbfile_path)
      if (file_test(mbfile_path) eq 0) then begin 
        errorlog,'CONTROL: File access error, No master bias file found!!',logonly=1
        Message,'CONTROL: File access error, No master bias file found!!' 
      endif else  begin 
        mbias_data=mrdfits(mbfile,0,mbhdr,/SILENT)
        r=float(SXPAR( mbhdr, 'RNOISE'))
        mbnxy=size(mbias)
        mbnx=mbnxy[1]
        mbny=mbnxy[2]
        FITS_INFO, mbfile,N_ext =numext
        if numext eq 1 then mb_dq=mrdfits(mbfile,1,hdr,/SILENT) else mb_dq=bytarr(mbnx,mbny)
        mbias={im:mbias_data,hdr:mbhdr,dq:mb_dq}
      endelse  
   endelse
  endif else begin
    if tag_exist(infile,'master_bias_file') ne 0 then begin
      mbfile_path=detectos(infile.master_bias_file)
      mbfile=file_search(mbfile_path)
      if (file_test(mbfile_path) eq 0) then begin
        errorlog,'CONTROL: File access error, No master bias file found!!',logonly=1
        Message,'CONTROL: File access error, No master bias file found!!'
      endif else begin 
        mbias_data=mrdfits(mbfile,0,mbhdr,/SILENT)
        r=float(SXPAR( mbhdr, 'RNOISE'))
        mbnxy=size(mbias)
        mbnx=mbnxy[1]
        mbny=mbnxy[2]
        FITS_INFO, mbfile,N_ext =numext
        if numext gt 0 then mb_dq=mrdfits(mbfile,1,hdr,/SILENT) else mb_dq=bytarr(mbnx,mbny)
        mbias={im:mbias_data,hdr:mbhdr,dq:mb_dq}
      endelse  
    endif  
  endelse  

  ;dark section
  if (crmd) then begin
    logprint,'CONTROL will create master dark now'
    if tag_exist(infile,'master_dark_file') eq 0 then begin
      if n_elements(darklist) eq 0 then begin
        errorlog,'CONTROL: File access error, No master dark or dark file found',logonly = 1
        logprint,'CONTROL: Verify data path or check dark file header. CONTROL checks for FILETYPE keyword in header to clasify files.'
        message,'CONTROL: File access error, No master dark or dark file found!'
      endif
        control_dark_combine,darklist,mdark,mbias,type=d_type,sat_value=sat_value,threshold=threshold
        mdark_file=inter_path+'mdark.fits'
        mwrfits,mdark.im,mdark_file,mdark.hdr,/create
        mwrfits,mdark.error,mdark_file,mdark.hdr, /SILENT
        mwrfits,mdark.dq,mdark_file,mdark.hdr, /SILENT
        sigma_d=mdark.error
        md_dq=mdark.dq
        dhdr=mdark.hdr
    endif else begin
      mdfile_path=detectos(infile.master_dark_file)
      mdfile=file_search(mdfile_path)
      if (file_test(mdfile_path) eq 0) then begin 
        errorlog,'CONTROL: File access error, No master dark file found!!',logonly=1
        message,'CONTROL: File access error, No master dark file found!!'
        return 
      endif else begin
        fits_info, mdfile, SILENT=silent,N_ext=extension
        mdark_data=mrdfits(mdfile,0,dhdr,/SILENT)
        nxyd=size(mdark_data)
        if (extension gt 0) then mdark_err=mrdfits(mdfile,1,dhdr2,/SILENT) else begin
          mdark_err=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master dark file does not contain error, assuming error to be zero'
        endelse
        if extension gt 1 then md_dq=mrdfits(mdfile,2,dhdr3,/SILENT) else md_dq=bytarr(nxyd[1],nxyd[2])  
      endelse 
      sigma_d= mdark_err
      mdark={im:mdark_data,error:mdark_err,hdr:dhdr,dq:md_dq} 
    
    endelse
  endif else begin
    if tag_exist(infile,'master_dark_file') ne 0 then begin
      mdfile_path=detectos(infile.master_dark_file)
      mdfile=file_search(mdfile_path)
      if (file_test(mdfile_path) eq 0) then begin
        errorlog,'CONTROL: File access error, No master dark file found!!',logonly=1
        message,'CONTROL: File access error, No master dark file found!!'
        return
      endif else begin
        fits_info, mdfile, SILENT=silent,N_ext=extension
        mdark_data=mrdfits(mdfile,0,dhdr,/SILENT)
        nxyd=size(mdark_data)
        if (extension gt 0) then mdark_err=mrdfits(mdfile,1,dhdr2,/SILENT) else begin
          mdark_err=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master dark file does not contain error, assuming error to be zero'
        endelse
        if extension gt 1 then md_dq=mrdfits(mdfile,2,dhdr3,/SILENT) else md_dq=bytarr(nxyd[1],nxyd[2]) 
      endelse
      sigma_d= mdark_err
      mdark={im:mdark_data,error:mdark_err,hdr:dhdr,dq:md_dq}
    endif
  endelse  

  ;flat section
  if (crmf) then begin
    logprint,'CONTROL will create master flat now'
    if tag_exist(infile,'master_flat_file') eq 0 then begin
      if n_elements(flatlist) eq 0 then begin
        errorlog,'CONTROL: File access error, No master flat or flat file found',logonly = 1
        logprint,'CONTROL: Verify data path or check flat file header. CONTROL checks for FILETYPE keyword in header to clasify files.'
        message,'CONTROL: File access error, No master flat or flat file found!'
      endif
      control_flat_combine,flatlist,mflat,mbias,mdark,type=f_type,sat_value=sat_value,threshold=threshold
      mflat_file=inter_path+'mflat.fits'
      mwrfits,mflat.im,mflat_file,mflat.hdr,/create
      mwrfits,mflat.error,mflat_file,mflat.hdr
      mwrfits,mflat.dq,mflat_file,mflat.hdr
      nflat=mflat.im
      sigma_fbdn=mflat.error
      mf_dq=mflat.dq  
    endif else begin
      mffile_path=detectos(infile.master_flat_file)
      mffile=file_search(mffile_path)
      if (file_test(mffile_path) eq 0) then begin 
        errorlog,'CONTROL: File access error, No master flat file found!!',logonly=1
        message,'CONTROL: File access error, No master flat file found!!'
        return 
      endif  else begin
        fits_info, mffile, SILENT=silent,N_ext=extension
        nflat=mrdfits(mffile,0,mfhdr,/SILENT)
        nxyd=size(nflat)
        if (extension gt 0) then sigma_fbdn=mrdfits(mffile,1,mdhdr,/SILENT) else begin
          sigma_fbdn=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master flat file does not contain error, assuming error to be zero'
        endelse
        if extension gt 1 then mf_dq=mrdfits(mdfile,2,hdr,/SILENT) else mf_dq=bytarr(nxyd[1],nxyd[2]) 
      endelse
    endelse
  endif else begin
    if tag_exist(infile,'master_flat_file') ne 0 then begin
      mffile_path=detectos(infile.master_flat_file)
      mffile=file_search(mffile_path)
      if (file_test(mffile_path) eq 0) then begin
        logprint,'CONTROL: Verify data path or check bias file header. CONTROL checks for FILETYPE keyword in header to clasify files.'
        errorlog,'CONTROL: File access error, No master flat file found!!',logonly=1
        message,'CONTROL: File access error, No master flat file found!!'
        return
      endif  else begin
        fits_info, mffile, SILENT=silent,N_ext=extension
        nflat=mrdfits(mffile,0,mfhdr,/SILENT)
        nxyd=size(mflat_data)
        if (extension gt 0) then sigma_fbdn=mrdfits(mffile,1,mdhdr,/SILENT) else begin
          sigma_fbdn=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master flat file does not contain error, assuming error to be zero'
        endelse
        if extension gt 1 then mf_dq=mrdfits(mdfile,2,hdr,/SILENT) else mf_dq=bytarr(nxyd[1],nxyd[2])
      endelse
    endif
  endelse 
  ;bias correction
  if datatype(mbias,2) ne 8 then begin
    logprint,'CONTROL could not find master bias file. Do you want to continue reducing this spectrum assuming bias to be zero'
    logprint,'Press q to exit. Press any key to use bias as zero.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master bias of zeros of size 2048x515 and assuming readout noise to be 12.25'
    sxaddpar,mbhdr,'RNOISE',12.25,'Read out noise'
    mbias={im:dblarr(2048,515),error:dblarr(2048,515),hdr:mbhdr,dq:bytarr(2048,515)}
    mb_dq=bytarr(2048,515)
    
  endif
  if (n_elements(spectra) eq 0) then begin
    logprint,'CONTROL: No science files found. CONTROL will exit now.',logonly = 1
    errorlog,'CONTROL: No science files found. CONTROL will exit now.',logonly = 1
    message,'CONTROL: No science files found. CONTROL will exit now.'
  endif
  ;bias correction carried out earlier as la cosmic requires bias corrected images
  dummy=mrdfits(spectra[0],0,hdr,/SILENT)
  nx=(size(dummy))[1]
  ny=(size(dummy))[2]
  ycut1=SXPAR( hdr, 'YCUT1')
  ycut2=SXPAR( hdr, 'YCUT2')
  ylen=ycut2-ycut1+1
  if (crb) then begin
    for i=0,n_elements(spectra)-1 do begin
      raw_im=mrdfits(spectra[i],0,hdr,/SILENT)
      nx=(size(raw_im))[1]
      ny=(size(raw_im))[2]
      FITS_INFO, spectra[i],N_ext =numext,/SILENT
      if numext eq 2 then dq_im=mrdfits(spectra[i],1,hdr,/SILENT) else dq_im=bytarr(nx,ny)
      if (SXPAR(hdr, 'BCFLG',MISSING=-1)) ne -1 then begin
        blfg_chk=SXPAR(hdr, 'BCFLG',MISSING=-1)
        if blfg_chk eq 1 then begin
          logprint,'CONTROL: The spectrum('+spectra[i]+') is already bias corrected. Skipping the bias correction for current spectrum.'
          goto,spectrum_bias_cr_end
        endif
      endif
      mnx=(size(mbias.im))[1]
      mny=(size(mbias.im))[2]
      ycut1=SXPAR( hdr, 'YCUT1')
      ycut2=SXPAR( hdr, 'YCUT2')
      ylen=ycut2-ycut1+1
      ccd_gain=SXPAR( hdr, 'GAIN')
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
      if nx eq mnx then begin
        if ny eq mny then begin
          raw_imb=raw_im-mbias.im
          sigma_imb=sqrt(raw_im+r^2)
          ;data quality
          dq_imb=dq_im+mb_dq
        endif else begin
          if (ylen ne ny) then mbias_nw=mbias.im[*,ycut1:ycut1+ny-1] else mbias_nw=mbias.im[*,ycut1:ycut2]
          if (ylen ne ny) then mbdq_nw=mb_dq[*,ycut1:ycut1+ny-1] else mbdq_nw=mb_dq[*,ycut1:ycut2]
          raw_imb=raw_im-mbias_nw
          sigma_imb=sqrt(raw_im+r^2)
          ;data quality
          dq_imb=dq_im+mbdq_nw
          logprint,'CONTROL: Bias correction carried out on spectrum('+spectra[i]+').'
        endelse 
      endif else begin
        errorlog,'CONTROL: X Dimension error with science image and bias image'
        return
      endelse
      ;add error as well
      sxaddpar, hdr, 'BCFLG', 1,'BIAS CORRECTION FLAG' ;bias correction flag
      mwrfits,raw_imb,inter_path+spectra_name[i]+'_b.fits',hdr,/create
      mwrfits,sigma_imb,inter_path+spectra_name[i]+'_b.fits',hdr
      ;data quality
      prb=where(dq_imb ge 1)
      if total(prb) ne -1 then dq_imb[prb]=1
      undefine,prb
      mwrfits,dq_imb,inter_path+spectra_name[i]+'_b.fits',hdr
    endfor
    rawb_list= inter_path+spectra_name+'_b.fits'
  endif else rawb_list = spectra  
  spectrum_bias_cr_end:
  if cscr_bfd eq 1 then begin
    if cr_cor eq 1 then begin
      rawbc_list= inter_path+spectra_name+'_bc.fits'
      mask_list = inter_path+spectra_name+'_bc_mask.fits'
      if (SXPAR(hdr, 'CRFLG',MISSING=-1)) ne -1 then begin
        crfg_chk=SXPAR(hdr, 'CRFLG',MISSING=-1)
        if crfg_chk eq 1 then begin
          logprint,'CONTROL: The spectrum('+spectra[i]+') is already corrected for cosmic rays. Skipping cosmic ray correction for current spectrum.'
          goto,spectrum_cr_cr_end
        endif  
      endif
      ;cosmic ray correction
      if ccd_gain eq 0 then ccd_gain=1
      logprint,'CONTROL: Carrying out cosmic ray correction with clip value of'+strtrim(string(crclip),2)+'.'
      la_cosmic,rawb_list,outlist=rawbc_list,masklist=mask_list,gain=ccd_gain,readn=r,sigclip=crclip; add other parameters after testing
      crflg=1
      for i=0,n_elements(rawbc_list)-1 do begin
        raw_imbc=mrdfits(rawbc_list[i],0,hdr)
        nxc=(size(raw_imbc))[1]
        nyc=(size(raw_imbc))[2]
        raw_imb=mrdfits(rawb_list[i],0,hdr)
        sxaddpar, hdr,'CRFLG',crflg,'COSMIC RAY CORRECTION FLAG'
        modfits,rawbc_list[i],0,hdr,EXTEN_NO=0
        FITS_INFO, rawb_list[i],N_ext=numext,/SILENT
        if numext gt 0 then sigma_imbc=mrdfits(rawb_list[i],1,hdr1,/SILENT) else sigma_imbc=sqrt(raw_imbc+r^2)
        if numext gt 1 then dq_imbc=mrdfits(rawb_list[i],2,hdr2,/SILENT) else  dq_imbc=bytarr(nxc,nyc)
        cc_mask=mrdfits(mask_list[i],0,mask_hdr,/SILENT)
        cr_loc=where(cc_mask eq 1)
        dq_imbc[cr_loc]=1
        prb=where(dq_imbc ge 1)
        if total(prb) ne -1 then dq_imbc[prb]=1
        undefine,prb
        mwrfits,sigma_imbc,rawbc_list[i],hdr
        mwrfits,dq_imbc,rawbc_list[i],hdr
      endfor
      logprint,'CONTROL: Corrections carried out for cosmic rays on spectrum using LA COSMIC, sigma clip used is '+STRTRIM(STRING(crclip),2)+'.'
    endif else begin
      logprint,'CONTROL: Cosmic ray correction not carried out as per user request.'
      rawbc_list=rawb_list
;     for i=0,n_elements(rawb_list)-1 do begin
;      mwrfits,sigma_imb,rawb_list[i],hdr
;      mwrfits,dq_imb,rawb_list[i],hdr
;     endfor
      crflg=0
    endelse
  endif else rawbc_list=rawb_list
  spectrum_cr_cr_end:
  ;dark and flat correction
  if datatype(mdark,2) ne 8 then begin
    logprint,'CONTROL could not find master dark file. Do you want to continue reducing this spectrum assuming dark to be zero'
    logprint,'Press q to exit. Press any key to use dark as zero.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master dark of zeros of size 2048x515' 
    mdark={im:dblarr(2048,515),error:dblarr(2048,515),dq:bytarr(2048,515)}
    sigma_d=dblarr(2048,515)
    md_dq=bytarr(2048,515)
  endif
  if datatype(mflat,2) ne 8 then begin
    logprint,'CONTROL could not find master flat file. Do you want to continue reducing this spectrum assuming flat to be one'
    logprint,'Press q to exit. Press any key to use flat as one.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master flat of ones of size 2048x515'
    sigma_fbdn=make_array(2048,515,value=1.0,/DOUBLE)
    nflat=make_array(2048,515,value=1.0,/DOUBLE)
    mf_dq=bytarr(2048,515)
  endif
  for i=0,n_elements(spectra)-1 do begin
    raw_imbc=mrdfits(rawbc_list[i],0,hdr,/SILENT)
    ;print,rawbc_list[i]
    FITS_INFO, rawbc_list[i],N_ext=numext,/SILENT
    bnx=(size(raw_imbc))[1]
    bny=(size(raw_imbc))[2]

    if numext gt 0 then sigma_imbc=mrdfits(rawbc_list[i],1,hdr1,/SILENT)  else sigma_imbc=dblarr(bnx,bny)    ; error after bias correction
    ;if n_elements(sigma_imbc) eq 1 then sigma_imbc=dblarr(bnx,bny)
    if numext gt 1 then dq_imbc=mrdfits(rawbc_list[i],2,hdr2,/SILENT)  else dq_imbc=bytarr(bnx,bny)
    if n_elements(dq_imbc) eq 1 then dq_imbc=dblarr(bnx,bny)
    if datatype(r,2) eq 0 then r=12.25
    sigma_imbc=sqrt(raw_imbc+r^2)
    if (SXPAR(hdr, 'DCFLG',MISSING=-1)) ne -1 then begin
      dcfg_chk=SXPAR(hdr, 'DCFLG',MISSING=-1)
      if dcfg_chk eq 1 then begin
        logprint,'CONTROL: The spectrum('+rawbc_list[i]+') is already dark corrected. Skipping dark correction for current spectrum.'
        goto,spectrum_dr_cr_end
      endif
    endif
    nx=(size(raw_imb))[1]
    ny=(size(raw_imb))[2]
    ycut1=SXPAR( hdr, 'YCUT1')
    ycut2=SXPAR( hdr, 'YCUT2')
    ylen=ycut2-ycut1+1
    if (crd) then begin
      logprint,'CONTROL will dark correct the science spectrum now'
      dnx=(size(mdark.im))[1]
      dny=(size(mdark.im))[2]
      if nx eq dnx then begin
        if ny eq dny then begin
          raw_imbd=raw_imbc-mdark.im
          sigma_imbd=sqrt((sigma_imbc)^2 + (sigma_d)^2) ;error after dark correction
          dsigma=sigma_d
          ;data quality
          dq_imbd=dq_imbc+md_dq
        endif else begin
          if (ylen ne ny) then begin
            mdark_nw=mdark.im[*,ycut1:ycut1+ny-1]
            sigma_dnw=sigma_d[*,ycut1:ycut1+ny-1]
            mddq_nw=md_dq[*,ycut1:ycut1+ny-1]
          endif else begin
            mdark_nw=mdark.im[*,ycut1:ycut2]
            sigma_dnw=sigma_d[*,ycut1:ycut2]
            mddq_nw=md_dq[*,ycut1:ycut2]
          endelse
          raw_imbd=raw_imbc-mdark_nw
          sigma_imbd=sqrt((sigma_imbc)^2 + (sigma_dnw)^2) ;error after dark correction
          dsigma=sigma_dnw
          ;data quality
          dq_imbd=dq_imbc+mddq_nw
        endelse
      endif else begin
        errorlog,'CONTROL: X Dimension error with science image and dark image'
        return
      endelse
      sxaddpar, hdr, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      sxaddpar, hdr1, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      sxaddpar, hdr2, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      ;hdr_d=update_header(hdr)
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',raw_imbd,hdr
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',sigma_imbd,hdr1, /APPEND
      ;data quality
      prb=where(dq_imbd ge 1)
      if total(prb) ne -1 then dq_imbd[prb]=1
      undefine,prb
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',dq_imbd,hdr2, /APPEND
      logprint,'CONTROL: Dark correction carried out on spectrum('+spectra_name[i]+').'
      rawbd_file=inter_path+spectra_name[i]+'_bcd.fits'
    endif else begin
      rawbd_file=rawbc_list[i]
      raw_imbd=raw_imbc
      sigma_imbd=sigma_imbc
      dq_imbd=dq_imbc
    endelse
    if cr_cor eq 1 then begin
     if (SXPAR(hdr, 'CRFLG',MISSING=-1)) ne 1 then begin
        logprint,'CONTROL: Cosmic Ray correction carried out after dark correction on spectrum('+spectra_name[i]+') with a clip value of '+STRTRIM(STRING(crclip),2)+'.'
        if ccd_gain eq 0 then ccd_gain=1
        rawbdc_file = inter_path+spectra_name[i]+'_bcdc.fits'
        rawbd_mask  = inter_path+spectra_name[i]+'_bdc_mask.fits'
        la_cosmic,rawbd_file,outlist=rawbdc_file,masklist=rawbd_mask,gain=ccd_gain,readn=r,sigclip=crclip; add other parameters after testing
        crflg=1
        sxaddpar, hdr,'CRFLG',crflg,'COSMIC RAY CORRECTION FLAG'
        modfits,rawbdc_file,0,hdr,EXTEN_NO=0
        raw_imbd=mrdfits(rawbdc_file,0,crhdr,/SILENT)
        cc_mask=mrdfits(rawbd_mask,0,mask_hdr,/SILENT)
        cr_loc=where(cc_mask eq 1)
        dq_imbd[cr_loc]=1
        prb=where(dq_imbd ge 1)
        if total(prb) ne -1 then dq_imbd[prb]=1
        undefine,prb
        if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdc.fits',sigma_imbd,hdr1, /APPEND
        if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdc.fits',dq_imbd,hdr2, /APPEND
      endif
    endif  

    spectrum_dr_cr_end:
    if(crf) then begin
      logprint,'CONTROL will flat correct the science spectrum now'
      if (SXPAR(hdr, 'FCFLG',MISSING=-1)) ne -1 then begin
        fcfg_chk=SXPAR(hdr, 'FCFLG',MISSING=-1)
        if fcfg_chk eq 1 then begin
          logprint,'CONTROL: The spectrum('+rawbc_list[i]+') is already flat corrected. Skipping flat correction for current spectrum.'
          goto,spectrum_fl_cr_end
        endif
      endif
      fnx=(size(nflat))[1]
      fny=(size(nflat))[2]
      if nx eq fnx then begin
        if ny eq fny then begin
          raw_imbdf=raw_imbd/nflat
          sigma_imbdf=sqrt((sigma_imbd/nflat)^2 + ((raw_imbd/nflat^2)* sigma_fbdn)^2) ;error
          ;data quality
          dq_imbdf=dq_imbd+mf_dq
        endif else begin
          if (ylen ne ny) then begin
            nflat_nw=nflat[*,ycut1:ycut1+ny-1]
            sigma_fnw=sigma_fbdn[*,ycut1:ycut1+ny-1]
            mfdq_nw=mf_dq[*,ycut1:ycut1+ny-1]
          endif else begin
            nflat_nw=nflat[*,ycut1:ycut2]
            sigma_fnw=sigma_fbdn[*,ycut1:ycut2]
            mfdq_nw=mf_dq[*,ycut1:ycut2]
          endelse
          raw_imbdf=raw_imbd/nflat_nw
          sigma_imbdf=sqrt((sigma_imbd/nflat_nw)^2 + ((raw_imbd/nflat_nw^2)* sigma_fnw)^2) ;error
          ;data quality
          dq_imbdf=dq_imbd+mfdq_nw 
        endelse
      endif else begin
        errorlog,'CONTROL: X Dimension error with science image and flat image'
        return
      endelse
      sxaddpar, hdr, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      sxaddpar, hdr1, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      sxaddpar, hdr2, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',raw_imbdf,hdr
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',sigma_imbdf,hdr1, /APPEND    
      prb=where(dq_imbdf ge 1)
      if total(prb) ne -1 then dq_imbdf[prb]=1
      undefine,prb
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',dq_imbdf,hdr2, /APPEND    
      logprint,'CONTROL: Flat correction carried out on spectrum('+rawbc_list[i]+').' 
      writefits,out_path+spectra_name[i]+'_2d.fits',raw_imbdf,hdr 
      writefits,out_path+spectra_name[i]+'_2d.fits',sigma_imbdf,hdr1, /APPEND 
      writefits,out_path+spectra_name[i]+'_2d.fits',dq_imbdf,hdr2, /APPEND
    endif else begin
      raw_imbdf=raw_imbd
      sigma_imbdf=sigma_imbd
      dq_imbdf=dq_imbd
    endelse
    spectrum_fl_cr_end:
    ;spectrum identification and extraction including background

    if(extr) then begin
      im_bflg=fix(SXPAR(hdr, 'BCFLG',MISSING=-1))
      im_dflg=fix(SXPAR(hdr, 'DCFLG',MISSING=-1))
      im_fflg=fix(SXPAR(hdr, 'FCFLG',MISSING=-1))
      logprint,'CONTROL will extract the science spectrum now'
      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1) then begin
        logprint,'CONTROL: The spectra('+spectra_name[i]+') requested for extraction is not bias or dark or flat corrected.'
        logprint,'Press q to skip the extraction for current spectrum. Press any key to continue extraction for the current spectrum.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the extraction for current spectrum as requested by the user.'
          goto,spectrum_loop_end
;          message,'CONTROL_DARK_COMBINE: Exiting as requested by the user.'
;          err = '%control_dark_combine: Insufficient number of parameters'
;          return
        endif
      endif
      if (fix(SXPAR(hdr, 'EXTFLF',MISSING=-1))) eq 1 then begin
        logprint,'CONTROL: The spectra('+rawbc_list[i]+') requested for extraction is already an extracted spectrum'
        logprint,'CONTROL: Skipping the extraction for current spectrum as requested by the user.'
          goto,extr_read
      endif
      in_image={data:raw_imbdf,error:sigma_imbdf,header:hdr,dq:dq_imbdf}
      if (tag_exist(infile,'trace_type') eq 1 )then t_type = infile.trace_type else begin 
        t_type = 'simple'
        logprint,'CONTROL: Type input for trace not defined using default: simple'
      endelse
      spectrum=control_trace(in_image,infile,t_type,spectra_name[i])
      if (idl_ver ge 8) then begin
        if isa(spectrum,'STRUCT') eq 0 then begin
          errorlog,'CONTROL: CONTROL_TRACE returned without extracted spectrum('+spectra_name[i]+'). Skipping further pipeline procedures for the corresponding file.'
          goto,spectrum_loop_end
        endif  
      endif else begin
        if datatype(spectrum,2) ne 8 then begin
          errorlog,'CONTROL: CONTROL_TRACE returned without extracted spectrum('+spectra_name[i]+'). Skipping further pipeline procedures for the corresponding file.'
          goto,spectrum_loop_end
        endif  
      endelse
      ext_spectrum=spectrum.data
      hdr=spectrum.header
      hdr_nw=update_header(hdr)
     ; sxaddpar, hdr_nw,'TFIELDS',tfileds,'Number of fields in each row'
      sxaddpar, hdr_nw,'TTYPE1  ','COUNTS','label for field   1'
      sxaddpar, hdr_nw,'TFORM1  ','17float','Data format of field 1'
      sxaddpar, hdr_nw,'TUNIT1  ','counts','Physical unit of field 1'
      sxaddpar, hdr_nw,'TDISP1  ','Float','Display format for column  1'
;     sxaddpar, hdr,'TNULL1  ',,'Undefined value for column  1'
;     sxaddpar, hdr,'TTYPE2  ',,'label for field   2'
;     sxaddpar, hdr,'TFORM2  ',,'Data format of field 2'
;     sxaddpar, hdr,'TUNIT2  ',,'Physical unit of field 2'
;     sxaddpar, hdr,'TDISP2  ',,'Display format for column  2'
;     sxaddpar, hdr,'TNULL2  ',,'Undefined value for column  2'
;     sxaddpar, hdr,'INHERIT ',,'Inherit the primary header'
;     sxaddpar, hdr,'EXTNAME ',,'Extension name'
;     sxaddpar, hdr,'EXTVER  ',,'Extension version number'
;     sxaddpar, hdr,'ROOTNAME',,'rootname of the observation set'
      sxaddpar, hdr_nw,'EXTFLF',1,'Extraction flag'
      ;sxdelpar,hdr,'NAXIS2'
      if save_temp_files eq 1 then mwrfits,spectrum,inter_path+spectra_name[i]+'_bcdfe.fits',hdr_nw,/create
    ;  sxaddpar, hdr_nw,'TFIELDS',tfileds,'Number of fields in each row'
;      if save_temp_files eq 1 then mwrfits,spectrum.error,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
;      if save_temp_files eq 1 then mwrfits,spectrum.background,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
;      if save_temp_files eq 1 then mwrfits,spectrum.bck_error,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
;      if save_temp_files eq 1 then mwrfits,spectrum.dq,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
;      if save_temp_files eq 1 then mwrfits,spectrum.dq_bg,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
      ;data:spectrum_val,error:noise,bck:background,bck_error:backg_error,centroid:centroid
      logprint,'CONTROL: Spectrum('+spectra_name[i]+') is extracted.' 
    endif ;else begin
      extr_read:
;      im_bflg=fix(SXPAR(hdr, 'BCFLG',MISSING=-1))
;      im_dflg=fix(SXPAR(hdr, 'DCFLG',MISSING=-1))
;      im_fflg=fix(SXPAR(hdr, 'FCFLG',MISSING=-1))
;      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1) then begin
;        logprint,'CONTROL: The spectra is not bias or dark or flat corrected with CONTROL.'
;        logprint,'Press q to skip the extraction for current spectrum. Press any key to continue assuming the spectrum is corrected.'
;        R = GET_KBRD()
;        if R eq 'q' then begin
;          logprint,'CONTROL: Skipping the extraction for current spectrum as requested by the user.'
;          goto,spectrum_loop_end
;        endif
;      endif
;      spectrum=mrdfits(rawbc_list[i],1,hdr_nw)
;      
;      if(tag_exist(spectrum,'data') eq 0) then begin 
;        logprint,'CONTROL: Could not find data in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif 
;      if(tag_exist(spectrum,'error') eq 0) then begin
;        logprint,'CONTROL: Could not find error in FITS file specified. Assuming it to be zero'
;        sp_error=dblarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'error', sp_error, after='data' 
;      endif 
;      if(tag_exist(spectrum,'dq') eq 0) then begin
;        logprint,'CONTROL: Could not find Data quality in FITS file specified. Assuming data to be good'
;        sp_dq=bytarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'dq', sp_dq, after='error'
;      endif   
;      if(tag_exist(spectrum,'background') eq 0) then begin
;        logprint,'CONTROL: Could not find background in FITS file specified. Assuming it to be zero'
;        sp_background=dblarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'background', sp_background, after='dq'
;      endif
;      if(tag_exist(spectrum,'bck_error') eq 0) then begin
;        logprint,'CONTROL: Could not find background error in FITS file specified. Assuming it to be zero'
;        sp_backgroundr=dblarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'bck_error', sp_backgroundr, after='background'
;      endif
;      if(tag_exist(spectrum,'dq_bg') eq 0) then begin
;        logprint,'CONTROL: Could not find Data qualiy for background in FITS file specified. Assuming it be all good'
;        sp_back_dq=dblarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'dq_bg', sp_back_dq, after='bck_error'
;      endif
;    endelse  
    if (bgs) then begin
      logprint,'CONTROL will carryout background substraction on extracted spectrum now'
      if (idl_ver ge 8) then begin
        if isa(spectrum,'STRUCT') eq 0 then begin 
          logprint,'Seems like extracted spectra were not avaliable assuming the extracted spectra to be in the data folder files.'
          spectrum=mrdfits(rawbc_list[i],1,hdr_nw,/SILENT)
        endif
      endif else begin
        if datatype(spectrum,2) ne 8 then begin
          logprint,'Seems like extracted spectra were not avaliable assuming the extracted spectra to be in the data folder files.'
          spectrum=mrdfits(rawbc_list[i],1,hdr_nw,/SILENT)        
        endif
      endelse  
      im_bflg=fix(SXPAR(hdr_nw, 'BCFLG',MISSING=-1))
      im_dflg=fix(SXPAR(hdr_nw, 'DCFLG',MISSING=-1))
      im_fflg=fix(SXPAR(hdr_nw, 'FCFLG',MISSING=-1))
      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1) then begin
        logprint,'CONTROL: The spectra('+spectra_name[i]+') is not bias or dark or flat corrected with CONTROL.'
        logprint,'Press q to skip the background subtraction for current spectrum. Press any key to continue assuming the spectrum is corrected.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the extraction for current spectrum as requested by the user.'
          goto,spectrum_loop_end
        endif
      endif
      if(tag_exist(spectrum,'data') eq 0) then begin
        errorlog,'CONTROL: Could not find data in FITS file specified. Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      if(tag_exist(spectrum,'error') eq 0) then begin
        errorlog,'CONTROL: Could not find error in FITS file specified. Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif  
      if(tag_exist(spectrum,'dq') eq 0) then begin
        logprint,'CONTROL: Could not find Data quality in FITS file specified. Assuming data to be good'
        errorlog,'CONTROL: Could not find Data quality in FITS file specified. Assuming data to be good',logonly=1
        sp_dq=bytarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'dq', sp_dq, after='error'
      endif
      if(tag_exist(spectrum,'background') eq 0) then begin
        logprint,'CONTROL: Could not find background in FITS file specified. Assuming it to be zero'
        errorlog,'CONTROL: Could not find background in FITS file specified. Assuming it to be zero',logonly=1
        sp_background=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'background', sp_background, after='dq'
      endif
      if(tag_exist(spectrum,'bck_error') eq 0) then begin
        logprint,'CONTROL: Could not find background error in FITS file specified. Assuming it to be zero'
        errorlog,'CONTROL: Could not find background error in FITS file specified. Assuming it to be zero',logonly=1
        sp_backgroundr=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'bck_error', sp_backgroundr, after='background'
      endif
      if(tag_exist(spectrum,'dq_bg') eq 0) then begin
        logprint,'CONTROL: Could not find Data qualiy for background in FITS file specified. Assuming it be all good'
        errorlog,'CONTROL: Could not find Data qualiy for background in FITS file specified. Assuming it be all good',logonly=1
        sp_back_dq=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'dq_bg', sp_back_dq, after='bck_error'
      endif
      if (size(spectrum.data))[0] ne (size(spectrum.background))[0] then begin
        errorlog,'CONTROL: Dimension missmatch for data and background in spectrum file. Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      spectra_1d=spectrum.data-spectrum.background
      if (((size(spectrum.error))[0]) ne ((size(spectrum.bck_error))[0])) then begin
        errorlog,'CONTROL: Dimension missmatch for data error and background error in spectrum file. Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      error_1d=sqrt(spectrum.error^2+spectrum.bck_error^2)
      bgflg=1
      ;data quality
      dq_1d=spectrum.dq+spectrum.dq_bg
      prb=where(dq_1d ge 1)
      if total(prb) ne -1 then dq_1d(prb)=1
      undefine,prb
      if (size(dq_1d))[0] ne 1 then begin
        logprint,'CONTROL: Data quality array is not one diamensional. Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif
      logprint,'CONTROL: Spectrum('+spectra_name[i]+') is now background subtracted.'
      sxaddpar, hdr_nw, 'BGFLG', bgflg ,'BACKGROUND CORRECTION FLAG' 
      hdr_1d=update_header(hdr_nw)
      spectrum_file={data:spectra_1d,error:error_1d,dq:dq_1d}
      if save_temp_files eq 1 then mwrfits,spectrum_file,out_path+spectra_name[i]+'_1d.fits',hdr_1d,/create
    
    endif ;else begin
;      im_bflg=fix(SXPAR(hdr_nw, 'BCFLG',MISSING=-1))
;      im_dflg=fix(SXPAR(hdr_nw, 'DCFLG',MISSING=-1))
;      im_fflg=fix(SXPAR(hdr_nw, 'FCFLG',MISSING=-1))
;      im_eflg=fix(SXPAR(hdr_nw, 'EXTFLF',MISSING=-1))
;      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1 or im_eflg eq -1) then begin
;        logprint,'CONTROL: The spectra('+rawbc_list[i]+') is not bias or dark or flat corrected or extracted with CONTROL.'
;        logprint,'Press q to skip the background subtraction for current spectrum ('+rawbc_list[i]+'). Press any key to continue background correction for the spectrum.'
;        R = GET_KBRD()
;        if R eq 'q' then begin
;          logprint,'CONTROL: Skipping the background subtraction for current spectrum as requested by the user.'
;          goto,spectrum_loop_end
;        endif
;      endif
;      if(tag_exist(spectrum,'data') eq 0) then begin
;        logprint,'CONTROL: Could not find data in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else spectra_1d=spectrum.data
;      if(tag_exist(spectrum,'error') eq 0) then begin
;        logprint,'CONTROL: Could not find error in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else error_1d=spectrum.error
;      if(tag_exist(spectrum,'dq') eq 0) then begin
;        logprint,'CONTROL: Could not find Data quality in FITS file specified. Assuming data to be good'
;        sp_dq=bytarr(n_elements(spectrum.data))
;        struct_add_field, spectrum, 'dq', sp_dq, after='error'
;      endif else dq_1d=spectrum.dq
;      bgflg=0
;    endelse
;    spectrum_file={data:spectra_1d,error:error_1d,dq:dq_1d}
;    sxaddpar, hdr_nw, 'BGFLG', bgflg ,'BACKGROUND CORRECTION FLAG'
;    if(tag_exist(spectrum,'wave') eq 1) then begin
;      struct_add_field, spectrum_file, 'wave', spectrum.wave, before='data'
;    endif  
    ;sxaddpar, hdr_nw,'TFIELDS',tfileds,'Number of fields in each row'
    ;defininig header
    ;sxaddpar, hdr,'TFIELDS',tfileds,'Number of fields in each row'
;    sxaddpar, hdr,'TTYPE1  ','','label for field   1'
;    sxaddpar, hdr,'TFORM1  ',,'Data format of field 1'
;    sxaddpar, hdr,'TUNIT1  ',,'Physical unit of field 1'
;    sxaddpar, hdr,'TDISP1  ',,'Display format for column  1'
;    sxaddpar, hdr,'TNULL1  ',,'Undefined value for column  1'
;    sxaddpar, hdr,'TTYPE2  ',,'label for field   1'
;    sxaddpar, hdr,'TFORM2  ',,'Data format of field 2'
;    sxaddpar, hdr,'TUNIT2  ',,'Physical unit of field 2'
;    sxaddpar, hdr,'TDISP2  ',,'Display format for column  2'
;    sxaddpar, hdr,'TNULL2  ',,'Undefined value for column  2'
;    sxaddpar, hdr,'INHERIT ',,'Inherit the primary header'
;    sxaddpar, hdr,'EXTNAME ',,'Extension name'
;    sxaddpar, hdr,'EXTVER  ',,'Extension version number'
;    sxaddpar, hdr,'ROOTNAME',,'rootname of the observation set'
;    sxaddpar, hdr,'EXPNAME ',,'Exposure identifier'
;    sxaddpar, hdr,'WCSAXES ',,'Number of World Coordinate System axes'
;    sxaddpar, hdr,'RA_APER ',,'RA of aperture reference position'
;    sxaddpar, hdr,'DEC_APER',,'Declination of aperture reference position'
;    sxaddpar, hdr,'PA_APER ',,'Position Angle of reference aperture center'
;    sxaddpar, hdr,'DISPAXIS',,'Dispersion axis;'
;    sxaddpar, hdr,'SHIFTA1 ',,'Spectrum shift in AXIS1 calculated from WAVECAL'
;    sxaddpar, hdr,'ORIENTAT',,'Position angle of image y axis (deg. e of n)'
;    sxaddpar, hdr,'SUNANGLE',,'Angle between sun and V1 axis'
;    sxaddpar, hdr,'MOONANGL',,'Angle between moon and V1 axis'
;    sxaddpar, hdr,'SUN_ALT ',,'Altitude of the sun above Earths limb'
;    sxaddpar, hdr,'FGSLOCK ',,'Commanded FGS lock (FINE,COARSE,GYROS,UNKNOWN)'
;    sxaddpar, hdr,'GYROMODE',,'Number of gyros scheduled, T=3+OBAD'
;    sxaddpar, hdr,'REFFRAME',,'Guide star catalog version'
;    sxaddpar, hdr,'DATE-OBS',,'UT date of start of observation (yyyy-mm-dd)'
;    sxaddpar, hdr,'TIME-OBS',,'UT time of start of observation (hh:mm:ss)'
;    sxaddpar, hdr,'EXPSTART',,'Exposure start time (Modified Julian Date)'
;    sxaddpar, hdr,'EXPEND  ',,'Exposure end time (Modified Julian Date)'
;    sxaddpar, hdr,'EXPTIME ',,'Exposure duration (seconds)--calculated'
;    sxaddpar, hdr,'EXPFLAG ',,'Exposure interruption indicator'
;    sxaddpar, hdr,'QUALCOM1',,'Comment'
;    sxaddpar, hdr,'QUALITY ',,'Data quality'
;    sxaddpar, hdr,'V_HELIO ',,'Heliocentric radial velocity (km/s)'
;    sxaddpar, hdr,'PATTSTEP',,'Position number of this point in the pattern'
;    sxaddpar, hdr,'NCOMBINE',,'Number of image sets combined during CR rejecti'
;    sxaddpar, hdr,'FILLCNT ',,'Number of segments containing fill'
;    sxaddpar, hdr,'ERRCNT  ',,'Number of segments containing errors'
;    sxaddpar, hdr,'PODPSFF ',,'Podps fill present (T/F)'
;    sxaddpar, hdr,'STDCFFF ',,'Science telemetry fill data present (T=1/F=0)'
;    sxaddpar, hdr,'STDCFFP ',,'Science telemetry fill pattern (hex)'
;    sxaddpar, hdr,'CCDTEMP ',,'CCD Temperature'
;    sxaddpar, hdr,'ONTEMP  ',,'On board temperature'
;    sxaddpar, hdr,'NGOODPIX',,'Number of good pixels'
;    sxaddpar, hdr,'SDQFLAGS',,'Serious data quality flags'
;    sxaddpar, hdr,'GOODMIN ',,'Minimum value of good pixels'
;    sxaddpar, hdr,'GOODMAX ',,'Maximum value of good pixels'
;    sxaddpar, hdr,'GOODMEAN',,'Mean value of good pixels'
;    sxaddpar, hdr,'SNRMIN  ',,'Minimum signal to noise of good pixels'
;    sxaddpar, hdr,'SNRMAX  ',,'Maximum signal to noise of good pixels'
;    sxaddpar, hdr,'SNRMEAN ',,'Mean value of signal to noise of good pixels'
;    sxaddpar, hdr,'SOFTERRS',,'Number of soft error pixels (DQF=1)'
;    sxaddpar, hdr,'MEANDARK',,'Average of the dark values subtracted'
;    sxaddpar, hdr,'SPORDER ',,'Spectral order'
;    sxaddpar, hdr,'CONT2EML',,'Intensity conversion: continuum -> emission'
;    sxaddpar, hdr,'SCALE_A1',,'Size of one pixel (arcsec) along dispersion axi'
;    sxaddpar, hdr,'OMEGAPIX',,'Solid angle (arcsec**2) subtended by one pixel'
;    sxaddpar, hdr,'CRSCROFF',,'Offset from 1-D extraction cross-corr.'
 
    ;wavelength calibration
    if (wcl) then begin
      logprint,'CONTROL will wavelength calibrate the 1D science spectrum now'
       if (idl_ver ge 8) then begin 
          if isa(spectrum_file,'STRUCT') eq 0 then spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endif else begin
        if datatype(spectrum,2) ne 8 then spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endelse  
      im_eflg=fix(SXPAR(hdr_1d, 'EXTFLF',MISSING=-1))
      im_bgflg=fix(SXPAR(hdr_1d, 'BGFLG',MISSING=-1))
       if (im_eflg eq -1 or im_bgflg eq -1) then begin
         logprint,'CONTROL: The specified FITS file is not extracted or background subtracted with CONTROL.'
        logprint,'Press q to skip the wavelength calibration for current spectrum ('+rawbc_list[i]+'). Press any key to continue with wavelength calibration.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the wavelength calibration for current spectrum as requested by the user.'
          goto,spectrum_loop_end
        endif
         goto,spectrum_loop_end
       endif
       if(tag_exist(spectrum_file,'data') eq 0) then begin
          logprint,'CONTROL: Could not find data in FITS file specified. Pipeline procedures for this file are aborted'
          goto,spectrum_loop_end
       endif else spectra_1d=spectrum.data
       if(tag_exist(spectrum_file,'error') eq 0) then begin
          logprint,'CONTROL: Could not find error in FITS file specified. Pipeline procedures for this file are aborted'
          goto,spectrum_loop_end
        endif else error_1d=spectrum.error
        if(tag_exist(spectrum_file,'dq') eq 0) then begin
           logprint,'CONTROL: Could not find Data quality in FITS file specified. Assuming data to be good'
           sp_dq=bytarr(n_elements(spectrum_file.data))
           struct_add_field, spectrum_file, 'dq', sp_dq, after='error'
        endif else dq_1d=spectrum_file.dq       
      ;print,spectra_name[i]
      if(tag_exist(infile,'wavecal_mode') eq 1) then wavecal_type=infile.wavecal_mode else wavecal_type='simple'
      if(tag_exist(spectrum_file,'data') eq 0) then begin
        logprint,'CONTROL: Could not find data in FITS file specified. Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif
      if (size(spectrum_file.data))[0] ne 1 then begin
        logprint,'CONTROL: FITS file specified ('+spectra_name[i]+') do not have a 1D spectrum. Pipeline procedures for this file will be aborted'
        goto,spectrum_loop_end
      endif
      wcl_spectrum=control_wavecal(spectrum_file.data,hdr_1d,infile,wavecal_type)
      sp_wavelength=wcl_spectrum.wavelength
      sp_data=wcl_spectrum.flux
      sp_error=spectrum_file.error
      sp_dq=spectrum_file.dq
      wl_hdr=update_header(wcl_spectrum.header)
      wcl_spectrum_file={wave:sp_wavelength,counts:sp_data,error:sp_error,dq:sp_dq}
      cgplot,sp_wavelength,sp_data,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,xtitle='wavelength [$\Angstrom$]',ytitle='flux [counts]'
      write_png,out_path+spectra_name[i]+'_1dw.png',TVRD(/TRUE)
      ;defininig header
;      sxaddpar, hdr,'TTYPE1  ','WAVELENGTH','label for field   1'
;      sxaddpar, hdr,'TFORM1  ',,'Data format of field 1'
;      sxaddpar, hdr,'TUNIT1  ','Angstroms','Physical unit of field 1'
;      sxaddpar, hdr,'TDISP1  ',,'Display format for column  3'
;      sxaddpar, hdr,'TNULL1  ',,'Undefined value for column  1'
;      sxaddpar, hdr,'TTYPE2  ','COUNTS','label for field   2'
;      sxaddpar, hdr,'TFORM2  ',,'Data format of field 2'
;      sxaddpar, hdr,'TUNIT2  ','counts','Physical unit of field 2'
;      sxaddpar, hdr,'TDISP2  ',,'Display format for column  2'
;      sxaddpar, hdr,'TNULL2  ',,'Undefined value for column  2'
;      sxaddpar, hdr,'TTYPE3  ','ERROR','label for field   3'
;      sxaddpar, hdr,'TFORM3  ',,'Data format of field 3'
;      sxaddpar, hdr,'TUNIT3  ','counts','Physical unit of field 3'
;      sxaddpar, hdr,'TDISP3  ',,'Display format for column  3'
;      sxaddpar, hdr,'TNULL3  ',,'Undefined value for column  3'

       mwrfits,wcl_spectrum_file,out_path+spectra_name[i]+'_1dw.fits',wl_hdr, /CREATE
       spectra_1dw  = sp_data
       error_1dw    = sp_error
       dq_1dw       = sp_dq
    endif ;else begin
;      spectrum_wl=mrdfits(rawbc_list[i],1,wl_hdr)
;      if(tag_exist(spectrum_wl,'wave') eq 0) then begin
;        logprint,'CONTROL: Could not find wavelength information in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else sp_wavelength=spectrum_wl.wave
;      if(tag_exist(spectrum_wl,'counts') eq 0) then begin
;        logprint,'CONTROL: Could not find count information in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else spectra_1dw=spectrum_wl.counts
;      if(tag_exist(spectrum_wl,'error') eq 0) then begin
;        logprint,'CONTROL: Could not find error information in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else error_1dw=spectrum_wl.error
;      if(tag_exist(spectrum_wl,'dq') eq 0) then begin
;        logprint,'CONTROL: Could not find Data quality information in FITS file specified. Pipeline procedures for this file is aborted'
;        goto,spectrum_loop_end
;      endif else dq_1dw=spectrum_wl.dq
;      
;    endelse
      ;flux calibration
      if (fcl) then begin
        logprint,'CONTROL will flux calibrate the wavelength calibrated 1D science spectrum now'
       if (idl_ver ge 8) then begin 
          if isa(wcl_spectrum_file,'STRUCT') eq 0 then wcl_spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
       endif else begin
          if datatype(spectrum,2) ne 8 then wcl_spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
       endelse
       if  fix(SXPAR(wl_hdr, 'WCALFLG',MISSING=-1)) eq -1 then begin
          logprint,'CONTROL: The FITS file specified('+rawbc_list[i]+') is not wavelength calibrated. Pipeline procedures for this file will be aborted'
         goto,spectrum_loop_end
       endif
       if tag_exist(infile,'flux_cal') eq 0 then begin
         logprint,'CONTROL:  Flux calibration file not found. Please re-run the pipeline with the file.',logonly=logonly
         message,'CONTROL: Flux calibration file not found. Please re-run the pipeline with actual file.'
       endif
        flux_calib_file=detectos(infile.flux_cal)
        if (file_test(flux_calib_file) eq 0) then begin
          logprint,'CONTROL: Flux calibration file not found. Please re-run the pipeline with actual file address',logonly=logonly
          message,'CONTROL: Flux calibration file not found. Please re-run the pipeline with actual file address'
        endif
        readcol,flux_calib_file,fcal_wave,fcal_value,F='D,D' ;flux calibration file file location
        if(tag_exist(wcl_spectrum_file,'wave') eq 0) then begin
           logprint,'CONTROL: Could not find wavelength information in FITS file specified. Pipeline procedures for this file are aborted'
           goto,spectrum_loop_end
        endif else sp_wavelength=wcl_spectrum_file.wave
        if(tag_exist(wcl_spectrum_file,'counts') eq 0) then begin
            logprint,'CONTROL: Could not find count information in FITS file specified. Pipeline procedures for this file are aborted'
             goto,spectrum_loop_end
        endif else spectra_1dw=wcl_spectrum_file.counts
        if(tag_exist(wcl_spectrum_file,'error') eq 0) then begin
            logprint,'CONTROL: Could not find error information in FITS file specified. Pipeline procedures for this file are aborted'
            goto,spectrum_loop_end
        endif else error_1dw=wcl_spectrum_file.error
        if(tag_exist(wcl_spectrum_file,'dq') eq 0) then begin
            logprint,'CONTROL: Could not find Data quality information in FITS file specified. Pipeline procedures for this file are aborted'
            goto,spectrum_loop_end
        endif else dq_1dw=wcl_spectrum_file.dq
        QE=interpol(fcal_value,fcal_wave,sp_wavelength,/SPLINE)
        spectra_1dwf=spectra_1dw*QE
        error_1dwf=error_1dw*QE
        wclf_spectrum_file={wave:sp_wavelength,flux:spectra_1dwf,error:error_1dwf,dq:dq_1dw}
        cgplot,sp_wavelength,spectra_1dwf,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,xtitle='wavelength [$\Angstrom$]',ytitle='flux [ergs s!u-1!n cm!u-2!n $\Angstrom$!u-1!n]'
        write_png,out_path+spectra_name[i]+'_1dwf.png',TVRD(/TRUE)
        spectrum_hdr=update_header(wl_hdr)
        sxaddpar, spectrum_hdr, 'FCALFLG', 1 ,'FLUX CORRECTION FLAG'
        ;defininig header
;        sxaddpar, hdr,'TTYPE1  ','WAVELENGTH','label for field   1'
;        sxaddpar, hdr,'TFORM1  ',,''
;        sxaddpar, hdr,'TUNIT1  ','Angstroms','Physical unit of field 1'
;        sxaddpar, hdr,'TDISP1  ',,'Display format for column  3'
;        sxaddpar, hdr,'TNULL1  ',,'Undefined value for column  1'
;        sxaddpar, hdr,'TTYPE2  ',,'FLUX'
;        sxaddpar, hdr,'TFORM2  ',,'Data format of field 2'
;        sxaddpar, hdr,'TUNIT2  ','erg/s/cm**2/Angstrom','Physical unit of field 2'
;        sxaddpar, hdr,'TDISP2  ',,'Display format for column  2'
;        sxaddpar, hdr,'TNULL2  ',,'Undefined value for column  2'
;        sxaddpar, hdr,'TTYPE3  ',,'FLUX ERROR'
;        sxaddpar, hdr,'TFORM3  ',,'Data format of field 3'
;        sxaddpar, hdr,'TUNIT3  ','erg/s/cm**2/Angstrom','Physical unit of field 3'
;        sxaddpar, hdr,'TDISP3  ',,'Display format for column  3'
;        sxaddpar, hdr,'TNULL3  ',,'Undefined value for column  3'
         mwrfits,wclf_spectrum_file,out_path+spectra_name[i]+'_1dwf.fits',spectrum_hdr, /CREATE
      endif 

spectrum_loop_end:
 endfor 
 
 spectrums=file_search(out_path+'*_1dwf.fits')
 if (science ne 0) then begin
 cute_light_curve,spectrums,lightcurve,wave_region=['MgI','MgII','FeII']
;
;           sxaddpar, lc_hdr,'TTYPE1  ','TIME','label for field   1'
;           sxaddpar, lc_hdr,'TFORM1  ',,''
;           sxaddpar, lc_hdr,'TUNIT1  ','Seconds','Physical unit of field 1'
;           sxaddpar, lc_hdr,'TDISP1  ',,'Display format for column  3'
;           sxaddpar, lc_hdr,'TNULL1  ',,'Undefined value for column  1'
;           sxaddpar, lc_hdr,'TTYPE2  ',,'FLUX'
;           sxaddpar, lc_hdr,'TFORM2  ',,'Data format of field 2'
;           sxaddpar, lc_hdr,'TUNIT2  ','erg/s/cm**2/Angstrom','Physical unit of field 2'
;           sxaddpar, lc_hdr,'TDISP2  ',,'Display format for column  2'
;           sxaddpar, lc_hdr,'TNULL2  ',,'Undefined value for column  2'
;           sxaddpar, lc_hdr,'TTYPE3  ',,'FLUX ERROR'
;           sxaddpar, lc_hdr,'TFORM3  ',,'Data format of field 3'
;           sxaddpar, lc_hdr,'TUNIT3  ','erg/s/cm**2/Angstrom','Physical unit of field 3'
;           sxaddpar, lc_hdr,'TTYPE4  ',,'FLUX'
;           sxaddpar, lc_hdr,'TFORM4  ',,'Data format of field 2'
;           sxaddpar, lc_hdr,'TUNIT4  ','erg/s/cm**2/Angstrom','Physical unit of field 2'
;           sxaddpar, lc_hdr,'TDISP4  ',,'Display format for column  2'
;           sxaddpar, lc_hdr,'TNULL4  ',,'Undefined value for column  2'
;           sxaddpar, lc_hdr,'TTYPE5  ',,'FLUX ERROR'
;           sxaddpar, lc_hdr,'TFORM5  ',,'Data format of field 3'
;           sxaddpar, lc_hdr,'TUNIT5  ','erg/s/cm**2/Angstrom','Physical unit of field 3'
;           sxaddpar, lc_hdr,'TTYPE6  ',,'FLUX'
;           sxaddpar, lc_hdr,'TFORM6  ',,'Data format of field 2'
;           sxaddpar, lc_hdr,'TUNIT6  ','erg/s/cm**2/Angstrom','Physical unit of field 2'
;           sxaddpar, lc_hdr,'TDISP6  ',,'Display format for column  2'
;           sxaddpar, lc_hdr,'TNULL6  ',,'Undefined value for column  2'
;           sxaddpar, lc_hdr,'TTYPE7  ',,'FLUX ERROR'
;           sxaddpar, lc_hdr,'TFORM7  ',,'Data format of field 3'
;           sxaddpar, lc_hdr,'TUNIT7  ','erg/s/cm**2/Angstrom','Physical unit of field 3'

   mwrfits,lightcurve,out_path+'light_curve.fits',lc_hdr, /CREATE

;   trlen=n_elements(spectra)
;   trlenby4=fix(trlen/4)
;   ph_st1=lc1[1,0:trlenby4]
;   ph_en1=lc1[1,nx-trlenby4:nx-1]
;   out_tran1=[ph_st1,ph_en1]
;   mean_ph1=mean(out_tran1)
;   tlc1=lc1[0,*]
;   nlc1=lc1[1,*]/mean_ph1
;   elc1=lc1[2,*]/mean_ph1
;   
;   ph_st2=lc2[1,0:trlenby4]
;   ph_en2=lc2[1,nx-trlenby4:nx-1]
;   out_tran2=[ph_st2,ph_en2]
;   mean_ph2=mean(out_tran2)
;   tlc2=lc2[0,*]
;   nlc2=lc2[1,*]/mean_ph2
;   elc2=lc2[2,*]/mean_ph2
;  
;   ph_st3=lc3[1,0:trlenby4]
;   ph_en3=lc3[1,nx-trlenby4:nx-1]
;   out_tran3=[ph_st3,ph_en3]
;   mean_ph3=mean(out_tran3)
;   tlc3=lc3[0,*]
;   nlc3=lc3[1,*]/mean_ph3
;   elc3=lc3[2,*]/mean_ph3
;   !P.Multi=[0,1,2]
;   cgplot,tlc1,nlc1,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc1,$
;     ERR_YHigh=elc1,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Short wavelengths'
;     
;     cgplot,tlc2,nlc2,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc2,$
;     ERR_YHigh=elc2,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Middle wavelengths'
;     
;   cgplot,tlc3,nlc3,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc3,$
;     ERR_YHigh=elc3,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Long wavelengths'
;     write_png,out_path+'light_curve.png',TVRD(/TRUE)
;     !P.Multi=0
 endif  
 code_end: 
 logprint,'CONTROL wrapping up the data reduction procedure.'
 logprint,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',close=close
close,/all
end
