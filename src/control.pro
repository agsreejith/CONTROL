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
;   Parameter file: Parameter file provided with the distribution providing details defining the 
;   pipeline operation.
;
; OUTPUT:
;   Produce science quality data products for CUTE mission (1D flux and wavelength calibrated 
;   spectra, 2D calibrated spectra, Transit light curves).
;
; PROCEDURE:
;
;data quality implemented till lc
;;#################################################################################################


pro control,file,help=help
  
  !quiet = 1
  close,/all
  
  ;display help for the program
  if keyword_set(help) then begin
    print,'                            CUTE AUTONOMOUS DATA REDUCTION PIPELINE'
    print,''
    print,'***************************************************************************************'
    print,''
    print,' This software is intented to be fully automated, aimed at producing science-quality' 
    print,' output with a single command line with zero user interference for CUTE data.'
    print,' It can be easily used for any single order spectral data in any wavelength without any'
    print,' modification.'
    print,''
    print,' Options avaliable in configuration file'
    print,'01. Path locations: Data path, intermediate file path and output file path.'
    print,'    If data path is not provided it is assumed to be the current directory.'
    print,'    If intermediate file path and output file path are not provided they are created in'
    print,'    the data path.'
    print,'02. The user has the option to select individual steps during reduction. The options are'
    print,'    as follows'
    print,'    hb_correction,create_mbias,create_mdark,create_mflat,cr_bias,cr_dark,cr_flat'
    print,'    cr_cosmic,extract,bg_sub,wcalib,fluxcalib,light_curve,retrieval,level3,all'
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
    print,'    light_curve: Create seven light curves form the data.'
    print,'    retrieval: Process the flux and wavelength calibrated spectra to obtain transmission'
    print,'    spectra.'
    print,'    level3: Carry out all processes required for generating Level 3 data of CUTE, ie., 
    print,'    flux and wavelength calibrated 1-D spectra.'
    print,'    all: Execute all the steps. This is also the default option if steps keyword is'
     print,'   absent.'
    print,'03. Location of hot and bad pixel map. The location of the map of cosmetic defects on '
    print,'    the CCD as a FITS file with the same size as the input CCD frame'
    print,'    byte-type pixels (8-bit) and values of 0B for good pixels and 1B for bad ones.'
    print,'04. Location of master calibration files if available.'
    print,'    This includes locations of master bias, master dark and/or master flat files '
    print,'    (as fits files).'
    print,'05. Set save_temp_files equal to 1 if you need all the intermediate files to be saved to'
    print,'    the intermediate file directory. Set it to 0 if you do not need them. Default is 1.'
    print,'06. Options for master file creation.'
    print,'    Statistical method to create master files. Options are: median, mean or mode'
    print,'    Threshold for rejection of pixels: sigma_deviation (how many sigmas from mean value'
    print,'    of the frame).'
    print,'    Saturation_level sets the saturation limit of CCD. Default value is 72000.'
    print,'07. Cosmic Ray correction options: Level of cr clipping for LA cosmic algorithm. 
    print,'    Default value used is 8.'
    print,'08. Trace parameters'
    print,'    Degree for centroid polynomial: Default is 1.'
    print,'    Trace type used for extraction. The options are: simple,fixed,variable,function.'
    print,'    For details of these options refer to the manual.'
    print,'    Additional parameters required for different trace options:'
    print,'     centroid  : Centroid will be calculated if not provided'
    print,'     slope     : Required for option: simple or fixed, defines the slope of the spectrum'
    print,'     width     : Required for option: simple, which defines the width of the spectrum'
    print,'     upper     : Required for option: fixed, defines the upper width of the spectrum'
    print,'     lower     : Required for option: fixed, defines the lower width of the spectrum'
    print,'     threshold : Required for option: variable or function, defines the threshold from '
    print,'                 the maximum/peak of the spectrum.'
    print,'    File location for information of where to extract background. The user could also '
    print,'    provide it as fixed number which corresponds to the shift in pixels from centroid.'
    print,'09. Wavelength calibration variables'
    print,'    Type of wavelength calibration required: Options are simple and crscor'
    print,'    Location of wavelength file. A text file with length equal to number of pixels in '
    print,'    cross dispersion direction representing wavelength to pixel mapping.'
    print,'    A location of look up table for stellar parameters. Required only for the option'
    print,'    crscor so as to compare with a model data.'
    print,'    If this look up table is not provided then the model file should be a two column'
    print,'    model file of wavelength vs flux.'
    print,'    Location of synthetic spectrum (model_file) for cross-correlation. CONTROL assumes'
    print,'    that folders are named by their temperature and files are named as model.flx or'
    print,'    a two column model file of wavelength vs flux, if stellar parameters look up table'
    print,'    is not provided.'
    print,'10. Flux calibration variables'
    print,'    Location of flux calibration file which provides wavelength vs response relation.'
    print,'To avoid problematic files from processing, create badfiles.txt in the data folder.'
    print,'***************************************************************************************'
    print,'Operation steps of pipeline are as follows'
    print,' The program works in a series of steps following standard CCD reduction techniques the
    print,' user can also select individual modules using the step function described above.'
    print,' A reduction log is created to help the users with processes carried out and mentioning'
    print,' the different parameters used if any.'
    print,' It also creates an errorlog, listen all errors that occured during the process.'
    print,' The steps involved in the program execution are as follows:'
    print,''
    print,'   1.  Prepare: Check for the different input parameters are set variables accordingly.'
    print,'       If data files are present classify them according to file type.'
    print,'   2.  Hot and bad pixel correction: Correct for hot and bad pixels in the frames.'
    print,'   3.  Create master bias: Create a master bias if master bias file is not found from a'
    print,'       set of bias files.'
    print,'   4.  Create master dark: Create a master dark if master dark file is not found from a '
    print,'       set of dark files.'
    print,'   5.  Create master flat: Create a master flat if master flat file is not found from a '
    print,'       set of flat files.'
    print,'   6.  Correct for bias: Correction of the effect of bias in science frames.'
    print,'   7.  Correct for dark: Correction of the effect of dark in science frames.'
    print,'   8.  Correct for flat: Correction of the effect of flat in science frames.'
    print,'   9.  Correct for cosmic rays: Correction for the effect of cosmic rays using LA cosmic.'
    print,'   10. Extract spectrum: Define spectrum trace and extract the spectrum.'
    print,'   11. Subtract Background: Subtract background from spectrum.'
    print,'   12. Wavelength calibration: Do wavelength calibration'
    print,'   13. Flux calibration: Do flux calibration'
    print,'   14. Default Light curve: Create 3 light curves (short, middle and long wavelength)'
    print,'   15. Transmission spectrum: Retrieve transmission spectrum.'
    print,'***************************************************************************************'
    print,'For additional help and more options, go to: https://github.com/agsreejith/CONTROL'
    return
  endif
  
  ;check if parameter file is provided with the procedure call.
  if keyword_defined(file) eq 0 then begin
    CD, Current=cur_dir                                ;check current directory for parameter file.
    file=cur_dir+'control_parameters.txt'
    print,'Configuration file assumed to be in current directory'
  endif
  
  ;display parameter file format if not provided by the user.
  if file_test(file) eq 0 then begin
    print,'CONTROL: Missing configuration file. Please re-run CONTROL with a configuration file'
    print,'CONTROL: Users should have received a configuration file with this software distribution.'
    print,'if not please refer to the configuration file format below.'
    print,'CONTROL: The configuration file should have the following parameters'
    print,'#######################################################################################'
    print,'#'
    print,'#      CONFIGURATION FILE FOR CONTROL v1.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE'
    print,'#'
    print,'#**************************************************************************************'
    print,'print,#Please note: Comments and empty lines have to start with #.'
    print,'#Use # also for parameters not used.'
    print,'#'
    print,'#**************************************************************************************'
    print,'#'
    print,'data_path       = '
    print,'temp_path      ='
    print,'out_path       ='
    print,'steps           = '
    print,'hb_map          = '
    print,'master_bias_file   ='
    print,'master_dark_file   ='
    print,'master_flat_file   ='
    print,'save_temp_files     = '
    print,'bias_combine_type   = '
    print,'dark_combine_type   = '
    print,'flat_combine_type   = '
    print,'saturation_limit  = '
    print,'sigma_deviation     = '
    print,'cosmic_ray_clip    ='
    print,'cent_poly_degree   ='
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
    print,'#**************************************************************************************'
    print,'Run CONTROL with the option help to get the details about configuration file keywords '
    print,'or please refer to the software manual.'
    message,'CONTROL: Exiting as configuration file is not avaliable'
  endif
  
  ;read the input parameter file.
  infile=gm_read_textstructure(file)
  
  ;set data path flag to 0.
  dpflg=0
  
  ;get system time of the the program run.
  t=systime(/UTC)
  tsec=systime(/seconds)
  
  ;check for the steps keyword in parameter file and set the necessary flags in the program.
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
    ;hb_correction,create_mbias,create_mdark,create_mflat,cr_bias,cr_dark,cr_flat,cr_cosmic,extract,
    ;wcalib,fluxcalib,light_curve,cr_systematics
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
  
  ;get idl version to enable switch between code segments when required.
  idl_ver=float(!Version.RELEASE)
  
  ;reading parameters and set defauts where required.
  if(tag_exist(infile,'data_path') eq 1) then data_path=infile.data_path
  if(tag_exist(infile,'temp_path') eq 1) then inter_path=infile.temp_path
  if(tag_exist(infile,'out_path') eq 1) then out_path=infile.out_path
  if(tag_exist(infile,'hbpix_corection_type') eq 1) then hb_type=infile.hbpix_corection_type $ 
  else hb_type='interpolate'
  if(tag_exist(infile,'bias_combine_type') eq 1) then b_type=infile.bias_combine_type
  if(tag_exist(infile,'dark_combine_type') eq 1) then d_type=infile.dark_combine_type
  if(tag_exist(infile,'flat_combine_type') eq 1) then f_type=infile.flat_combine_type
  if(tag_exist(infile,'cosmic_ray_clip') eq 1) then crclip=infile.cosmic_ray_clip else crclip = 8
  if(tag_exist(infile,'cosmic_before_dark') eq 1) then cscr_bfd =fix(infile.cosmic_before_dark)$ 
  else cscr_bfd = 0
  if(tag_exist(infile,'save_temp_files') eq 1) then save_temp_files=infile.save_temp_files $
  else save_temp_files = 1
  if(tag_exist(infile,'hb_map') eq 1) then hbmask=infile.hb_map
  if(tag_exist(infile,'saturation_limit') eq 1) then sat_value=infile.saturation_limit
  if(tag_exist(infile,'sigma_deviation') eq 1) then threshold=infile.sigma_deviation
  
  ;flag for step that determine light curve is reassigned to variable science.
  science = sci_out
  
  ;check if path location for data is specifed.
  if n_elements(data_path) eq 0 then begin
    dpflg=1
    CD, Current=data_path                ;Assuming data to be current directory. 
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': data_path=data_path+'\' ;WINDOWS
      'UNIX': data_path=data_path+'/'    ; UNIX.
    ENDCASE
  endif else begin
    data_path=detectos(data_path) 
    if (file_test(data_path,/DIRECTORY)eq 0) then begin
      dpflg=2
      FILE_MKDIR,data_path                ;create data path directory if not found.
    endif
  endelse
  
  ;checks for output path.
  if n_elements(out_path) eq 0 then begin
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': out_path=data_path+'output\' ;WINDOWS
      'UNIX': out_path=data_path+'output/'; UNIX.
    ENDCASE
    dpflg=3
   endif else out_path=detectos(out_path)

  if (file_test(out_path,/DIRECTORY)eq 0) then begin
    dpflg=4
    FILE_MKDIR,out_path
  endif
  
  logprint,' _______  _______  __    _  _______  ______    _______  ___',$
  logfile=out_path+'control_log'+string(tsec,Format='(I12)')+'.txt'
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

  errorlog,' _______  _______  __    _  _______  ______    _______  ___',$
  logfile=out_path+'control_error_log'+string(tsec,Format='(I12)')+'.txt',logonly=1
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
  
  ;Basic error and log information for missing/unavaliable files
  if dpflg eq 1 then $
    errorlog,'CONTROL: No data file path found. Data files assumed to be in current directory'
  if dpflg eq 2 then $
    logprint,'CONTROL: No data file path found. Creating data file path based on confi file'
  if dpflg eq 3 then $
    logprint,'CONTROL: No output file path found. Output file directory created in data directory'
  if dpflg eq 4 then $
    errorlog,'CONTROL: No output file path found. Creating output file path based on config file'
  
  ;path verifications
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
    errorlog,'CONTROL: No temporary file path found.'$
            +' Temporary file directory created in data directory'
  endif else inter_path=detectos(inter_path)

  if n_elements(out_path) eq 0 then begin
    CASE StrUpCase(!Version.OS_Family) OF
      'WINDOWS': out_path=data_path+'output\' ;WINDOWS
      'UNIX': out_path=data_path+'output/'; UNIX.
    ENDCASE
    errorlog,'CONTROL: No output file path found.'$
            +' Output file directory created in data directory'
  endif else out_path=detectos(out_path)

  ;creates paths if they dont exist.
  if (file_test(data_path,/DIRECTORY)eq 0) then begin
    errorlog,'CONTROL: No data file path found.'$
             +' Creating data file path based on configuration file'
    FILE_MKDIR,data_path
  endif
  if (file_test(inter_path,/DIRECTORY)eq 0) then begin
    errorlog,'CONTROL: No temporary file path found.'$
             +' Creating temporary file path based on configuration file'
    FILE_MKDIR,inter_path
  endif
  if (file_test(out_path,/DIRECTORY)eq 0) then begin
    errorlog,'CONTROL: No output file path found.'$
             +' Creating output file path based on configuration file'
    FILE_MKDIR,out_path
  endif
  if save_temp_files lt 0 or save_temp_files gt 1 then begin
    logprint,'Allowed values for save temporary files are 1 or 0.'$
             +' For any other values default of 1 will be assumed.'
    save_temp_files=1
  endif
  
  ;Get file list from files in data path.
  all_file_list=file_search(data_path+'*.fits') 
  ;unique files?
  if (n_elements(all_file_list) eq 0) then begin
    errorlog,'CONTROL: File acess error, No fits file found!!',logonly = logonly
    errorlog,'CONTROL: Verify data path.'
    message,'File acess error, No fits file found for reduction!!'
    return
  endif
  ;check for a badfiles from a list
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

  ;carry out bad and hot pixel corrections on good files
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
        sxaddpar, hdr, 'HBFLG', 1.,'Hot and bad pixel correction flag'
        sxaddpar, hdr, 'HBMASK', hbmask,'Location of hot and bad pixel mask file'
        sxaddpar, hdr, 'HBTYPE',hb_type,'Type of hot and bad pixel correction employed.
        sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE v1.0 .',hdr,/COMMENT
        writefits,inter_path+file_names[i]+'_hb.fits',out_im,hdr ;hotbad file
        file_type[i]=SXPAR(hdr, 'FILETYPE',MISSING=-1);Find object description
      endfor
      file_list_hb=file_search(inter_path+'*_hb.fits')
    endif else begin
      logprint,'No input mask file found for hot and bad pixel correction,'$
               +' proceeding without hot and bad pixel correction.'
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
    ;Logging just as a error, will not cause program to stop.
  endelse
  
  ;Compile a list of DARK frames
  darklist_ind=where(strpos(file_type,'DARK') ge 0)
  if darklist_ind[0] ne -1 then darklist=file_list_hb[darklist_ind] else begin
    errorlog,'CONTROL: No DARK files found.'
    ;message,'CONTROL: No DARK files found. CONTROL will exit now.' 
    ;Logging just as a error, will not cause program to stop.
  endelse
  
  ;Compile a list of FLAT frames
  flatlist_ind=where(strpos(file_type,'FLAT') ge 0)
  if flatlist_ind[0] ne -1 then flatlist=file_list_hb[flatlist_ind] else begin
    errorlog,'CONTROL: No FLAT files found.'
    ;message,'CONTROL: No FLAT files found. CONTROL will exit now.' 
    ; Logging just as a error, will not cause program to stop.
  endelse
  
  
  ;Get list of science frames
  spectra_ind=where(strpos(file_type, 'OBJECT') ge 0)
  if spectra_ind[0] ne -1 then begin
    spectra=file_list_hb[spectra_ind]
    spectra_name=file_basename(spectra,'.fits')
    spectra_bname=file_basename(spectra,'_hb.fits')
  endif else begin
    errorlog,'CONTROL: No SCIENCE files found.'
    ;message,'CONTROL: No SCIENCE files found. CONTROL will exit now.' 
    ; Logging just as a error, will not cause program to stop.
  endelse
  
  ;calculate readout noise and gain
;  if n_elements(biaslist) eq 0 then $
;    logprint,'CONTROL: Cannot calculate readnoise and gain without bias files.'$
;             +' These values will be read from headers of master files'$
;    else begin
;    b1=mrdfits(biaslist[0],0,hdr,/SILENT)
;    b2=mrdfits(biaslist[1],0,hdr,/SILENT)
;    if n_elements(flatlist) eq 0 then begin
;      logprint,'CONTROL: Cannot calculate gain without flat files.'$
;               +' These values will be read from headers of master files'
;      if tag_exist(infile,'master_flat_file') eq 0 then $
;        logprint,'CONTROL: No master flat file found to read gain values,'$
;                 +' reading them from bias files'
;      g_calc=SXPAR( hdr, 'CCDGAIN')
;    endif else begin
;      f1=mrdfits(flatlist[0],0,hdr,/SILENT)
;      f2=mrdfits(flatlist[1],0,hdr,/SILENT)
;      g_calc=gain(f1,f2,b1,b2)
;      g_calc=1
;    endelse
;    R=rnoise(b1,b2,g_calc)
;  endelse
;  for i=0,n_elements(biaslist)-1 do begin
;    h = headfits(biaslist[i])
;    sxaddpar,h,'RNOISE',R
;    sxaddpar,h,'CCDGAIN',g_calc
;    modfits,biaslist[i],0,h
;  endfor

  ;bias section
  if (crmb) then begin
    logprint,'CONTROL will create master bias now'
    if tag_exist(infile,'master_bias_file') eq 0 then begin
      if n_elements(biaslist) eq 0 then begin
        errorlog,'CONTROL: File access error, No master bias or bias file found',logonly = logonly
        logprint,'CONTROL: Verify data path or check bias file header.'$
                 +' CONTROL checks for FILETYPE keyword in header to clasify files'
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
        logprint,'CONTROL: Verify data path or check dark file header.'$
                 +' CONTROL checks for FILETYPE keyword in header to clasify files.'
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
        if extension gt 1 then md_dq=mrdfits(mdfile,2,dhdr3,/SILENT) $
        else md_dq=bytarr(nxyd[1],nxyd[2])
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
        if extension gt 1 then md_dq=mrdfits(mdfile,2,dhdr3,/SILENT) $
        else md_dq=bytarr(nxyd[1],nxyd[2])
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
        logprint,'CONTROL: Verify data path or check flat file header.'$
                 +' CONTROL checks for FILETYPE keyword in header to clasify files.'
        message,'CONTROL: File access error, No master flat or flat file found!'
      endif
      control_flat_combine,flatlist,mflat,mbias,mdark,$
                           type=f_type,sat_value=sat_value,threshold=threshold
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
        if extension gt 1 then mf_dq=mrdfits(mdfile,2,hdr,/SILENT) $
        else mf_dq=bytarr(nxyd[1],nxyd[2])
      endelse
    endelse
  endif else begin
    if tag_exist(infile,'master_flat_file') ne 0 then begin
      mffile_path=detectos(infile.master_flat_file)
      mffile=file_search(mffile_path)
      if (file_test(mffile_path) eq 0) then begin
        logprint,'CONTROL: Verify data path or check bias file header.'$
                 +' CONTROL checks for FILETYPE keyword in header to clasify files.'
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
        if extension gt 1 then mf_dq=mrdfits(mdfile,2,hdr,/SILENT) $
        else mf_dq=bytarr(nxyd[1],nxyd[2])
      endelse
    endif
  endelse
  ;bias correction
  if datatype(mbias,2) ne 8 then begin
    logprint,'CONTROL could not find master bias file.'$
             +' Do you want to continue reducing this spectrum assuming bias to be zero'
    logprint,'Press q to exit. Press any key to use bias as zero.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master bias of zeros (2048x512) and assuming readout noise to be 12.25'
    sxaddpar,mbhdr,'RNOISE',12.25,'Read out noise'
    mbias={im:dblarr(2048,512),error:dblarr(2048,512),hdr:mbhdr,dq:bytarr(2048,512)}
    mb_dq=bytarr(2048,512)

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
          logprint,'CONTROL: The spectrum('+spectra[i]+') is already bias corrected.'$
                   +' Skipping the bias correction for current spectrum.'
          goto,spectrum_bias_cr_end
        endif
      endif
      mnx=(size(mbias.im))[1]
      mny=(size(mbias.im))[2]
      ycut1=SXPAR( hdr, 'YCUT1')
      ycut2=SXPAR( hdr, 'YCUT2')
      ylen=ycut2-ycut1+1
      ccd_gain=SXPAR( hdr, 'GAIN')
      if ccd_gain le 0 then ccd_gain=1
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
      if nx eq mnx then begin
        if ny eq mny then begin
          raw_imb=(raw_im-mbias.im)/ccd_gain
          sigma_imb=sqrt(raw_im+r^2)/ccd_gain
          ;data quality
          dq_imb=dq_im+mb_dq
          logprint,'CONTROL: Bias correction carried out on spectrum('+spectra[i]+').'
        endif else begin
          if (ylen ne ny) then mbias_nw=mbias.im[*,ycut1:ycut1+ny-1] $
          else mbias_nw=mbias.im[*,ycut1:ycut2]
          if (ylen ne ny) then mbdq_nw=mb_dq[*,ycut1:ycut1+ny-1] $
          else mbdq_nw=mb_dq[*,ycut1:ycut2]
          raw_imb=(raw_im-mbias_nw)/ccd_gain
          sigma_imb=sqrt(raw_im+r^2)/ccd_gain
          ;data quality
          dq_imb=dq_im+mbdq_nw
          logprint,'CONTROL: Bias correction carried out on spectrum('+spectra[i]+').'
        endelse
      endif else begin
        errorlog,'CONTROL: X Dimension error with science image and bias image'
        return
      endelse
      ;data quality
      prb=where(dq_imb ge 1)
      if total(prb) ne -1 then dq_imb[prb]=1
      nbad=n_elements(prb)
      dqfactr=float(nbad/n_elements(dq_imb))
      ;headers
      sxaddpar, hdr, 'BCFLG', 1,'BIAS CORRECTION FLAG' ;bias correction flag
      sxaddpar, hdr, 'SIGMBIAS', SXPAR(mbias.hdr,'SIGMBIAS'), $
                'Standard deviation of the master bias frame'
      sxaddpar, hdr, 'MEANBIAS', SXPAR(mbias.hdr,'MEANBIAS'), $
                'Mean value of the master bias frame'
      sxaddpar, hdr, 'MDNBIAS ', SXPAR(mbias.hdr,'MDNBIAS'), $
                'Median value of the master bias frame'
      sxaddpar, hdr, 'MAXBIAS ', SXPAR(mbias.hdr,'MAXBIAS'), $
                'Maximum value of the master bias frame'
      sxaddpar, hdr, 'MINBIAS ', SXPAR(mbias.hdr,'MINBIAS'), $
                'Minimum value of the master bias frame'
      sxaddpar, hdr, 'NFRAMEMB', SXPAR(mbias.hdr,'NFRAMES'), $
                'Number of frames used in bias combine'
      sxaddpar, hdr, 'BIASTYP ', SXPAR(mbias.hdr,'BIASTYP'), $
                'Type of bias combine employed'
      sxaddpar, hdr, 'BIASSAT ', SXPAR(mbias.hdr,'BIASSAT'), $
                'Saturation limit used in bias frames'
      sxaddpar, hdr, 'BIASSIG ', SXPAR(mbias.hdr,'BIASSIG'), $
                'Deviation limit for good bias frames'
      sxaddpar, hdr, 'NBADPIX ',nbad,'Number of bad data quality pixels'
      sxaddpar, hdr, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
                ' fraction of total number of pixels'
      sxdelpar,hdr,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0.',hdr,/COMMENT
             
      mwrfits,raw_imb,inter_path+spectra_name[i]+'_b.fits',hdr,/create
      mwrfits,sigma_imb,inter_path+spectra_name[i]+'_b.fits',hdr
      mwrfits,dq_imb,inter_path+spectra_name[i]+'_b.fits',hdr
      undefine,prb,nbad,dqfactr
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
          logprint,'CONTROL: The spectrum('+spectra[i]+') is already corrected for cosmic rays.'$
                   +' Skipping cosmic ray correction for current spectrum.'
          goto,spectrum_cr_cr_end
        endif
      endif
      ;cosmic ray correction
      logprint,'CONTROL: Carrying out cosmic ray correction with clip value of'$
                +strtrim(string(crclip),2)+'.'
      la_cosmic,rawb_list,outlist=rawbc_list,masklist=mask_list,$
                gain=ccd_gain,readn=r,sigclip=crclip; add other parameters after testing
      crflg=1
      for i=0,n_elements(rawbc_list)-1 do begin
        raw_imbc=mrdfits(rawbc_list[i],0,hdr)
        nxc=(size(raw_imbc))[1]
        nyc=(size(raw_imbc))[2]
        raw_imb=mrdfits(rawb_list[i],0,hdr)
        sxaddpar, hdr,'CRFLG',crflg,'COSMIC RAY CORRECTION FLAG'
        sxaddpar, hdr,'CSRYWEN', 'before dark','When cosmic ray correction was carried out'
        sxaddpar, hdr,'CRCYCLP', crclip,'Clip value used for cosmic ray'
        sxdelpar, hdr,'COMMENT'
        sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .',hdr,/COMMENT
        modfits,rawbc_list[i],0,hdr,EXTEN_NO=0
        FITS_INFO, rawb_list[i],N_ext=numext,/SILENT
        if numext gt 0 then sigma_imbc=mrdfits(rawb_list[i],1,hdr1,/SILENT) $
        else sigma_imbc=sqrt(raw_imbc+r^2)
        if numext gt 1 then dq_imbc=mrdfits(rawb_list[i],2,hdr2,/SILENT) $
        else  dq_imbc=bytarr(nxc,nyc)
        cc_mask=mrdfits(mask_list[i],0,mask_hdr,/SILENT)
        cr_loc=where(cc_mask eq 1)
        dq_imbc[cr_loc]=1
        prb=where(dq_imbc ge 1)
        if total(prb) ne -1 then dq_imbc[prb]=1
        nbad=n_elements(prb)
        dqfactr=float(nbad/n_elements(dq_imb))
        sxaddpar, hdr, 'NBADPIX ',nbad,'Number of bad data quality pixels'
        sxaddpar, hdr, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
          ' fraction of total number of pixels'
        mwrfits,sigma_imbc,rawbc_list[i],hdr
        mwrfits,dq_imbc,rawbc_list[i],hdr
        undefine,prb,nbad,dqfactr
      endfor
      logprint,'CONTROL: Corrections carried out for cosmic rays on spectrum using LA COSMIC,'$
               +' sigma clip used is '+STRTRIM(STRING(crclip),2)+'.'
    endif else begin
      logprint,'CONTROL: Cosmic ray correction not carried out as per user request.'
      rawbc_list=rawb_list
      crflg=0
    endelse
  endif else rawbc_list=rawb_list
  spectrum_cr_cr_end:
  ;dark and flat correction
  if datatype(mdark,2) ne 8 then begin
    logprint,'CONTROL could not find master dark file.'$
             +' Do you want to continue reducing this spectrum assuming dark to be zero'
    logprint,'Press q to exit. Press any key to use dark as zero.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master dark of zeros of size 2048x512'
    mdark={im:dblarr(2048,512),error:dblarr(2048,512),dq:bytarr(2048,512)}
    sigma_d=dblarr(2048,512)
    md_dq=bytarr(2048,512)
  endif
  if datatype(mflat,2) ne 8 then begin
    logprint,'CONTROL could not find master flat file.'$
             +' Do you want to continue reducing this spectrum assuming flat to be one'
    logprint,'Press q to exit. Press any key to use flat as one.'
    Rin = GET_KBRD()
    if Rin eq 'q' then begin
      logprint,'CONTROL: Exiting as requested by the user.'
      goto,code_end
    endif
    logprint,'Creating master flat of ones of size 2048x512'
    sigma_fbdn=make_array(2048,512,value=1.0,/DOUBLE)
    nflat=make_array(2048,512,value=1.0,/DOUBLE)
    mf_dq=bytarr(2048,512)
  endif
  for i=0,n_elements(spectra)-1 do begin
    raw_imbc=mrdfits(rawbc_list[i],0,hdr,/SILENT)
    ;print,rawbc_list[i]
    FITS_INFO, rawbc_list[i],N_ext=numext,/SILENT
    bnx=(size(raw_imbc))[1]
    bny=(size(raw_imbc))[2]

    if numext gt 0 then sigma_imbc=mrdfits(rawbc_list[i],1,hdr1,/SILENT) $
    else sigma_imbc=dblarr(bnx,bny)    ; error after bias correction
    ;if n_elements(sigma_imbc) eq 1 then sigma_imbc=dblarr(bnx,bny)
    if numext gt 1 then dq_imbc=mrdfits(rawbc_list[i],2,hdr2,/SILENT) $
    else dq_imbc=bytarr(bnx,bny)
    if n_elements(dq_imbc) eq 1 then dq_imbc=dblarr(bnx,bny)
    if datatype(r,2) eq 0 then r=12.25
    sigma_imbc=sqrt(raw_imbc+r^2)
    if (SXPAR(hdr, 'DCFLG',MISSING=-1)) ne -1 then begin
      dcfg_chk=SXPAR(hdr, 'DCFLG',MISSING=-1)
      if dcfg_chk eq 1 then begin
        logprint,'CONTROL: The spectrum('+rawbc_list[i]+$
                  ') is already dark corrected. Skipping dark correction for current spectrum.'
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
      ;data quality
      prb=where(dq_imbd ge 1)
      if total(prb) ne -1 then dq_imbd[prb]=1
      nbad=n_elements(prb)
      dqfactr=float(nbad/n_elements(dq_imb))
      ;headers
      sxaddpar, hdr, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      sxaddpar, hdr, 'SIGMDARK',SXPAR(mdark.hdr,'SIGMDARK'), $
                'Standard deviation of the master dark frame'
      sxaddpar, hdr, 'MEANDARK',SXPAR(mdark.hdr,'MEANDARK'), $
                'Mean value of the master dark frame'
      sxaddpar, hdr, 'MDNDARK ',SXPAR(mdark.hdr,'MDNDARK'), $
                'Median value of the master dark frame'
      sxaddpar, hdr, 'MAXDARK', SXPAR(mdark.hdr,'MAXDARK'), $
                'Maximum value of the master dark frame'
      sxaddpar, hdr, 'MINDARK', SXPAR(mdark.hdr,'MINDARK'), $
                'Minimum value of the master bias frame'
      sxaddpar, hdr, 'NFRAMEMD',SXPAR(mdark.hdr,'NFRAMES'), $
                'Number of frames used in dark combine'
      sxaddpar, hdr, 'DARKTYP ',SXPAR(mdark.hdr,'DARKTYP'), $
                'Type of dark combine employed'
      sxaddpar, hdr, 'DARKSAT ',SXPAR(mdark.hdr,'DARKSAT'), $
                'Saturation limit used in dark frames'
      sxaddpar, hdr, 'DARKSIG ',SXPAR(mdark.hdr,'DARKSIG'), $
                'Deviation limit for good dark frames'
      sxaddpar, hdr, 'NBADPIX ',nbad,'Number of bad data quality pixels'
      sxaddpar, hdr, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
                ' fraction of total number of pixels'
      sxaddpar, hdr1, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      sxaddpar, hdr2, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      sxdelpar, hdr,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .',hdr,/COMMENT
      ;hdr_d=update_header(hdr)
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',raw_imbd,hdr
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',sigma_imbd,$
                                             hdr1, /APPEND
      
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcd.fits',dq_imbd,$
                                             hdr2, /APPEND
      logprint,'CONTROL: Dark correction carried out on spectrum('+spectra_name[i]+').'
      rawbd_file=inter_path+spectra_name[i]+'_bcd.fits'
       undefine,prb,nbad,dqfactr
    endif else begin
      rawbd_file=rawbc_list[i]
      raw_imbd=raw_imbc
      sigma_imbd=sigma_imbc
      dq_imbd=dq_imbc
    endelse
    if cr_cor eq 1 then begin
      if (SXPAR(hdr, 'CRFLG',MISSING=-1)) ne 1 then begin
        logprint,'CONTROL: Cosmic Ray correction carried out after dark correction on spectrum('$
                  +spectra_name[i]+') with a clip value of '+STRTRIM(STRING(crclip),2)+'.'
        if ccd_gain eq 0 then ccd_gain=1
        rawbdc_file = inter_path+spectra_name[i]+'_bcdc.fits'
        rawbd_mask  = inter_path+spectra_name[i]+'_bdc_mask.fits'
        la_cosmic,rawbd_file,outlist=rawbdc_file,masklist=rawbd_mask,$
                  gain=ccd_gain,readn=r,sigclip=crclip ; add other parameters after testing
        crflg=1
        sxaddpar, hdr,'CRFLG  ', crflg,'COSMIC RAY CORRECTION FLAG'
        sxaddpar, hdr,'CSRYWEN', 'after dark','When cosmic ray correction was carried out'
        sxaddpar, hdr,'CRCYCLP', crclip,'Clip value used for cosmic ray'
        modfits,rawbdc_file,0,hdr,EXTEN_NO=0
        raw_imbd=mrdfits(rawbdc_file,0,crhdr,/SILENT)
        cc_mask=mrdfits(rawbd_mask,0,mask_hdr,/SILENT)
        cr_loc=where(cc_mask eq 1)
        dq_imbd[cr_loc]=1
        prb=where(dq_imbd ge 1)
        if total(prb) ne -1 then dq_imbd[prb]=1
        nbad=n_elements(prb)
        dqfactr=float(nbad/n_elements(dq_imb))
        sxaddpar, hdr, 'NBADPIX ',nbad,'Number of bad data quality pixels'
        sxaddpar, hdr, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
          ' fraction of total number of pixels'
        
        if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdc.fits',sigma_imbd,$
                                               hdr1, /APPEND
        if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdc.fits',dq_imbd,$
                                               hdr2, /APPEND
      undefine,prb,nbad,dqfactr
      endif
    endif

    spectrum_dr_cr_end:
    if(crf) then begin
      logprint,'CONTROL will flat correct the science spectrum now'
      if (SXPAR(hdr, 'FCFLG',MISSING=-1)) ne -1 then begin
        fcfg_chk=SXPAR(hdr, 'FCFLG',MISSING=-1)
        if fcfg_chk eq 1 then begin
          logprint,'CONTROL: The spectrum('+rawbc_list[i]$
                   +') is already flat corrected. Skipping flat correction for current spectrum.'
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
      ;data quality
      prb=where(dq_imbdf ge 1)
      if total(prb) ne -1 then dq_imbdf[prb]=1
      nbad=n_elements(prb)
      dqfactr=float(nbad/n_elements(dq_imb))
      ;headers
      sxaddpar, hdr, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      sxaddpar, hdr, 'SIGMFLAT',SXPAR(mflat.hdr,'SIGMFLAT'), $
        'Standard deviation of the master flat frame'
      sxaddpar, hdr, 'MEANFLAT',SXPAR(mflat.hdr,'MEANFLAT'), $
        'Mean value of the master flat frame'
      sxaddpar, hdr, 'MDNFLAT ',SXPAR(mflat.hdr,'MDNFLAT'), $
        'Median value of the master dark frame'
      sxaddpar, hdr, 'MAXFLAT', SXPAR(mflat.hdr,'MAXFLAT'), $
        'Maximum value of the master dark frame'
      sxaddpar, hdr, 'MINFLAT', SXPAR(mflat.hdr,'MINFLAT'), $
        'Minimum value of the master bias frame'
      sxaddpar, hdr, 'NFRAMEMF',SXPAR(mflat.hdr,'NFRAMES'), $
        'Number of frames used in dark combine'
      sxaddpar, hdr, 'FLATTYP ',SXPAR(mflat.hdr,'FLATTYP'), $
        'Type of dark combine employed'
      sxaddpar, hdr, 'FLATSAT ',SXPAR(mflat.hdr,'FLATSAT'), $
        'Saturation limit used in dark frames'
      sxaddpar, hdr, 'FLATSIG ',SXPAR(mflat.hdr,'FLATSIG'), $
        'Deviation limit for good dark frames'
      sxaddpar, hdr, 'NBADPIX ',nbad,'Number of bad data quality pixels'
      sxaddpar, hdr, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
        ' fraction of total number of pixels'
      sxaddpar, hdr1, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      sxaddpar, hdr2, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      sxdelpar, hdr,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .',hdr,/COMMENT
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',raw_imbdf,hdr
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',sigma_imbdf,$
                                             hdr1, /APPEND
      
      undefine,prb,nbad,dqfactr
      if save_temp_files eq 1 then writefits,inter_path+spectra_name[i]+'_bcdf.fits',dq_imbdf,hdr2 $
                                              ,/APPEND
      logprint,'CONTROL: Flat correction carried out on spectrum('+rawbc_list[i]+').'
      writefits,out_path+spectra_bname[i]+'_2d.fits',raw_imbdf,hdr
      writefits,out_path+spectra_bname[i]+'_2d.fits',sigma_imbdf,hdr1, /APPEND
      writefits,out_path+spectra_bname[i]+'_2d.fits',dq_imbdf,hdr2, /APPEND
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
        logprint,'CONTROL: The spectra('+spectra_name[i]+$
                 ') requested for extraction is not bias or dark or flat corrected.'
        logprint,'Press q to skip the extraction for current spectrum.'
        logprint,'Press any key to continue extraction for the current spectrum.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the extraction for current spectrum as requested.'
          goto,spectrum_loop_end
        endif
      endif
      if (fix(SXPAR(hdr, 'EXTFLF',MISSING=-1))) eq 1 then begin
        logprint,'CONTROL: The spectra('+rawbc_list[i]+$
                 ') requested for extraction is already an extracted spectrum'
        logprint,'CONTROL: Skipping the extraction for current spectrum as requested.'
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
          errorlog,'CONTROL: CONTROL_TRACE returned without extracted spectrum('+spectra_name[i]+$
                   '). Skipping further pipeline procedures for the corresponding file.'
          goto,spectrum_loop_end
        endif
      endif else begin
        if datatype(spectrum,2) ne 8 then begin
          errorlog,'CONTROL: CONTROL_TRACE returned without extracted spectrum('+spectra_name[i]+$
                   '). Skipping further pipeline procedures for the corresponding file.'
          goto,spectrum_loop_end
        endif
      endelse
      ext_spectrum=spectrum.data
      hdr=spectrum.header
      spectrum_save={data:spectrum.data,error:spectrum.error,dq:spectrum.dq,$
                     background:spectrum.background,bg_error:spectrum.bck_error,$
                     bg_dq:spectrum.dq_bg}
      hdr_nw=update_header(hdr)
      ;update header
      sxdelpar, hdr, 'NBADPIX '
      sxdelpar, hdr, 'DQFACTR '
      sxaddpar, hdr_nw,'TTYPE1  ','data','Extracted 1D spectrum (not background corrected).'
      sxaddpar, hdr_nw,'TFORM1  ','17float','Data format of field 1'
      sxaddpar, hdr_nw,'TUNIT1  ','counts','Physical unit of field 1'
      sxaddpar, hdr_nw,'TDISP1  ','Float','Display format for column 1'
      sxaddpar, hdr_nw,'TNULL1  ','NAN','Undefined value for column 1'
      sxaddpar, hdr_nw,'TTYPE2  ','error','Error of extracted 1D spectrum.'
      sxaddpar, hdr_nw,'TFORM2  ','17float','Data format of field 2'
      sxaddpar, hdr_nw,'TUNIT2  ','counts','Physical unit of field 2'
      sxaddpar, hdr_nw,'TDISP2  ','Float','Display format for column 2'
      sxaddpar, hdr_nw,'TNULL2  ','NAN','Undefined value for column 2'
      sxaddpar, hdr_nw,'TTYPE3  ','dq','Data quality of extracted 1D spectrum.'
      sxaddpar, hdr_nw,'TFORM3  ','BYTE','Data format of field 3'
      sxaddpar, hdr_nw,'TUNIT3  ','BYTE','Physical unit of field 3'
      sxaddpar, hdr_nw,'TDISP3  ','BYTE','Display format for column 3'
      sxaddpar, hdr_nw,'TNULL3  ','0','Undefined value for column 3'
      sxaddpar, hdr_nw,'TTYPE4  ','background','Extracted 1D background spectrum.'
      sxaddpar, hdr_nw,'TFORM4  ','17float','Data format of field 4'
      sxaddpar, hdr_nw,'TUNIT4  ','counts','Physical unit of field 4'
      sxaddpar, hdr_nw,'TDISP4  ','Float','Display format for column 4'
      sxaddpar, hdr_nw,'TNULL4  ','NAN','Undefined value for column 4'
      sxaddpar, hdr_nw,'TTYPE5  ','bg_error','Error of extracted 1D background spectrum.'
      sxaddpar, hdr_nw,'TFORM5  ','17float','Data format of field 5'
      sxaddpar, hdr_nw,'TUNIT5  ','counts','Physical unit of field 5'
      sxaddpar, hdr_nw,'TDISP5  ','Float','Display format for column 5'
      sxaddpar, hdr_nw,'TNULL5  ','NAN','Undefined value for column 5'
      sxaddpar, hdr_nw,'TTYPE6  ','bg_dq','Data quality of extracted 1D background spectrum.'
      sxaddpar, hdr_nw,'TFORM6  ','BYTE','Data format of field 6'
      sxaddpar, hdr_nw,'TUNIT6  ','BYTE','Physical unit of field 6'
      sxaddpar, hdr_nw,'TDISP6  ','BYTE','Display format for column 6'
      sxaddpar, hdr_nw,'TNULL6  ','0','Undefined value for column 6'
      sxdelpar, hdr_nw,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .',hdr_nw,/COMMENT
      
      if save_temp_files eq 1 then mwrfits,spectrum_save,inter_path+spectra_name[i]+'_bcdfe.fits',$
                                           hdr_nw,/create
      logprint,'CONTROL: Spectrum('+spectra_name[i]+') is extracted.'
    endif ;else begin
    extr_read:
    ;im_bflg=fix(SXPAR(hdr, 'BCFLG',MISSING=-1))
    ;im_dflg=fix(SXPAR(hdr, 'DCFLG',MISSING=-1))
    ;im_fflg=fix(SXPAR(hdr, 'FCFLG',MISSING=-1))
    ;if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1) then begin
    ;  logprint,'CONTROL: The spectra is not bias or dark or flat corrected with CONTROL.'
    ;        logprint,'Press q to skip the extraction for current spectrum. 
    ;        logprint,'Press any key to continue assuming the spectrum is corrected.'
    ;        R = GET_KBRD()
    ;        if R eq 'q' then begin
    ;          logprint,'CONTROL: Skipping the extraction for current spectrum as '$
    ;                  +'requested by the user.'
    ;          goto,spectrum_loop_end
    ;        endif
    ;      endif
    ;      spectrum=mrdfits(rawbc_list[i],1,hdr_nw)
    ;
    ;      if(tag_exist(spectrum,'data') eq 0) then begin
    ;        logprint,'CONTROL: Could not find data in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif
    ;      if(tag_exist(spectrum,'error') eq 0) then begin
    ;        logprint,'CONTROL: Could not find error in FITS file specified.'$
    ;                +' Assuming it to be zero'
    ;        sp_error=dblarr(n_elements(spectrum.data))
    ;        struct_add_field, spectrum, 'error', sp_error, after='data'
    ;      endif
    ;      if(tag_exist(spectrum,'dq') eq 0) then begin
    ;        logprint,'CONTROL: Could not find Data quality in FITS file specified.'$
    ;                +' Assuming data to be good'
    ;        sp_dq=bytarr(n_elements(spectrum.data))
    ;        struct_add_field, spectrum, 'dq', sp_dq, after='error'
    ;      endif
    ;      if(tag_exist(spectrum,'background') eq 0) then begin
    ;        logprint,'CONTROL: Could not find background in FITS file specified.'$
    ;                +' Assuming it to be zero'
    ;        sp_background=dblarr(n_elements(spectrum.data))
    ;        struct_add_field, spectrum, 'background', sp_background, after='dq'
    ;      endif
    ;      if(tag_exist(spectrum,'bck_error') eq 0) then begin
    ;        logprint,'CONTROL: Could not find background error in FITS file specified.'$
    ;                +' Assuming it to be zero'
    ;        sp_backgroundr=dblarr(n_elements(spectrum.data))
    ;        struct_add_field, spectrum, 'bck_error', sp_backgroundr, after='background'
    ;      endif
    ;      if(tag_exist(spectrum,'dq_bg') eq 0) then begin
    ;        logprint,'CONTROL: Could not find Data qualiy for background in FITS file specified.'$
    ;                +' Assuming it be all good'
    ;        sp_back_dq=dblarr(n_elements(spectrum.data))
    ;        struct_add_field, spectrum, 'dq_bg', sp_back_dq, after='bck_error'
    ;      endif
    ;    endelse
    if (bgs) then begin
      logprint,'CONTROL will carryout background substraction on extracted spectrum now'
      if (idl_ver ge 8) then begin
        if isa(spectrum,'STRUCT') eq 0 then begin
          logprint,'Seems like extracted spectra were not avaliable.'$
                   +' Assuming the extracted spectra to be in the data folder files.'
          spectrum=mrdfits(rawbc_list[i],1,hdr_nw,/SILENT)
        endif
      endif else begin
        if datatype(spectrum,2) ne 8 then begin
          logprint,'Seems like extracted spectra were not avaliable.'$
                   +' Assuming the extracted spectra to be in the data folder files.'
          spectrum=mrdfits(rawbc_list[i],1,hdr_nw,/SILENT)
        endif
      endelse
      im_bflg=fix(SXPAR(hdr_nw, 'BCFLG',MISSING=-1))
      im_dflg=fix(SXPAR(hdr_nw, 'DCFLG',MISSING=-1))
      im_fflg=fix(SXPAR(hdr_nw, 'FCFLG',MISSING=-1))
      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1) then begin
        logprint,'CONTROL: The spectra('+spectra_name[i]$
                 +') is not bias or dark or flat corrected with CONTROL.'
        logprint,'Press q to skip the background subtraction for current spectrum.'$
                 +' Press any key to continue assuming the spectrum is corrected.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the extraction for current spectrum as requested.'
          goto,spectrum_loop_end
        endif
      endif
      if(tag_exist(spectrum,'data') eq 0) then begin
        errorlog,'CONTROL: Could not find data in FITS file specified.'$
                 +' Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      if(tag_exist(spectrum,'error') eq 0) then begin
        errorlog,'CONTROL: Could not find error in FITS file specified.'$
                 +' Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      if(tag_exist(spectrum,'dq') eq 0) then begin
        logprint,'CONTROL: Could not find Data quality in FITS file specified.'$
                 +' Assuming data to be good'
        errorlog,'CONTROL: Could not find Data quality in FITS file specified.'$
                 +' Assuming data to be good',logonly=1
        sp_dq=bytarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'dq', sp_dq, after='error'
      endif
      if(tag_exist(spectrum,'background') eq 0) then begin
        logprint,'CONTROL: Could not find background in FITS file specified.'$
                 +' Assuming it to be zero'
        errorlog,'CONTROL: Could not find background in FITS file specified.'$
                 +' Assuming it to be zero',logonly=1
        sp_background=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'background', sp_background, after='dq'
      endif
      if(tag_exist(spectrum,'bck_error') eq 0) then begin
        logprint,'CONTROL: Could not find background error in FITS file specified.'$
                 +' Assuming it to be zero'
        errorlog,'CONTROL: Could not find background error in FITS file specified.'$
                 +' Assuming it to be zero',logonly=1
        sp_backgroundr=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'bck_error', sp_backgroundr, after='background'
      endif
      if(tag_exist(spectrum,'dq_bg') eq 0) then begin
        logprint,'CONTROL: Could not find Data qualiy for background in FITS file specified.'$
                 +' Assuming it be all good'
        errorlog,'CONTROL: Could not find Data qualiy for background in FITS file specified.'$
                 +' Assuming it be all good',logonly=1
        sp_back_dq=dblarr(n_elements(spectrum.data))
        struct_add_field, spectrum, 'dq_bg', sp_back_dq, after='bck_error'
      endif
      if (size(spectrum.data))[0] ne (size(spectrum.background))[0] then begin
        errorlog,'CONTROL: Dimension missmatch for data and background in spectrum file.'$
                 +' Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      spectra_1d=spectrum.data-spectrum.background
      if (((size(spectrum.error))[0]) ne ((size(spectrum.bck_error))[0])) then begin
        errorlog,'CONTROL: Dimension missmatch for data error and background error in spectrum.'$
                 +' Pipeline procedures for this file is aborted'
        goto,spectrum_loop_end
      endif
      error_1d=sqrt(spectrum.error^2+spectrum.bck_error^2)
      bgflg=1
      ;data quality
      dq_1d=spectrum.dq+spectrum.dq_bg
      prb=where(dq_1d ge 1)
      if total(prb) ne -1 then dq_1d(prb)=1
      nbad=n_elements(prb)
      dqfactr=float(nbad/n_elements(dq_imb))
      if (size(dq_1d))[0] ne 1 then begin
        logprint,'CONTROL: Data quality array is not one diamensional.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif
      logprint,'CONTROL: Spectrum('+spectra_name[i]+') is now background subtracted.'
      sxaddpar, hdr_nw, 'BGFLG', bgflg ,'BACKGROUND CORRECTION FLAG'
      hdr_1d=update_header(hdr_nw)
      ;update header
      sxaddpar, hdr_1d,'TTYPE1  ','data','Extracted 1D spectrum.'
      sxaddpar, hdr_1d,'TFORM1  ','17float','Data format of field 1'
      sxaddpar, hdr_1d,'TUNIT1  ','counts','Physical unit of field 1'
      sxaddpar, hdr_1d,'TDISP1  ','float','Display format for column 1'
      sxaddpar, hdr_1d,'TNULL1  ','NAN','Undefined value for column 1'
      sxaddpar, hdr_1d,'TTYPE2  ','error','Error of extracted 1D spectrum.'
      sxaddpar, hdr_1d,'TFORM2  ','17float','Data format of field 2'
      sxaddpar, hdr_1d,'TUNIT2  ','counts','Physical unit of field 2'
      sxaddpar, hdr_1d,'TDISP2  ','float','Display format for column 2'
      sxaddpar, hdr_1d,'TNULL2  ','NAN','Undefined value for column 2'
      sxaddpar, hdr_1d,'TTYPE3  ','dq','Data quality of extracted 1D spectrum.'
      sxaddpar, hdr_1d,'TFORM3  ','byte','Data format of field 3'
      sxaddpar, hdr_1d,'TUNIT3  ','bytes','Physical unit of field 3'
      sxaddpar, hdr_1d,'TDISP3  ','bytes','Display format for column 3'
      sxaddpar, hdr_1d,'TNULL3  ','0','Undefined value for column 3'
      sxaddpar, hdr_1d, 'NBADPIX ',nbad,'Number of bad data quality pixels'
      sxaddpar, hdr_1d, 'DQFACTR ',dqfactr,'Number of bad data quality pixels as a'+ $
        ' fraction of total number of pixels'
      sxdelpar, hdr_1d,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .',hdr_1d,/COMMENT
        
      spectrum_file={data:spectra_1d,error:error_1d,dq:dq_1d}
      if save_temp_files eq 1 then mwrfits,spectrum_file,out_path+spectra_bname[i]+'_1d.fits',$
                                           hdr_1d,/create
      undefine,prb,nbad,dqfactr
    endif ;else begin
    ;      im_bflg=fix(SXPAR(hdr_nw, 'BCFLG',MISSING=-1))
    ;      im_dflg=fix(SXPAR(hdr_nw, 'DCFLG',MISSING=-1))
    ;      im_fflg=fix(SXPAR(hdr_nw, 'FCFLG',MISSING=-1))
    ;      im_eflg=fix(SXPAR(hdr_nw, 'EXTFLF',MISSING=-1))
    ;      if (im_bflg eq -1 or im_dflg eq -1 or im_fflg eq -1 or im_eflg eq -1) then begin
    ;        logprint,'CONTROL: The spectra('+rawbc_list[i]+') is not bias or dark or flat'$
    ;                +'corrected or extracted with CONTROL.'
    ;        logprint,'Press q to skip the background subtraction for current spectrum ('$
    ;        +rawbc_list[i]+'). Press any key to continue background correction for the spectrum.'
    ;        R = GET_KBRD()
    ;        if R eq 'q' then begin
    ;          logprint,'CONTROL: Skipping the background subtraction for current spectrum as'$
    ;                  +'requested by the user.'
    ;          goto,spectrum_loop_end
    ;        endif
    ;      endif
    ;      if(tag_exist(spectrum,'data') eq 0) then begin
    ;        logprint,'CONTROL: Could not find data in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else spectra_1d=spectrum.data
    ;      if(tag_exist(spectrum,'error') eq 0) then begin
    ;        logprint,'CONTROL: Could not find error in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else error_1d=spectrum.error
    ;      if(tag_exist(spectrum,'dq') eq 0) then begin
    ;        logprint,'CONTROL: Could not find Data quality in FITS file specified.'$
    ;                +' Assuming data to be good'
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
    ;wavelength calibration
    if (wcl) then begin
      logprint,'CONTROL will wavelength calibrate the 1D science spectrum now'
      if (idl_ver ge 8) then begin
        if isa(spectrum_file,'STRUCT') eq 0 then $
           spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endif else begin
        if datatype(spectrum,2) ne 8 then spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endelse
      im_eflg=fix(SXPAR(hdr_1d, 'EXTFLF',MISSING=-1))
      im_bgflg=fix(SXPAR(hdr_1d, 'BGFLG',MISSING=-1))
      if (im_eflg eq -1 or im_bgflg eq -1) then begin
        logprint,'CONTROL: The specified FITS file is not extracted'$
                 +' or background subtracted with CONTROL.'
        logprint,'Press q to skip the wavelength calibration for current spectrum ('$
                 +rawbc_list[i]+'). Press any key to continue with wavelength calibration.'
        Rin = GET_KBRD()
        if Rin eq 'q' then begin
          logprint,'CONTROL: Skipping the wavelength calibration for current spectrum'$
                   +' as requested by the user.'
          goto,spectrum_loop_end
        endif
        goto,spectrum_loop_end
      endif
      if(tag_exist(spectrum_file,'data') eq 0) then begin
        logprint,'CONTROL: Could not find data in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else spectra_1d=spectrum.data
      if(tag_exist(spectrum_file,'error') eq 0) then begin
        logprint,'CONTROL: Could not find error in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else error_1d=spectrum.error
      if(tag_exist(spectrum_file,'dq') eq 0) then begin
        logprint,'CONTROL: Could not find Data quality in FITS file specified.'$
                 +' Assuming data to be good'
        sp_dq=bytarr(n_elements(spectrum_file.data))
        struct_add_field, spectrum_file, 'dq', sp_dq, after='error'
      endif else dq_1d=spectrum_file.dq
      ;print,spectra_name[i]
      if(tag_exist(infile,'wavecal_mode') eq 1) then wavecal_type=infile.wavecal_mode else $
                                                     wavecal_type='simple'
      if(tag_exist(spectrum_file,'data') eq 0) then begin
        logprint,'CONTROL: Could not find data in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif
      if (size(spectrum_file.data))[0] ne 1 then begin
        logprint,'CONTROL: FITS file specified ('+spectra_name[i]+') do not have a 1D spectrum.'$
                 +' Pipeline procedures for this file will be aborted'
        goto,spectrum_loop_end
      endif
      wcl_spectrum=control_wavecal(spectrum_file.data,hdr_1d,infile,wavecal_type)
      sp_wavelength=wcl_spectrum.wavelength
      sp_data=wcl_spectrum.flux
      sp_error=spectrum_file.error
      sp_dq=spectrum_file.dq
      wl_hdr=update_header(wcl_spectrum.header)
      wcl_spectrum_file={wave:sp_wavelength,flux:sp_data,error:sp_error,dq:sp_dq}
      window,xsize=1300,ysize=800
      cgplot,sp_wavelength,sp_data,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,$
             xtitle='wavelength [$\Angstrom$]',ytitle='flux [counts]'
      write_png,out_path+spectra_bname[i]+'_1dw.png',TVRD(/TRUE)
      ;defininig header
      sxaddpar, wl_hdr,'TTYPE1  ','wave','Wavelength of the spectrum. '
      sxaddpar, wl_hdr,'TFORM1  ','float','Data format of field 1'
      sxaddpar, wl_hdr,'TUNIT1  ','Angstroms','Physical unit of field 1'
      sxaddpar, wl_hdr,'TDISP1  ','float','Display format for column  3'
      sxaddpar, wl_hdr,'TNULL1  ','NAN','Undefined value for column  1'
      sxaddpar, wl_hdr,'TTYPE2  ','flux','Wavelength claibrated 1D flux in counts.'
      sxaddpar, wl_hdr,'TFORM2  ','float','Data format of field 2'
      sxaddpar, wl_hdr,'TUNIT2  ','counts','Physical unit of field 2'
      sxaddpar, wl_hdr,'TDISP2  ','float','Display format for column  2'
      sxaddpar, wl_hdr,'TNULL2  ','NAN','Undefined value for column  2'
      sxaddpar, wl_hdr,'TTYPE3  ','error','Error of wavelength claibrated 1D flux in counts.'
      sxaddpar, wl_hdr,'TFORM3  ','float','Data format of field 3'
      sxaddpar, wl_hdr,'TUNIT3  ','counts','Physical unit of field 3'
      sxaddpar, wl_hdr,'TDISP3  ','float','Display format for column  3'
      sxaddpar, wl_hdr,'TNULL3  ','NAN','Undefined value for column  3'
      sxaddpar, wl_hdr,'TTYPE4  ','dq','Data quality of wavelength calibrated 1D spectrum. '
      sxaddpar, wl_hdr,'TFORM4  ','bytes','Data format of field 1'
      sxaddpar, wl_hdr,'TUNIT4  ','bytes','Physical unit of field 1'
      sxaddpar, wl_hdr,'TDISP4  ','bytes','Display format for column  3'
      sxaddpar, wl_hdr,'TNULL4  ','0','Undefined value for column  1'
      sxdelpar, wl_hdr,'COMMENT'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0.',wl_hdr,/COMMENT

      mwrfits,wcl_spectrum_file,out_path+spectra_bname[i]+'_1dw.fits',wl_hdr, /CREATE
      spectra_1dw  = sp_data
      error_1dw    = sp_error
      dq_1dw       = sp_dq
    endif ;else begin
    ;      spectrum_wl=mrdfits(rawbc_list[i],1,wl_hdr)
    ;      if(tag_exist(spectrum_wl,'wave') eq 0) then begin
    ;        logprint,'CONTROL: Could not find wavelength information in FITS file specified.'$
    ;                +'Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else sp_wavelength=spectrum_wl.wave
    ;      if(tag_exist(spectrum_wl,'counts') eq 0) then begin
    ;        logprint,'CONTROL: Could not find count information in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else spectra_1dw=spectrum_wl.counts
    ;      if(tag_exist(spectrum_wl,'error') eq 0) then begin
    ;        logprint,'CONTROL: Could not find error information in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else error_1dw=spectrum_wl.error
    ;      if(tag_exist(spectrum_wl,'dq') eq 0) then begin
    ;        logprint,'CONTROL: Could not find Data quality information in FITS file specified.'$
    ;                +' Pipeline procedures for this file is aborted'
    ;        goto,spectrum_loop_end
    ;      endif else dq_1dw=spectrum_wl.dq
    ;
    ;    endelse
    ;flux calibration
    if (fcl) then begin
      logprint,'CONTROL will flux calibrate the wavelength calibrated 1D science spectrum now'
      if (idl_ver ge 8) then begin
        if isa(wcl_spectrum_file,'STRUCT') eq 0 then $
               wcl_spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endif else begin
        if datatype(spectrum,2) ne 8 then wcl_spectrum_file=mrdfits(rawbc_list[i],1,hdr_1d,/SILENT)
      endelse
      if  fix(SXPAR(wl_hdr, 'WCALFLG',MISSING=-1)) eq -1 then begin
        logprint,'CONTROL: The file specified('+rawbc_list[i]+') is not wavelength calibrated.'$
                 +' Pipeline procedures for this file will be aborted'
        goto,spectrum_loop_end
      endif
      if tag_exist(infile,'flux_cal') eq 0 then begin
        logprint,'CONTROL:  Flux calibration file not found.'$
                 +' Please re-run the pipeline with the file.',logonly=logonly
        message,'CONTROL: Flux calibration file not found.'$
                 +' Please re-run the pipeline with actual file.'
      endif
      flux_calib_file=detectos(infile.flux_cal)
      if (file_test(flux_calib_file) eq 0) then begin
        new_flux_calib_file=detectos('calibration/flux_cal.txt')
        if (file_test(new_flux_calib_file) eq 0) then begin
          logprint,'CONTROL: Flux calibration file not found.'$
                   +' Please re-run the pipeline with actual file address',logonly=logonly
          message,'CONTROL: Flux calibration file not found.'$
                   +' Please re-run the pipeline with actual file address'
        endif else begin
          logprint,'CONTROL: Flux calibration file found in calibration folder.'
          flux_calib_file =detectos(new_flux_calib_file)
        endelse
      endif
      readcol,flux_calib_file,fcal_wave,fcal_value,F='D,D' ;flux calibration file file location
      if(tag_exist(wcl_spectrum_file,'wave') eq 0) then begin
        logprint,'CONTROL: Could not find wavelength information in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else sp_wavelength=wcl_spectrum_file.wave
      if(tag_exist(wcl_spectrum_file,'flux') eq 0) then begin
        logprint,'CONTROL: Could not find count information in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else spectra_1dw=wcl_spectrum_file.flux
      if(tag_exist(wcl_spectrum_file,'error') eq 0) then begin
        logprint,'CONTROL: Could not find error information in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else error_1dw=wcl_spectrum_file.error
      if(tag_exist(wcl_spectrum_file,'dq') eq 0) then begin
        logprint,'CONTROL: Could not find Data quality information in FITS file specified.'$
                 +' Pipeline procedures for this file are aborted'
        goto,spectrum_loop_end
      endif else dq_1dw=wcl_spectrum_file.dq
      QE=interpol(fcal_value,fcal_wave,sp_wavelength,/SPLINE)
      exptime=double(sxpar(wl_hdr,'EXPTIME'))
      spectra_1dwf=spectra_1dw*QE/exptime
      error_1dwf=error_1dw*QE/exptime
      wclf_spectrum_file={wave:sp_wavelength,flux:spectra_1dwf,error:error_1dwf,dq:dq_1dw}
      window,xsize=1300,ysize=800
      cgplot,sp_wavelength,spectra_1dwf,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,$
             xtitle='wavelength [$\Angstrom$]',ytitle='flux [ergs s!u-2!n cm!u-2!n $\Angstrom$!u-1!n]'
        write_png,out_path+spectra_bname[i]+'_1dwf.png',TVRD(/TRUE)
      spectrum_hdr=update_header(wl_hdr)
      sxaddpar, spectrum_hdr, 'FCALFLG', 1 ,'Flux correction flag'
      sxaddpar, spectrum_hdr, 'FCALFLE', flux_calib_file,'Location of flux response file'
      ;defininig header
      sxaddpar, spectrum_hdr,'TTYPE1  ','wave','label for field   1'
      sxaddpar, spectrum_hdr,'TFORM1  ','float','Data format of field 1'
      sxaddpar, spectrum_hdr,'TUNIT1  ','Angstroms','Physical unit of field 1'
      sxaddpar, spectrum_hdr,'TDISP1  ','float','Display format for column  3'
      sxaddpar, spectrum_hdr,'TNULL1  ','NAN','Undefined value for column  1'
      sxaddpar, spectrum_hdr,'TTYPE2  ','flux','label for field   2'
      sxaddpar, spectrum_hdr,'TFORM2  ','float','Data format of field 2'
      sxaddpar, spectrum_hdr,'TUNIT2  ','erg/cm**2/Angstrom','Physical unit of field 2'
      sxaddpar, spectrum_hdr,'TDISP2  ','float','Display format for column  2'
      sxaddpar, spectrum_hdr,'TNULL2  ','NAN','Undefined value for column  2'
      sxaddpar, spectrum_hdr,'TTYPE3  ','error','label for field   3'
      sxaddpar, spectrum_hdr,'TFORM3  ','float','Data format of field 3'
      sxaddpar, spectrum_hdr,'TUNIT3  ','erg/cm**2/Angstrom','Physical unit of field 3'
      sxaddpar, spectrum_hdr,'TDISP3  ','float','Display format for column  3'
      sxaddpar, spectrum_hdr,'TNULL3  ','NAN','Undefined value for column  3'
      sxaddpar, spectrum_hdr,'TTYPE1  ','dq','label for field   1'
      sxaddpar, spectrum_hdr,'TFORM1  ','byte','Data format of field 1'
      sxaddpar, spectrum_hdr,'TUNIT1  ','byte','Physical unit of field 1'
      sxaddpar, spectrum_hdr,'TDISP1  ','byte','Display format for column  3'
      sxaddpar, spectrum_hdr,'TNULL1  ','0','Undefined value for column  1'
      sxdelpar, spectrum_hdr,'COMMENT'
        sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .'$
                 ,spectrum_hdr,/COMMENT
      mwrfits,wclf_spectrum_file,out_path+spectra_bname[i]+'_1dwf.fits',spectrum_hdr, /CREATE
    endif

    spectrum_loop_end:
  endfor

  spectrums=file_search(out_path+'*_1dwf.fits')
  if (science ne 0) then begin
      control_light_curve,spectrums,lightcurve,wave_region=['MgI','MgII','FeII']
    
      sxaddpar, lc_hdr, 'TELESCOP', SXPAR(hdr,'TELESCOP'),'Telescope name'
      sxaddpar, lc_hdr, 'ROOTNAME', SXPAR(hdr,'ROOTNAME'),'Root directory'
      sxaddpar, lc_hdr, 'PRGRM_ID', SXPAR(hdr,'PRGRM_ID'),'Program ID'
      sxaddpar, lc_hdr, 'TARGT_ID', SXPAR(hdr,'TARGT_ID'),'Target ID'
      sxaddpar, lc_hdr, 'EXP_ID  ', SXPAR(hdr,'EXP_ID'),'Exposure ID'
      sxaddpar, lc_hdr, 'OBS_ID  ', SXPAR(hdr,'OBS_ID'),'Observation ID'
      sxaddpar, lc_hdr, 'TTYPE1  ','TIME','label for field   1'
      sxaddpar, lc_hdr, 'TFORM1  ','float','Data format of field 2'
      sxaddpar, lc_hdr, 'TUNIT1  ','Julian Date','Physical unit of field 1'
      sxaddpar, lc_hdr, 'TDISP1  ','float','Display format for column  1'
      sxaddpar, lc_hdr, 'TNULL1  ','NAN','Undefined value for column  1'
      sxaddpar, lc_hdr, 'TTYPE2  ','FULL_DATA','label for field   2'
      sxaddpar, lc_hdr, 'TFORM2  ','float','Data format of field 2'
      sxaddpar, lc_hdr, 'TUNIT2  ','erg/s/cm**2/Angstrom','Physical unit of field 2'
      sxaddpar, lc_hdr, 'TDISP2  ','float','Display format for column  2'
      sxaddpar, lc_hdr, 'TNULL2  ','NAN','Undefined value for column  2'
      sxaddpar, lc_hdr, 'TTYPE2  ','FULL_ERROR','label for field   3'
      sxaddpar, lc_hdr, 'TFORM3  ','float','Data format of field 3'
      sxaddpar, lc_hdr, 'TUNIT3  ','erg/s/cm**2/Angstrom','Physical unit of field 3'
      sxaddpar, lc_hdr, 'TDISP3  ','float','Display format for column  3'
      sxaddpar, lc_hdr, 'TNULL3  ','NAN','Undefined value for column  3'
      sxaddpar, lc_hdr, 'TTYPE4  ','SHORT_DATA','label for field   4'
      sxaddpar, lc_hdr, 'TFORM4  ','float','Data format of field 4'
      sxaddpar, lc_hdr, 'TUNIT4  ','erg/s/cm**2/Angstrom','Physical unit of field 4'
      sxaddpar, lc_hdr, 'TDISP4  ','float','Display format for column  4'
      sxaddpar, lc_hdr, 'TNULL4  ','NAN','Undefined value for column  4'
      sxaddpar, lc_hdr, 'TTYPE5  ','SHORT_ERROR','label for field   5'
      sxaddpar, lc_hdr, 'TFORM5  ','float','Data format of field 5'
      sxaddpar, lc_hdr, 'TUNIT5  ','erg/s/cm**2/Angstrom','Physical unit of field 5'
      sxaddpar, lc_hdr, 'TDISP5  ','float','Display format for column  5'
      sxaddpar, lc_hdr, 'TNULL5  ','NAN','Undefined value for column  5'
      sxaddpar, lc_hdr, 'TTYPE6  ','MIDDLE_DATA','label for field   6'
      sxaddpar, lc_hdr, 'TFORM6  ','float','Data format of field 6'
      sxaddpar, lc_hdr, 'TUNIT6  ','erg/s/cm**2/Angstrom','Physical unit of field 6'
      sxaddpar, lc_hdr, 'TDISP6  ','float','Display format for column  6'
      sxaddpar, lc_hdr, 'TNULL6  ','NAN','Undefined value for column  6'
      sxaddpar, lc_hdr, 'TTYPE7  ','MIDDLE_ERROR','label for field   7'
      sxaddpar, lc_hdr, 'TFORM7  ','float','Data format of field 7'
      sxaddpar, lc_hdr, 'TUNIT7  ','erg/s/cm**2/Angstrom','Physical unit of field 7'
      sxaddpar, lc_hdr, 'TDISP7  ','float','Display format for column  7'
      sxaddpar, lc_hdr, 'TNULL7  ','NAN','Undefined value for column  7'
      sxaddpar, lc_hdr, 'TTYPE8  ','LONG_DATA','label for field   8'
      sxaddpar, lc_hdr, 'TFORM8  ','float','Data format of field 8'
      sxaddpar, lc_hdr, 'TUNIT8  ','erg/s/cm**2/Angstrom','Physical unit of field 8'
      sxaddpar, lc_hdr, 'TDISP8  ','float','Display format for column  8'
      sxaddpar, lc_hdr, 'TNULL8  ','NAN','Undefined value for column  8'
      sxaddpar, lc_hdr, 'TTYPE9  ','LONG_ERROR','label for field   9'
      sxaddpar, lc_hdr, 'TFORM9  ','float','Data format of field 9'
      sxaddpar, lc_hdr, 'TUNIT9  ','erg/s/cm**2/Angstrom','Physical unit of field 9'
      sxaddpar, lc_hdr, 'TDISP9  ','float','Display format for column  9'
      sxaddpar, lc_hdr, 'TNULL9  ','NAN','Undefined value for column  9'
      sxaddpar, lc_hdr, 'TTYPE10  ','MGII_DATA','label for field   10'
      sxaddpar, lc_hdr, 'TFORM10  ','float','Data format of field 10'
      sxaddpar, lc_hdr, 'TUNIT10  ','erg/s/cm**2/Angstrom','Physical unit of field 10'
      sxaddpar, lc_hdr, 'TDISP10  ','float','Display format for column  10'
      sxaddpar, lc_hdr, 'TNULL10  ','NAN','Undefined value for column  10'
      sxaddpar, lc_hdr, 'TTYPE11  ','MGII_ERROR','label for field   11'
      sxaddpar, lc_hdr, 'TFORM11  ','float','Data format of field 11'
      sxaddpar, lc_hdr, 'TUNIT11  ','erg/s/cm**2/Angstrom','Physical unit of field 11'
      sxaddpar, lc_hdr, 'TDISP11  ','float','Display format for column  11'
      sxaddpar, lc_hdr, 'TNULL11  ','NAN','Undefined value for column  11'
      sxaddpar, lc_hdr, 'TTYPE12  ','MGI_DATA','label for field   12'
      sxaddpar, lc_hdr, 'TFORM12  ','float','Data format of field 12'
      sxaddpar, lc_hdr, 'TUNIT12  ','erg/s/cm**2/Angstrom','Physical unit of field 12'
      sxaddpar, lc_hdr, 'TDISP12  ','float','Display format for column  12'
      sxaddpar, lc_hdr, 'TNULL12  ','NAN','Undefined value for column  12'
      sxaddpar, lc_hdr, 'TTYPE13  ','MGI_ERROR','label for field   13'
      sxaddpar, lc_hdr, 'TFORM13  ','float','Data format of field 13'
      sxaddpar, lc_hdr, 'TUNIT13  ','erg/s/cm**2/Angstrom','Physical unit of field 13'
      sxaddpar, lc_hdr, 'TDISP13  ','float','Display format for column  13'
      sxaddpar, lc_hdr, 'TNULL13  ','NAN','Undefined value for column  13'
      sxaddpar, lc_hdr, 'TTYPE14  ','FEII_DATA','label for field   14'
      sxaddpar, lc_hdr, 'TFORM14  ','float','Data format of field 14'
      sxaddpar, lc_hdr, 'TUNIT14  ','erg/s/cm**2/Angstrom','Physical unit of field 14'
      sxaddpar, lc_hdr, 'TDISP14  ','float','Display format for column  14'
      sxaddpar, lc_hdr, 'TNULL14  ','NAN','Undefined value for column  14'
      sxaddpar, lc_hdr, 'TTYPE15  ','FEII_ERROR','label for field   15'
      sxaddpar, lc_hdr, 'TFORM15  ','float','Data format of field 15'
      sxaddpar, lc_hdr, 'TUNIT15  ','erg/s/cm**2/Angstrom','Physical unit of field 15'
      sxaddpar, lc_hdr, 'TDISP15  ','float','Display format for column  15'
      sxaddpar, lc_hdr, 'TNULL15  ','NAN','Undefined value for column  15'
      sxaddhist,'File processed with CUTE AUTONOMOUS DATA REDUCTION PIPELINE V1.0 .'$
                ,lc_hdr,/COMMENT
              mwrfits,lightcurve,out_path+'light_curve.fits',lc_hdr, /CREATE

    ;light curve plot:add here
    ;!P.Multi=[0,1,2]
    ;cgplot,tlc1,nlc1,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc1,$
    ;ERR_YHigh=elc1,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Short wavelengths'
    ;
    ;cgplot,tlc2,nlc2,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc2,$
    ;ERR_YHigh=elc2,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Middle wavelengths'
    ;
    ;cgplot,tlc3,nlc3,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc3,$
    ;ERR_YHigh=elc3,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Long wavelengths'
    ;write_png,out_path+'light_curve.png',TVRD(/TRUE)
    ;!P.Multi=0
  endif
  code_end:
  logprint,'CONTROL wrapping up the data reduction procedure.'
  logprint,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',close=close
  close,/all
end