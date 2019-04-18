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
  print,'    create_mflat: Creatze master flat.'
  print,'    cr_bias: Correct for bias in science frames.'
  print,'    cr_dark:Correct for dark in science frames.'
  print,'    cr_flat:Correct for flat in science frames.'
  print,'    cr_cosmic: correct for cosmic rays in science frames.'
  print,'    extract: Define the trace and extract the spectrum.'
  print,'    bg_sub: subtract background from the spectrum.'
  print,'    wcalib: Do wavelength calibration.'
  print,'    fluxcalib: Do Flux calibration.'
  print,'    light_curve: Create three light curves form the data.'
  print,'    retrieval: Process the flux and wavelength calibrated spectra to obtain transmission spectra.'
  print,'    level3: Carry out all process required for generating Level 3 data of CUTE, ie.,  flux and wavelength calibrated 1-D spectra.'
  print,'    all: Execute all the steps. This is also the default option if steps keyword is absent.'
  print,'03. Location of hot and bad pixel map. The location of map of cosmetic defects on the CCD as a FITS file with the same sizes as the input CCD frame'
  print,'   byte-type pixels (8-bit) and values of 1B for good pixels and 0B for bad once.'
  print,'04. Location of master calibration files if available.'
  print,'   This include locations of master bias, master dark and/or master flat files (as fits files).'
  print,'05. Set save_temp_files equal to 1 if you need all the intermediate files to be saved to the intermediate file directory.'
  print,'06. Options for master file creation.'
  print,'   Statistical method to create master files. Options are: median, mean or mode'
  print,'   Threshold for rejection of pixels: sigma_deviation(how many sigmas from mean value of the frame.'
  print,'   Saturation_level sets the saturation limit of CCD. if not set default of 72000 is used.'
  print,'07. Cosmic Ray correction options: Level of cr clipping for LA cosmic algorithm. Default value set the code is 800'
  print,'08. Trace parameters'
  print,'   Degree for centroid polynomial: Default is 1.'
  print,'   Trace type used for extraction. The options are: simple,fixed,variable,function. for details of these options refer to the manual.'
  print,'   Additional parameters required for different trace options:'
  print,'     centroid  : Centroid will be calculated if not provided'
  print,'     slope   : Required for option: simple and fixed, which defines the slope of the spectrum'
  print,'     width     : Required for option: simple, which defines the width of the spectrum'
  print,'     upper     : Required for option: fixed, which defines the upper width of the spectrum'
  print,'     lower     : Required for option: fixed, which defines the lower width of the spectrum'
  print,'     threshold : Required for option: variable or function, which defines the threshold from the maximum/peak of the spectrum.'
  print,'   File location for information of where to extract background. The user could also provide it as fixed number which corresponds to the shift in pixels from centroid.'
  print,'09. Wavelength calibration variables'
  print,'   Type of wavelength calibration required: Options are simple and crscor'
  print,'   Location of wavelength file. A text file with length equal to number of pixels in cross dispersion direction representing wavelength to pixel mapping.'
  print,'   A location of look up table for stellar parameters. Required only for the option crscor so as to compare with a model data.'
  print,'   If this look up table is not provided then the model file should be the a two column model file of wavelength vs flux.'
  print,'   Location of synthetic spectrum (model_file) for cross-correlation. CONTROL assumes that folder are named by their temperature and file are named as model.flx or'
  print,'   a two column model file of wavelength vs flux, if stellar parameters look up table is not provided.'
  print,'10. Flux calibration variables'
  print,'   Location of flux calibration file which provides wavelength vs response relation.'
  print,'  To remove problematic files from processing, create the file badfiles.txt in the data folder'
  print,'*************************************************************************************************************************************************************************'
  
  print,'Operation steps of pipeline is as follows'
  print,' The program works in a series of steps following standard CCD reduction techniques user can also select individual modules using the step function described above.'
  print,' A reduction log is created to help users with processes carried out and mentioning the different parameters used if any.'
  print,' The steps involved in the program execution are as follows:'
  print,''
  print,'   1. Prepare: Check for the different input parameters are set variables accordingly. If data files are present classify them according to file type.'
  print,'   2. Hot and bad pixel correction: Correct for hot and bad pixels in the frames.'
  print,'   3. Create master bias: Create a master bias if master bias file is not found from a set of bias files.'
  print,'   4. Create master dark: Create a master dark if master dark file is not found from a set of dark files.'
  print,'   5. Create master flat: Create a master flat if master flat file is not found from a set of flat files.'
  print,'   6. Correct for bias: Correct of the effect of bias in science frames.'
  print,'   7. Correct for dark: Correct of the effect of dark in science frames.'
  print,'   8. Correct for flat: Correct of the effect of flat in science frames.'
  print,'   9. Correct for cosmic rays: Correct for the effect of cosmic rays using LA cosmic algorithm.'
  print,'   10.Extract spectrum: Define spectrum trace and extract the spectrum.'
  print,'   11.Subtract Background: Subtract background from spectrum.'
  print,'   12.Wavelength calibration: Do wavelength calibration'
  print,'   13.Flux calibration: Do flux calibration'
  print,'   14.Default Light curve: Create 3 light curves (short, middle and long wavelength)'
  print,'   15.Transmission spectrum: Retrie transmission spectrum.'
  print,'***************************************************************************************************************************************************************************'
  print, '  For additional help and more options, go to website:'
  return
endif

if keyword_defined(file) eq 0 then begin
  CD, Current=cur_dir
  file=cur_dir+'control_parameters.txt'
  print,'Configuration file assumed to be in current directory'
endif
if file_test(file) eq 0 then begin
  print,'CONTROL: Missing configuration file. Please re-run CONTROL with a configuration file'
  print,'CONTROL: Users should have received a configuration file with this software distribution if not please refer to the configuratioj file format below.'
  print,'CONTROL: The configuration file should have the folowwing parameters'
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

;reading parameters
if(tag_exist(infile,'data_path') eq 1) then data_path=infile.data_path
if(tag_exist(infile,'temp_path') eq 1) then inter_path=infile.temp_path
if(tag_exist(infile,'out_path') eq 1) then out_path=infile.out_path
if(tag_exist(infile,'bias_combine_type') eq 1) then b_type=infile.bias_combine_type
if(tag_exist(infile,'dark_combine_type') eq 1) then d_type=infile.dark_combine_type
if(tag_exist(infile,'flat_combine_type') eq 1) then f_type=infile.flat_combine_type
if(tag_exist(infile,'cosmic_ray_clip') eq 1) then crclip=infile.cosmic_ray_clip else crclip = 800
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
endif else data_path=detectos(data_path)

logprint,' _______  _______  __    _  _______  ______    _______  ___',logfile=data_path+'control_log'+string(tsec,Format='(I12)')+'.txt'     
logprint,'|       ||       ||  |  | ||       ||    _ |  |       ||   |    
logprint,'|       ||   _   ||   |_| ||_     _||   | ||  |   _   ||   |    
logprint,'|       ||  | |  ||       |  |   |  |   |_||_ |  | |  ||   |    
logprint,'|      _||  |_|  ||  _    |  |   |  |    __  ||  |_|  ||   |___ 
logprint,'|     |_ |       || | |   |  |   |  |   |  | ||       ||       |
logprint,'|_______||_______||_|  |__|  |___|  |___|  |_||_______||_______|
logprint,'---------------------------------------------------------------'
logprint,'    CONTROL v1.0: CUTE AUTONOMOUS DATA REDUCTION PIPELINE      '
logprint,'---------------------------------------------------------------' 
logprint,' 
logprint,'Control execution log dated:',logonly=1
logprint,t,logonly=1  

                                                                                                                 
if dpflg eq 1 then logprint,'CONTROL: No data file path found. data files assumed to be in current directory''
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
    logprint,'CONTROL: No temporary file path found. Temporary file directory created in data directory'
  endif else inter_path=detectos(inter_path)

  
  if n_elements(out_path) eq 0 then begin
        CASE StrUpCase(!Version.OS_Family) OF
          'WINDOWS': out_path=data_path+'output\' ;WINDOWS
          'UNIX': out_path=data_path+'output/'; UNIX. 
        ENDCASE 
    logprint,'CONTROL: No output file path found. Output file directory created in data directory'
  endif else out_path=detectos(out_path)
  
  ;creates paths if they dont exist
  if (file_test(data_path,/DIRECTORY)eq 0) then begin
    logprint,'CONTROL: No data file path found. Creating data file path based on configuration file'
    FILE_MKDIR,data_path
  endif  
  if (file_test(inter_path,/DIRECTORY)eq 0) then begin 
    logprint,'CONTROL: No temporary file path found. Creating temporary file path based on configuration file'
    FILE_MKDIR,inter_path
  endif
  if (file_test(out_path,/DIRECTORY)eq 0) then begin
    logprint,'CONTROL: No output file path found. Creating output file path based on configuration file'
    FILE_MKDIR,out_path
  endif
  
  all_file_list=file_search(data_path+'*.fits') ;Get file list
  ;unique files?
  if (n_elements(all_file_list) eq 0) then begin
    logprint,'CONTROL: File acess error, No fits file found!!',logonly = logonly
    logprint,'CONTROL: Verify data path.'
    message,'File acess error, No fits file found for reduction!!'
    return
  endif
  badfilefile =data_path+'badfiles.txt'  
  if file_test(badfilefile) then readcol, badfilefile, badfiles, format='a' else badfiles=string(0)
  file_list_loc=where(all_file_list ne badfiles)
  file_list=all_file_list[file_list_loc]
  file_names=file_basename(file_list,'.fits') ;Get file names
  file_type=strarr(n_elements(file_list))
  if n_elements(file_list) eq 0 then begin
    logprint,'CONTROL: File acess error, No good fits file found!!',logonly = logonly
    logprint,'CONTROL: Verify data path.'
    message,'File acess error, No good fits file found for reduction!!'
  endif  
  ;carry out bad and hot pixel corrections
  if(hb) then begin
    if n_elements(hbmask) eq 0 then begin
      logprint,'CONTROL: No mask file found for hot and bad pixel correction',logonly = logonly
      message,'CONTROL: No mask file found for hot and bad pixel correction'
    endif  
    hbmask=detectos(hbmask) 
    hbmask=string(hbmask)
    if file_test(hbmask) then begin
      mask=mrdfits(hbmask,0,mask_hdr)
      logprint,'Applying hot and bad pixel correction based on input mask file.'
      for i=0,n_elements(file_list)-1 do begin 
        control_hotbad,file_list[i],mask,out
        hdr=out.header
        out_im=out.data
        sxaddpar, hdr, 'HBFLG', 1. ;Hot and bad pixel correction flag
        writefits,inter_path+file_names[i]+'_hb.fits',out_im,hdr ;hotbad file 
        file_type[i]=SXPAR(hdr, 'FILETYPE');Find object description
      endfor
      file_list_hb=file_search(inter_path+'*_hb.fits')
    endif else begin
      logprint,'No input mask file found for hot and bad pixel correction, proceeding without hot and bad pixel correction.'
      for i=0,n_elements(file_list)-1 do begin
        hdr=HEADFITS(file_list[i],exten=0) 
        file_type[i]=SXPAR(hdr, 'FILETYPE')
        file_list_hb=file_search(data_path+'*.fits')
      endfor
    endelse    
  endif else begin
    logprint,'Proceeding without hot and bad pixel correction as requested by the user.'
    for i=0,n_elements(file_list)-1 do begin
      hdr=HEADFITS(file_list[i],exten=0) 
      file_type[i]=SXPAR(hdr, 'FILETYPE')
      file_list_hb=file_search(data_path+'*.fits')
    endfor  
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
  spectra=file_list_hb[where(strpos(file_type, 'OBJECT') ge 0)]
  spectra_name=file_basename(spectra,'.fits')

  ;bias section
  if (crmb) then begin
    if tag_exist(infile,'master_bias_file') eq 0 then begin
      if n_elements(biaslist) eq 0 then begin
        logprint,'CONTROL: File acess error, No master bias or bias fial found',logonly = logonly
        logprint,'CONTROL: Verify data path.'
        message,'CONTROL: File acess error, No master bias or bias fial found!'
      endif
      control_bias_combine,biaslist,mbias,b_type,sat_value,threshold
      mbias_file=inter_path+'mbias.fits' 
      mwrfits,mbias.im,mbias_file,mbias.hdr,/create, /SILENT
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
    endif else begin
      mbfile_path=detectos(infile.master_bias_file)
      mbfile=file_search(mbfile_path)
      if n_elements(mbfile eq 0) then begin 
        logprint,'CONTROL: File acess error, No master bias file found!!',logonly=logonly
        Message,'CONTROL: File acess error, No master bias file found!!' 
      endif else mbias_data=mrdfits(mbfile,0,mbhdr)
      mbias={im:mbias_data,hdr:mbhdr}
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
   endelse
  endif else begin
    if tag_exist(infile,'master_bias_file') ne 0 then begin
      mbfile_path=detectos(infile.master_bias_file)
      mbfile=file_search(mbfile_path)
      if n_elements(mbfile eq 0) then begin
        logprint,'CONTROL: File acess error, No master bias file found!!',logonly=logonly
        Message,'CONTROL: File acess error, No master bias file found!!'
      endif else mbias_data=mrdfits(mbfile,0,mbhdr)
      mbias={im:mbias_data,hdr:mbhdr}
      r=float(SXPAR( mbias.hdr, 'RNOISE'))
    endif  
  endelse  
     
  ;dark section
  if (crmd) then begin
    if tag_exist(infile,'master_dark_file') eq 0 then begin
      if n_elements(darklist) eq 0 then begin
        logprint,'CONTROL: File acess error, No master dark or dark fial found',logonly = logonly
        logprint,'CONTROL: Verify data path.'
        message,'CONTROL: File acess error, No master dark or dark fial found!'
      endif
        control_dark_combine,darklist,d_type,mdark,mbias,sat_value,threshold
        mdark_file=inter_path+'mdark.fits'
        mwrfits,mdark.im,mdark_file,mdark.hdr,/create
        mwrfits,mdark.error,mdark_file,mdark.hdr, /SILENT
        sigma_d=mdark.error
    endif else begin
      mdfile_path=detectos(infile.master_dark_file)
      mdfile=file_search(mdfile_path)
      if n_elements(mdfile eq 0) then begin 
        logprint,'CONTROL: File acess error, No master dark file found!!',logonly=logonly
        message,'CONTROL: File acess error, No master dark file found!!'
        return 
      endif else begin
        fits_info, mdfile, SILENT=silent,N_ext=extension
        mdark_data=mrdfits(mdfile,0,mdhdr)
        nxyd=size(mdark_data)
        if (extension eq 2) then mdark_err=mrdfits(mdfile,1,mdhdr) else begin
          mdark_err=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master dark file does not contain error, assuming error to be zero'
        endelse  
      endelse 
      mdark={im:mdark_data,error:mdark_err,hdr:dhdr} 
    endelse
  endif else begin
    if tag_exist(infile,'master_dark_file') ne 0 then begin
      mdfile_path=detectos(infile.master_dark_file)
      mdfile=file_search(mdfile_path)
      if n_elements(mdfile eq 0) then begin
        logprint,'CONTROL: File acess error, No master dark file found!!',logonly=logonly
        message,'CONTROL: File acess error, No master dark file found!!'
        return
      endif else begin
        fits_info, mdfile, SILENT=silent,N_ext=extension
        mdark_data=mrdfits(mdfile,0,mdhdr)
        nxyd=size(mdark_data)
        if (extension eq 2) then mdark_err=mrdfits(mdfile,1,mdhdr) else begin
          mdark_err=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master dark file does not contain error, assuming error to be zero'
        endelse
      endelse
      mdark={im:mdark_data,error:mdark_err,hdr:dhdr}
    endif
  endelse  
    
  ;flat section
  if (crmf) then begin
    if tag_exist(infile,'master_flat_file') eq 0 then begin
      if n_elements(flatlist) eq 0 then begin
        logprint,'CONTROL: File acess error, No master flat or flat fial found',logonly = logonly
        logprint,'CONTROL: Verify data path.'
        message,'CONTROL: File acess error, No master flat or flat fial found!'
      endif
      control_flat_combine,flatlist,f_type,mflat,mbias,mdark,sat_value,threshold
      mflat_file=inter_path+'mflat.fits'
      mwrfits,mflat.im,mflat_file,mflat.hdr,/create
      mwrfits,mflat.error,mflat_file,mflat.hdr
      nflat=mflat.im
      sigma_fbdn=mflat.error  
    endif else begin
      mffile_path=detectos(infile.master_flat_file)
      mffile=file_search(mffile_path)
      if n_elements(mffile eq 0) then begin 
        print,'CONTROL: File acess error, No master flat file found!!',logonly=logonly
        message,'CONTROL: File acess error, No master flat file found!!'
        return 
      endif  else begin
        fits_info, mffile, SILENT=silent,N_ext=extension
        nflat=mrdfits(mffile,0,mfhdr)
        nxyd=size(mflat_data)
        if (extension eq 2) then sigma_fbdn=mrdfits(mffile,1,mdhdr) else begin
          sigma_fbdn=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master flat file does not contain error, assuming error to be zero'
        endelse
      endelse
    endelse
  endif else begin
    if tag_exist(infile,'master_flat_file') ne 0 then begin
      mffile_path=detectos(infile.master_flat_file)
      mffile=file_search(mffile_path)
      if n_elements(mffile eq 0) then begin
        print,'CONTROL: File acess error, No master flat file found!!',logonly=logonly
        message,'CONTROL: File acess error, No master flat file found!!'
        return
      endif  else begin
        fits_info, mffile, SILENT=silent,N_ext=extension
        nflat=mrdfits(mffile,0,mfhdr)
        nxyd=size(mflat_data)
        if (extension eq 2) then sigma_fbdn=mrdfits(mffile,1,mdhdr) else begin
          sigma_fbdn=make_array(nxyd[1],nxyd[2],value=0.0,/DOUBLE)
          logprint,'CONTROL: Master flat file does not contain error, assuming error to be zero'
        endelse
      endelse
    endif
  endelse 
     
  if (n_elements(spectra) eq 0) then begin
    logprint,'CONTROL: No science files found. CONTROL will exit now.',logonly = logonly
    message,'CONTROL: No science files found. CONTROL will exit now.'
  endif
  ;bias correction carried out earlier as la cosmic requires bias corrected images
  dummy=mrdfits(spectra[0],0,hdr)
  nx=(size(dummy))[1]
  ny=(size(dummy))[2]
  ycut1=SXPAR( hdr, 'YCUT1')
  ycut2=SXPAR( hdr, 'YCUT2')
  ylen=ycut2-ycut1+1
  if (crb) then begin
    for i=0,n_elements(spectra)-1 do begin
      raw_im=mrdfits(spectra[i],0,hdr)
      nx=(size(raw_im))[1]
      ny=(size(raw_im))[2]
      mnx=(size(mbias.im))[1]
      mny=(size(mbias.im))[2]
      ycut1=SXPAR( hdr, 'YCUT1')
      ycut2=SXPAR( hdr, 'YCUT2')
      ylen=ycut2-ycut1+1
      ccd_gain=SXPAR( hdr, 'GAIN')
      if nx eq mnx then begin
        if ny eq mny then begin
          raw_imb=raw_im-mbias.im
          sigma_imb=sqrt(raw_im+r^2)
        endif else begin
          if (ylen ne ny) then mbias_nw=mbias.im[*,ycut1:ycut1+ny-1] else mbias_nw=mbias.im[*,ycut1:ycut2]
          raw_imb=raw_im-mbias_nw
          sigma_imb=sqrt(raw_im+r^2)
        endelse 
      endif else begin
        logprint,'CONTROL: X Diamnesion error with science image and bias image'
        return
      endelse
      ;add error as well
      sxaddpar, hdr, 'BCFLG', 1,'BIAS CORRECTION FLAG' ;bias correction flag
      mwrfits,raw_imb,inter_path+spectra_name[i]+'_b.fits',hdr,/create
      mwrfits,sigma_imb,inter_path+spectra_name[i]+'_b.fits',hdr
    endfor
    rawb_list= inter_path+spectra_name+'_b.fits'
  endif else rawb_list = spectra  
  
  if cr_cor eq 1 then begin
    rawbc_list= inter_path+spectra_name+'_bc.fits'
    ;cosmic ray correction
    if ccd_gain eq 0 then ccd_gain=1
    la_cosmic,rawb_list,outlist=rawbc_list,gain=ccd_gain,readn=r,sigclip=crclip; add other parameters after testing
    crflg=1
  endif else begin
    logprint,'CONTROL: Cosmic ray correction not carried out as per user request.'
    rawbc_list=rawb_list
    crflg=0
  endelse
  

  ;light curve stuff
  lenby3=fix(nx/3)
  lenby33=fix(nx-2*lenby3)
  
  lc1=dblarr(3,n_elements(spectra))
  lc2=dblarr(3,n_elements(spectra))
  lc3=dblarr(3,n_elements(spectra))
  
  ;dark and flat correction
  for i=0,n_elements(spectra)-1 do begin
    raw_imb=mrdfits(rawbc_list[i],0,hdr)
    print,rawbc_list[i]
    sigma_imb=mrdfits(rawb_list[i],1,hdr) ; error after bias correction
    if sigma_imb eq 0 then sigma_imb=dblarr(nx,ny)
    sxaddpar, hdr,'CRFLG',crflg,'COSMIC RAY CORRECTION FLAG'
    nx=(size(raw_imb))[1]
    ny=(size(raw_imb))[2]
    ycut1=SXPAR( hdr, 'YCUT1')
    ycut2=SXPAR( hdr, 'YCUT2')
    ylen=ycut2-ycut1+1
    if (crd) then begin
      dnx=(size(mdark.im))[1]
      dny=(size(mdark.im))[2]
      if nx eq dnx then begin
        if ny eq dny then begin
          raw_imbd=raw_imb-mdark.im
          sigma_imbd=sqrt((sigma_imb)^2 + (sigma_d)^2) ;error after dark correction
        endif else begin
          if (ylen ne ny) then begin
            mdark_nw=mdark.im[*,ycut1:ycut1+ny-1]
            sigma_dnw=sigma_d[*,ycut1:ycut1+ny-1]
          endif else begin
            mdark_nw=mdark.im[*,ycut1:ycut2]
            sigma_dnw=sigma_d[*,ycut1:ycut2]
          endelse
          raw_imbd=raw_imb-mdark_nw
          sigma_imbd=sqrt((sigma_imb)^2 + (sigma_dnw)^2) ;error after dark correction
        endelse
      endif else begin
        logprint,'CONTROL: X Diamnesion error with science image and dark image'
        return
      endelse
      sxaddpar, hdr, 'DCFLG', 1,'DARK CORRECTION FLAG'  ;dark correction flag
      if save_temp_files eq 1 then mwrfits,raw_imbd,inter_path+spectra_name[i]+'_bcd.fits',hdr,/create
      if save_temp_files eq 1 then mwrfits,sigma_imbd,inter_path+spectra_name[i]+'_bcd.fits',hdr
    endif else begin
      raw_imbd=raw_imb
      sigma_imbd=sigma_imb
    endelse
    if(crf) then begin
      fnx=(size(nflat))[1]
      fny=(size(nflat))[2]
      if nx eq fnx then begin
        if ny eq fny then begin
          raw_imbdf=raw_imbd/nflat
          sigma_imbdf=sqrt((sigma_imbd/nflat)^2 + ((raw_imbd/nflat^2)* sigma_fbdn)^2) ;error
        endif else begin
          if (ylen ne ny) then begin
            nflat_nw=nflat[*,ycut1:ycut1+ny-1]
            sigma_fnw=sigma_fbdn[*,ycut1:ycut1+ny-1]
          endif else begin
            nflat_nw=nflat[*,ycut1:ycut2]
            sigma_fnw=sigma_fbdn[*,ycut1:ycut2]
          endelse
          raw_imbdf=raw_imbd/nflat_nw
          sigma_imbdf=sqrt((sigma_imbd/nflat_nw)^2 + ((raw_imbd/nflat_nw^2)* sigma_fnw)^2) ;error
        endelse
      endif else begin
        logprint,'CONTROL: X Diamnesion error with science image and flat image'
        return
      endelse
      sxaddpar, hdr, 'FCFLG', 1 ,'FLAT CORRECTION FLAG' ;flat correction flag
      if save_temp_files eq 1 then mwrfits,raw_imbdf,inter_path+spectra_name[i]+'_bcdf.fits',hdr,/create
      if save_temp_files eq 1 then mwrfits,sigma_imbdf,inter_path+spectra_name[i]+'_bcdf.fits',hdr    
    endif else begin
      raw_imbdf=raw_imbd
      sigma_imbdf=sigma_imb
    endelse
    ;spectrum identification and extraction including background
    in_image={data:raw_imbdf,error:sigma_imbdf,header:hdr}
   
    if (tag_exist(infile,'trace_type') eq 1 )then t_type = infile.trace_type else t_type = 'simple'
    spectrum=control_trace(in_image,infile,t_type,spectra_name[i])
    ext_spectrum=spectrum.data
    hdr=spectrum.header
    tfileds=4
    sxaddpar, hdr,'TFIELDS',tfileds,'Number of fields in each row'
    sxaddpar, hdr,'TTYPE1  ','COUNTS','label for field   1'
    sxaddpar, hdr,'TFORM1  ','17float','Data format of field 1'
    sxaddpar, hdr,'TUNIT1  ','counts','Physical unit of field 1'
    sxaddpar, hdr,'TDISP1  ','Float','Display format for column  1'
;   sxaddpar, hdr,'TNULL1  ',,'Undefined value for column  1'
;   sxaddpar, hdr,'TTYPE2  ',,'label for field   2'
;   sxaddpar, hdr,'TFORM2  ',,'Data format of field 2'
;   sxaddpar, hdr,'TUNIT2  ',,'Physical unit of field 2'
;   sxaddpar, hdr,'TDISP2  ',,'Display format for column  2'
;   sxaddpar, hdr,'TNULL2  ',,'Undefined value for column  2'
;   sxaddpar, hdr,'INHERIT ',,'Inherit the primary header'
;   sxaddpar, hdr,'EXTNAME ',,'Extension name'
;   sxaddpar, hdr,'EXTVER  ',,'Extension version number'
;   sxaddpar, hdr,'ROOTNAME',,'rootname of the observation set'
    if save_temp_files eq 1 then mwrfits,ext_spectrum,inter_path+spectra_name[i]+'_bcdfe.fits',hdr,/create
    if save_temp_files eq 1 then mwrfits,spectrum.error,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
    if save_temp_files eq 1 then mwrfits,spectrum.background,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
    if save_temp_files eq 1 then mwrfits,spectrum.bck_error,inter_path+spectra_name[i]+'_bcdfe.fits',hdr
    ;data:spectrum_val,error:noise,bck:background,bck_error:backg_error,centroid:centroid
    if (bgs) then begin
      spectra_1d=spectrum.data-spectrum.background
      error_1d=sqrt(spectrum.error^2+spectrum.bck_error^2)
      bgflg=1
    endif else begin
      spectra_1d=spectrum.data
      error_1d=sqrt(spectrum.error^2)
      bgflg=0
    endelse
    spectrum_file={data:spectra_1d,error:error_1d}
    sxaddpar, hdr, 'BGFLG', bgflg ,'BACKGROUND CORRECTION FLAG'
    ;defininig header
    ;PCOUNT  =   Size of special data area
    ;GCOUNT  = One data group
;    sxaddpar, hdr,'TFIELDS',tfileds,'Number of fields in each row'
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

       
    if save_temp_files eq 1 then mwrfits,spectrum_file,out_path+spectra_name[i]+'_1d.fits',hdr,/create
    
    ;wavelength calibration
    if (wcl) then begin
      print,spectra_name[i]
      wavecal_type=infile.wavecal_mode
      wcl_spectrum=control_wavecal(spectrum_file.data,hdr,infile,wavecal_type)
      sp_wavelength=wcl_spectrum.wavelength
      sp_data=wcl_spectrum.flux
      sp_error=spectrum_file.error
      sp_hdr=wcl_spectrum.header
      wcl_spectrum_file={wave:sp_wavelength,counts:sp_data,error:sp_error}
      hdr=sp_hdr
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
       mwrfits,wcl_spectrum_file,out_path+spectra_name[i]+'_1dw.fits',hdr
      ;flux calibration
      if (fcl) then begin
        flux_calib_file=detectos(infile.flux_cal)
        if (file_test(flux_calib_file) eq 0) then begin
          logprint,'CONTROL: Flux calibration file not found. Please re-run the pipeline with actual file address',logonly=logonly
          message,'CONTROL: Flux calibration file not found. Please re-run the pipeline with actual file address'
        endif
        readcol,flux_calib_file,fcal_wave,fcal_value,F='D,D' ;flux calibration file file location
        QE=interpol(fcal_value,fcal_wave,sp_wavelength,/SPLINE)
        spectra_1dwf=spectra_1dw*QE
        error_1dwf=error_1dw*QE
        wclf_spectrum_file={wave:sp_wavelength,flux:spectra_1dwf,error:error_1dwf}
        cgplot,sp_wavelength,spectra_1dwf,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,xtitle='wavelength [$\Angstrom$]',ytitle='flux [ergs s!u-1!n cm!u-2!n $\Angstrom$!u-1!n]'
        write_png,out_path+spectra_name[i]+'_1dwf.png',TVRD(/TRUE)
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
         mwrfits,wclf_spectrum_file,out_path+spectra_name[i]+'_1dwf.fits',spectrum_hdr
      endif
    endif
   if (science) then begin
     ;three light curves
     st=0
     en=lenby3-1
     n_w=en-st+1
     x=make_array(n_w, /INDEX, /NOZERO, START=st)
     y=spectra_1dwf[st:en]
     dy=error_1dwf[st:en]
     tp=trapz_error(x,y,dy)
     lc1[1,i]=tp[0]
     lc1[2,i]=tp[1]
     st=lenby3
     en=2*lenby3-1
     n_w=en-st+1
     x=make_array(n_w, /INDEX, /NOZERO, START=st)
     y=spectra_1dwf[st:en]
     dy=error_1dwf[st:en]
     tp=trapz_error(x,y,dy)
     lc2[1,i]=tp[0]
     lc2[2,i]=tp[1]
     st=2*lenby3
     en=nx-1
     n_w=en-st+1
     x=make_array(n_w, /INDEX, /NOZERO, START=st)
     y=spectra_1dwf[st:en]
     dy=error_1dwf[st:en]
     tp=trapz_error(x,y,dy)
     lc3[1,i]=tp[0]
     lc3[2,i]=tp[1]
     lc1[0,i]=float(SXPAR( hdr, 'TIME'))
     lc2[0,i]=float(SXPAR( hdr, 'TIME'))
     lc3[0,i]=float(SXPAR( hdr, 'TIME'))
   endif
   ;three light curves
;   st=0
;   en=lenby3-1
;   n_w=en-st+1
;   x=make_array(n_w, /INDEX, /NOZERO, START=st)
;   y=spectra_1dwf[st:en]
;   dy=error_1dwf[st:en]
;   tp=trapz_error(x,y,dy)
;   lc1[1,i]=tp[0]
;   lc1[2,i]=tp[1]
;   st=lenby3
;   en=2*lenby3-1
;   n_w=en-st+1
;   x=make_array(n_w, /INDEX, /NOZERO, START=st)
;   y=spectra_1dwf[st:en]
;   dy=error_1dwf[st:en]
;   tp=trapz_error(x,y,dy)
;   lc2[1,i]=tp[0]
;   lc2[2,i]=tp[1]
;   st=2*lenby3
;   en=nx-1
;   n_w=en-st+1
;   x=make_array(n_w, /INDEX, /NOZERO, START=st)
;   y=spectra_1dwf[st:en]
;   dy=error_1dwf[st:en]
;   tp=trapz_error(x,y,dy)
;   lc3[1,i]=tp[0]
;   lc3[2,i]=tp[1]
;   lc1[0,i]=float(SXPAR( hdr, 'TIME'))
;   lc2[0,i]=float(SXPAR( hdr, 'TIME'))
;   lc2[0,i]=float(SXPAR( hdr, 'TIME'))
;   
 endfor 
 if (science ne 0) then begin
   light_curve={time:lc1[0,*],lcs_data:lc1[1,*],lcs_error:lc1[2,*],lcm_data:lc2[1,*],lcm_error:lc2[2,*],lcl_data:lc3[1,*],lcl_error:lc3[2,*]}
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
   mwrfits,light_curve,out_path+'light_curve.fits',lc_hdr
   trlen=n_elements(spectra)
   trlenby4=fix(trlen/4)
   ph_st1=lc1[1,0:trlenby4]
   ph_en1=lc1[1,nx-trlenby4:nx-1]
   out_tran1=[ph_st1,ph_en1]
   mean_ph1=mean(out_tran1)
   tlc1=lc1[0,*]
   nlc1=lc1[1,*]/mean_ph1
   elc1=lc1[2,*]/mean_ph1
   
   ph_st2=lc2[1,0:trlenby4]
   ph_en2=lc2[1,nx-trlenby4:nx-1]
   out_tran2=[ph_st2,ph_en2]
   mean_ph2=mean(out_tran2)
   tlc2=lc2[0,*]
   nlc2=lc2[1,*]/mean_ph2
   elc2=lc2[2,*]/mean_ph2
  
   ph_st3=lc3[1,0:trlenby4]
   ph_en3=lc3[1,nx-trlenby4:nx-1]
   out_tran3=[ph_st3,ph_en3]
   mean_ph3=mean(out_tran3)
   tlc3=lc3[0,*]
   nlc3=lc3[1,*]/mean_ph3
   elc3=lc3[2,*]/mean_ph3
   !P.Multi=[0,1,2]
   cgplot,tlc1,nlc1,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc1,$
     ERR_YHigh=elc1,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Short wavelengths'
     
     cgplot,tlc2,nlc2,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc2,$
     ERR_YHigh=elc2,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Middle wavelengths'
     
   cgplot,tlc3,nlc3,psym=16,symsize=2,charsize=2,charthick=1.5,xthick=1.5,ythick=1.5,ERR_YLow=elc3,$
     ERR_YHigh=elc3,yrange=[0.93,1.03],xtitle='time',ytitle='relative flux',Label='Long wavelengths'
     write_png,out_path+'light_curve.png',TVRD(/TRUE)
     !P.Multi=0
 endif  
  
 logprint,'CONTROL wrapping up the data reduction procedure.'
 logprint,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',close=close
close,/all
end