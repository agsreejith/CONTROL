;
;  NAME:
;
;     CROSS_CORRELATE
;     
;  CALLING SEQUENCE:
;        res = cross_correlate(lam_ref, flux_ref, lam_test, flux_test, lmin, lmax, CCF=ccf, GAUSS=gauss$
;                              ,PLOT=plot)
;
;  INPUTS:
;        lam/flux_ref  =  lambda/flux of reference spectrum
;        lam/flux_test =  lambda/flux of test spectrum
;        lmin/lmax     =  min/max wavelength for the cross-correlation
;        
;  OPTIONAL INPUTS: 
;         GAUSS        = Calculates the maximum of the cross-correlation function by a Gaussian 
;                        fitting (default is by linear finite differences)' 
;             
;  OUTPUT:
;         res          = Apparent shift between test and reference spectrum in wavelength.
; 
; OPTIONAL OUTPUT:
;        CCF           = Optional keyword for outputting the normalized cross correlation function'
;

;
;
; PURPOSE:
;
;     Derives the normalized cross correlation function for two spectra and
;     calculate its maximum by Gaussian fitting or finite differences. The
;     location of the maximum represents the apparent shift of spectrum 1 with
;     respect to spectrum 2.
;
;       Based on:
;
;       http://archive.stsci.edu/iue/manual/dacguide/node78.html
;
;
;##################################################################################################



function cross_correlate, lam_ref, flux_ref, lam_test, flux_test, lmin, lmax, CCF=ccf, GAUSS=gauss $
                        , PLOT=plot


if n_params() lt 6 then begin
    print, ''
    print, ' CALLING SEQUENCE: '
    print, ''
    print, ' res = cross_correlate( lam_ref, flux_ref, lam_test, flux_test, lmin, lmax'$
         +' [,  CCF=ccf, /GAUSS, /PLOT ])'
    print, ''
    print, ''
    print, 'lam/flux_ref:  lambda/flux of reference spectrum'
    print, ''
    print, 'lam/flux_test:  lambda/flux of test spectrum'
    print, ''
    print, 'lmin/lmax:  min/max wavelength for the cross-correlation'
    print, ''
    print, 'CCF:  Optional keyword for outputting the normalized cross correlation function'
    print, ''
    print, '/GAUSS: Calculates the maximum of the cross-correlation function by a Gaussian fitting'
    print, '              (default is by linear finite differences)'
    print, ''
    goto, fin
endif


f_st = flux_ref
f_nd = flux_test

lam_st = lam_ref
lam_nd = lam_test

mm1 = where(lam_st ge lmin and lam_st le lmax)
mm2 = where(lam_nd ge lmin and lam_nd le lmax)

lam1 = lam_st[mm1]
lam2 = lam_nd[mm2]

flux1 = f_st[mm1]
flux2 = f_nd[mm2]


;-------------------------- Compute a normalized cross correlation function

nspr = 15
deltaw = lam2[1] - lam2[0]

mi   = n_elements(flux1)
mi2  = n_elements(flux2)
avg1 = total(flux1)/mi & ff= flux1 - avg1   ; subtract average fluxes 
avg2 = total(flux2)/mi & fs= flux2 - avg2
ntot = nspr+nspr+1

cross= fltarr(ntot)
temp = fs
for l=0,ntot-1 do begin
    ns = nspr - l
    temp = fs*shift(ff,ns)
    ls = ns > 0
    us = mi  - 1 + (ns < 0)
    nele = us - ls + 1
    cross(l) = total(temp(ls:us)) / nele
endfor 
crmax= max(cross) & crmin= min(cross)
cross= (cross -crmin)/(crmax-crmin)       ; normalized function 

vspace = deltaw * (indgen(n_elements(cross)) - nspr)

if keyword_set(plot) then begin
    window, 1, xp=0, yp=500
    plot,vspace,cross,psym=6,xrange=[min(vspace),max(vspace)],yrange=[-0.10,1.15], $
      title=' Normalized Cross Correlation Function '
endif



;-------------------------- Find maximum

if keyword_set(gauss) then begin                     ;Gaussian fit

    yfit=gaussfit(float(vspace),cross,a,nterms=3)
    delv=a[1]
    ;sigy=sqrt(total((cross-yfit)*(cross-yfit))/(inp*3+1))

endif else begin                                     ; Max by first differences

    ntot = n_elements(cross)
    j=max(cross)
    ind = where(cross eq j)
    kb=ind(0)-3>0 & ke=ind(0)+3<(ntot-1)
    temp = cross - shift(cross,-1)
    diff = temp(kb:ke-2)
    temp = (vspace + shift(vspace,-1))/2.
    vsub = temp(kb:ke-2)
    nt = n_elements(diff)
    wt = fltarr(nt)
    xmid = float(nt)/2.0 - 0.5
    lfit = linfit(vsub,diff,sigma=si)
    a = lfit[0]
    b = lfit[1]
    delv = -a/b

endelse

;-------------------------- Apply shift + interpolation


if keyword_set(plot) then begin

    print, ''
    ;print, '  Resulting shift is: ',spc_str(delv)+' Angstroms'
    print, ''

    lam3 = lam2 - delv
    flux3 = interpol(flux2,lam2,lam3)

    wset,1
    plots, [delv,delv], [-0.10,1.15]
    if keyword_set(gauss) then oplot, vspace, yfit
    
    window,0
    ll = (max(lam1)-min(lam1))*0.1
    xmin = min(lam1)-ll
    xmax = max(lam1)+ll
    plot, lam1, flux1, xr=[xmin,xmax], xs=1
    oplot, lam2, flux2, color=250
    oplot, lam2, flux3, color=50

endif

ccf = cross

return, delv
fin:
end
