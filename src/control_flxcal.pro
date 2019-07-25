pro control_flxcal

data=mrdfits(file,0,hdr)
w1=data.wave
f1=data.counts
ref_file=detectos(ref_file)
if readfits(ref_file,h,0) eq -1 then readcol,ref_file,refw,reff,'D,D' else begin
  reference=mrdfits(ref_file,0,hdr)
  refw=reference.wave
  reff=reference.response
endelse
;if (total(w1(0:1)-refw(0:1)) ne 0.0) then begin
;       f_new=interpol,f1,w1,refw
;       print,'warning: 1st spectrum resampled to match 2nd spectrum'
;    endif else f_new=f1
l=fix(w1[n_elements(w1)-1])-w1[0]
new_wave=make_array(l,/INDEX, /NOZERO,start=fix(w1[0])

f_spline=spline(w1,f1,new_wave)
;add remove lines
f_new=interpol(f_spline,new_wave,refw,/SPLINE)
response=reff/f_new
writecol,flux_cal,refw,response,fmt='(2(F17.7,1x))
;*1.988e-8/wavelength

end