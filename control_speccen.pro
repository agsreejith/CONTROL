pro control_speccen


m=max(im,loc,DIMENSION=2)

;x=nx
;y=loc
fit_out=linfit(nx,loc)
A=fit_out[0]
B=fit_out[1]

columnsum=total(im,1)
max_column=max(columnsum,loc_col)


if boxext eq 1 then begin
  l=loc[1023]-width
  u=loc[1023]-width
  spectrum= total((im)[*,l:u],1)
endif

end

sfdfr