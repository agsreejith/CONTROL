pro errorlog, p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, _extra=ex, logfileset=logfileset, close=close, logonly=logonly

  ; prints to a terminal AND to a log file.

  common log, logfile, loglun

  loglun = 79
  fil = 1
  if (fstat(loglun)).open eq 0 then begin
    if n_elements(logfileset) gt 0 then begin
      logfile = logfileset
      openw, loglun, logfile, /append
    endif else begin
      fil = 0
    endelse
  endif
  scr = keyword_set(logonly) eq 0

  case n_params() of
    0: begin
      if scr then print
      if fil then printf, loglun
    end
    1: begin
      if scr then print, p1, _extra=ex
      if fil then printf, loglun, p1, _extra=ex
    end
    2: begin
      if scr then print, p1, p2, _extra=ex
      if fil then printf, loglun, p1, p2, _extra=ex
    end
    3: begin
      if scr then print, p1, p2, p3, _extra=ex
      if fil then printf, loglun, p1, p2, p3, _extra=ex
    end
    4: begin
      if scr then print, p1, p2, p3, p4, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, _extra=ex
    end
    5: begin
      if scr then print, p1, p2, p3, p4, p5, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, _extra=ex
    end
    6: begin
      if scr then print, p1, p2, p3, p4, p5, p6, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, _extra=ex
    end
    7: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, _extra=ex
    end
    8: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, _extra=ex
    end
    9: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, _extra=ex
    end
    10: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, _extra=ex
    end
    11: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, _extra=ex
    end
    12: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,_extra=ex
    end
    13: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11,p12,p13, _extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,p13, _extra=ex
    end
    14: begin
      if scr then print, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,p13,p14,_extra=ex
      if fil then printf, loglun, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11,p12,p13,p14, _extra=ex
    end
  endcase

  if keyword_set(close) then close, loglun

end

