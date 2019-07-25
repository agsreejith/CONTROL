;+++++++++++++++++++++++++++++++++++++++++++++
; NAME:
;   detectos
;
; PURPOSE:
;   To check the OS and make sure a path is correct for that OS
;
; CALLING SEQUENCE:
;   path = detectos(path)
;
; INPUTS:
;   A string
;
; KEYWORD PARAMETERS:
;
;   path: This is the string which we want to be correctly formated for the OS
;   
; OUTPUTS:
;   The same string as inputed with the correct / or \ 
;   
; COMMON BLOCKS:
;   None
;
; CALLS TO:
;   None
;
; MODIFICATION HISTORY:
;     Written 28/10/2011. S M Elvidge 
;
;---------------------------------------------------------

FUNCTION detectos,path

; !version.os_family will output either 'unix' or 'Windows'

IF (!version.os_family EQ 'Windows') THEN BEGIN
  WHILE ((i = STRPOS(path, '/')) NE -1) DO BEGIN
    STRPUT, path, '\', i
  ENDWHILE
ENDIF ELSE BEGIN
  WHILE ((i = STRPOS(path, '\')) NE -1) DO BEGIN
    STRPUT, path, '/', i
  ENDWHILE
ENDELSE

RETURN,path

END

