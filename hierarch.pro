FUNCTION MAKE_CORRECT_TAG,NAME

  GOOD_CHARS=string([(byte('a'))(0)+bindgen((byte('z')-byte('a'))(0)+1B), $
                     (byte('A'))(0)+bindgen((byte('Z')-byte('A'))(0)+1B), $
                     (byte('0'))(0)+bindgen(10B), byte('_')])
  TAG=BYTE(NAME)
  FOR I=0,N_ELEMENTS(TAG)-1 DO $
    IF(STRPOS(GOOD_CHARS,TAG[I]) LT 0) THEN TAG[I]=byte('_')
  TAG=STRING(TAG)

  RETURN,TAG
END

FUNCTION HIERARCH, HDR, NAME, ABORT, COUNT=MATCHES, COMMENT=COMMENT, ERRMES=ERRMES
;+
;
; Name        :
;   HIERARCH()
; Purpose     :
;   Obtain the value of a parameter in a FITS header.
; Explanation :
;   This program is an extension of FXPAR to handle unconventional length
;   keywords in the FITS headers. For more comments see FXPAR.
; Use         :
;   Result = HIERACH( HDR, NAME  [, ABORT [, COUNT=COUNT [, COMMENT=COMMENT [, ERRMES=ERRMES ]]]] )
;
;   Result = HIERACH(HEADER,'DATE')       ;Finds the value of DATE
;   Result = HIERACH(HEADER,'NAXIS*')     ;Returns array dimensions as vector
;
; Inputs      :
;   HDR = FITS header string array. Each
;         element should have a length of 80 characters
;   NAME    = String name of the parameter to return.  If NAME is of the
;         form 'keyword*' then a structure is returned (see below). If
;         the keyword is 'HISTORY' or 'COMMENT' a string array is
;         returned.
; Opt. Inputs :
;   ABORT   = String specifying that FXPAR should do a RETALL if a
;         parameter is not found.  ABORT should contain a string to be
;         printed if the keyword parameter is not found.  If not
;         supplied, HIERACH will return with a negative !err if a keyword
;         is not found.
;   COMMENT = string containg comment, only for a single parameter.
;   COUNT   = Integer specifying the number of parameters returned. If COUNT
;         is 1 the scalar is returned. If larger than 1 the result is packed in
;         a structure (see below). If COUNT is zero, not such keyword was found
;   ERRMES  = optional string containg an error message
; Outputs     :
;   The returned value of the function is the value(s) associated with the
;   requested keyword in the header array.
;
;   If the parameter is complex, double precision, floating point, long or
;   string, then the result is of that type.  Apostrophes are stripped from
;   strings.  If the parameter is logical, 1 is returned for T, and 0 is
;   returned for F.
;
;   If NAME was of form 'keyword*' then a structure in the form:
;   {keyword,value,'C'+keyword,comment,...} is returned. That is different
;   from FXPAR convention but allows to retrieve similar keywords with
;   different value types.
;
;   If the keyword is 'COMMENT' or 'HISTORY' a string array is returned.
;
; Keywords    :
;   COUNT   = Optional keyword to return a value equal to the number of
;         parameters found by HIERARCH;
;   COMMENTS= Array of comments associated with the returned values.
; Calls       :
;   GETTOK, VALID_NUM, MAKE_CORRECT_TAG
; Common      :
;   None.
; Restrictions:
;   None.
; Side effects:
;   Keyword COUNT returns the number of parameters found.
;
;   The system variable !err is set to -1 if parameter not found, 0 for a
;   scalar value returned.  If a vector is returned it is set to the number
;   of keyword matches found.
;
;   If a keyword occurs more than once in a header, a warning is given,
;   and the first occurence is used.  However, if the keyword is "HISTORY",
;   "COMMENT", or "        " (blank), then multiple values are returned.
; Category    :
;   Data Handling, I/O, FITS, Generic.
; Written     :
;   Nikolai Piskunov, UAO, based on FXPAR from the standard FITS library
; Modified    :
;
;-
;------------------------------------------------------------------------------
;
;  Check the number of parameters.
;
    ERRMES = ''
    IF(N_PARAMS() LT 2) THEN BEGIN
      PRINT,'Syntax:  result =  HIERACH(HDR, NAME [,COUNT=COUNT[,ABORT[,COMMENT=COMMENT[,ERRMES=ERRMES]]]])'
      ERRMES = 'Not all mandatory parameters were set'
      RETURN, -1
    ENDIF
;
;  Determine the abort condition.
;
    VALUE = 0
    IF(N_PARAMS() LE 2) THEN BEGIN
      ABORT_RETURN = 0
      ABORT = 'FITS Header'
    END ELSE ABORT_RETURN = 1
    IF(ABORT_RETURN) THEN ON_ERROR,1 ELSE ON_ERROR,2
;
;  Check for valid header.  Check header for proper attributes.
;
    S = SIZE(HDR)
    IF((S[0] NE 1) OR (S[2] NE 7)) THEN $
      ERRMES = 'FITS Header (first parameter) must be a string array'
;
;  Convert the selected keyword NAME to uppercase.
;
    NAM = STRTRIM(STRUPCASE(NAME))
;
;  Determine if NAME is of form 'keyword*'.  If so, then set the VECTOR flag.
;  One must consider the possibility that NAM is an empty string.
;
    NAMELENGTH1 = (STRLEN(NAM) - 1) > 1
    IF(STRPOS(NAM, '*')) EQ NAMELENGTH1 THEN BEGIN
      NAM1 = STRMID( NAM, 0, NAMELENGTH1)
      VECTOR = 1              ;Flag for vector output
;
;  Otherwise, extend NAME with blanks to eight characters.
;
    ENDIF ELSE BEGIN
      VECTOR = 0
    ENDELSE
;
;  If of the form 'keyword*', then find all instances of 'keyword' followed by
;  a number.  Store the positions of the located keywords in NFOUND, and the
;  value of the number field in NUMBER.
;
    IND_KEYWORDS = WHERE(STRPOS(HDR,'=') GT 0, N_KEY)
    IF(N_KEY GT 0) THEN BEGIN
      KEYWORD = HDR(IND_KEYWORDS)
      HHDR = KEYWORD
      FOR I = 0,N_KEY-1 DO KEYWORD[I] = STRMID(KEYWORD[I], 0, STRPOS(KEYWORD[I],'='))
    ENDIF
    IF(NAM EQ 'HISTORY' OR NAM EQ 'COMMENT') THEN BEGIN
      NFOUND = WHERE(STRPOS(HDR,NAM) EQ 0, MATCHES)
      IF(MATCHES GT 0) THEN HHDR=HDR[NFOUND]
    ENDIF ELSE IF(VECTOR) THEN BEGIN
      NFOUND = WHERE(STRPOS(KEYWORD,NAM1) GE 0, MATCHES)
      if(MATCHES gt 0) THEN BEGIN
        K=SORT(KEYWORD[NFOUND])
        KWRD=KEYWORD[NFOUND[K]]
;        HHDR=HHDR[NFOUND[K]]
        KK=UNIQ(KWRD,K)
        IF(N_ELEMENTS(KK) LT MATCHES) THEN BEGIN  ; Found identical keywords.
          J=0L                                    ; Take the first one and complain
          FOR II = 0,N_ELEMENTS(KK)-1 DO BEGIN
            III=WHERE(KEYWORD EQ KWRD[II], NIII)
            IF(NIII GT 1) THEN BEGIN
              MESSAGE,/INFORMATIONAL, 'WARNING- Keyword ' + KWRD[II] + $
                                      ' located more than once in ' + ABORT
              J=[J,III[0]]
            ENDIF ELSE J=[J,III[0]]
          ENDFOR
          MATCHES=N_ELEMENTS(J)-1
          NFOUND=J[1:MATCHES-1]
        ENDIF
;        IF(MATCHES GT 0) THEN NFOUND = IND_KEYWORDS[NFOUND]
      ENDIF
;
;  Otherwise, find all the instances of the requested keyword.  If more than
;  one is found, and NAME is not one of the special cases, then print an error
;  message.
;
    ENDIF ELSE BEGIN
      IF(STRLEN(NAM) LT 8) THEN $
        NFOUND = WHERE(STRPOS(KEYWORD,NAM+' ') EQ 0, MATCHES) $
      ELSE $
        NFOUND = WHERE(STRPOS(KEYWORD,NAM) EQ 0, MATCHES)
      IF((MATCHES GT 1) AND (NAM NE '')) THEN    $
      MESSAGE,/INFORMATIONAL, 'WARNING- Keyword ' +   $
                              NAM + ' located more than once in ' + ABORT
    ENDELSE
;
;  Extract the parameter field from the specified header lines.  If one of the
;  special cases, then done.
;
    IF(MATCHES GT 0) THEN BEGIN
      LINE = HHDR[NFOUND]
      SVALUE = LINE
      FOR I = 0, MATCHES-1 DO SVALUE[I] = STRMID(LINE[I], STRPOS(LINE[I],'=')+1)
      SVALUE = STRTRIM(SVALUE,2)
      IF(NAM EQ 'HISTORY ' OR NAM EQ 'COMMENT ' OR NAM EQ '        ') THEN BEGIN
        VALUE = STRTRIM(STRMID(LINE,8,72),2)
        COMMENTS = STRARR(N_ELEMENTS(VALUE))
;
;  Otherwise, test to see if the parameter contains a string, signalled by
;  beginning with a single quote character (') (apostrophe).
;
      ENDIF ELSE FOR I = 0,MATCHES-1 DO BEGIN
        IF(STRMID(SVALUE[I],0,1) EQ "'" ) THEN BEGIN
          TEST = STRMID(SVALUE[I],1,STRLEN( SVALUE[I] )-1)
          NEXT_CHAR = 0
          VALUE = ''
;
;  Find the next apostrophe.
;
NEXT_APOST:
          ENDAP = STRPOS(TEST, "'", NEXT_CHAR)
          IF(ENDAP LT 0) THEN MESSAGE,     $
            'WARNING: Value of '+NAME+' invalid in '+ABORT+ " (no trailing ')", /info
            VALUE = VALUE + STRMID(TEST, NEXT_CHAR, ENDAP-NEXT_CHAR )
;
;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.
;
          IF(STRMID(TEST,ENDAP+1,1) EQ "'") THEN BEGIN
            VALUE = VALUE + "'"
            NEXT_CHAR = ENDAP+2
            GOTO, NEXT_APOST
          ENDIF
;
;  Extract the comment, if any.
;
          SLASH = STRPOS(TEST, "/", ENDAP)
          IF(SLASH LT 0) THEN COMMENT = '' $
          ELSE                COMMENT = STRMID(TEST, SLASH+1, STRLEN(TEST)-SLASH-1)
;
;  If not a string, then separate the parameter field from the comment field.
;
        ENDIF ELSE BEGIN
          TEST = SVALUE[I]
          SLASH = STRPOS(TEST, "/")
          IF(SLASH GT 0) THEN BEGIN
            COMMENT = STRMID(TEST, SLASH+1, STRLEN(TEST)-SLASH-1)
            TEST = STRMID(TEST, 0, SLASH)
          ENDIF ELSE COMMENT = ''
;
;  Find the first word in TEST.  Is it a logical value ('T' or 'F')?
;
          TEST2 = TEST
          VALUE = GETTOK(TEST2,' ')
          TEST2 = STRTRIM(TEST2,2)
          IF(VALUE EQ 'T') THEN BEGIN
            VALUE = 1
          ENDIF ELSE IF(VALUE EQ 'F') THEN BEGIN
            VALUE = 0
          ENDIF ELSE BEGIN
;
;  Test to see if a complex number.  It's a complex number if the value and the
;  next word, if any, both are valid numbers.
;
            IF(STRLEN(TEST2) EQ 0) THEN GOTO, NOT_COMPLEX
            VALUE2 = GETTOK(TEST2,' ')
            IF(VALID_NUM(VALUE,VAL1) AND VALID_NUM(VALUE2,VAL2)) THEN BEGIN
              VALUE = COMPLEX(VAL1,VAL2)
              GOTO, GOT_VALUE
            ENDIF
;
;  Not a complex number.  Decide if it is a floating point, double precision,
;  or integer number.  If an error occurs, then a string value is returned.
;  If the integer is not within the range of a valid long value, then it will
;  be converted to a double.
;
NOT_COMPLEX:
            ON_IOERROR, GOT_VALUE
            VALUE = TEST
            IF(NOT VALID_NUM(VALUE)) THEN GOTO, GOT_VALUE
            IF(STRPOS(VALUE,'.') GE 0 OR STRPOS(VALUE,'E') GE 0 OR $
               STRPOS(VALUE,'D') GE 0) THEN BEGIN
                IF(STRPOS(VALUE,'D') GT 0 OR STRLEN(VALUE) GE 8) THEN BEGIN
                  VALUE = DOUBLE(VALUE)
                ENDIF ELSE VALUE = FLOAT(VALUE)
            ENDIF ELSE BEGIN
              LMAX = 2.0D^31 - 1.0D
              LMIN = -2.0D31
              VALUE = DOUBLE(VALUE)
              IF(VALUE GE LMIN) AND (VALUE LE LMAX) THEN VALUE = LONG(VALUE)
            ENDELSE

;
GOT_VALUE:
            ON_IOERROR, NULL
          ENDELSE
        ENDELSE     ; if string
;
;  Add to vector if required.
;
        IF(VECTOR) THEN BEGIN
          IF(I GT 0) THEN BEGIN
            VTAG=STRTRIM(KEYWORD[NFOUND[I]],2)
            VTAG=MAKE_CORRECT_TAG(VTAG)
            RESULT = CREATE_STRUCT(RESULT,VTAG,VALUE,'C'+VTAG,COMMENT)
          ENDIF ELSE BEGIN
            VTAG=STRTRIM(KEYWORD[NFOUND[I]],2)
            VTAG=MAKE_CORRECT_TAG(VTAG)
            RESULT = CREATE_STRUCT(VTAG,VALUE,'C'+VTAG,COMMENT)
          ENDELSE
        ENDIF ELSE BEGIN
          RESULT = VALUE
        ENDELSE
      ENDFOR
;
;  Set the value of !ERR for the number of matches for vectors, or simply 0
;  otherwise.
;
      IF(VECTOR) THEN BEGIN
        ERRMES=''
        !ERR = MATCHES
        RETURN, RESULT
      ENDIF ELSE BEGIN
        ERRMES = 'Keyword '+NAM+' not found'
        !ERR = 0
      ENDELSE
;
;  Error point for keyword not found.
;
    ENDIF ELSE BEGIN
      IF ABORT_RETURN THEN MESSAGE,'Keyword '+NAM+' not found in '+ABORT
      ERRMES = 'Keyword '+NAM+' not found'
      !ERR = -1
    ENDELSE
;
  RETURN, VALUE
END
