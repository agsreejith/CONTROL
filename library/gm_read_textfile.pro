

function gm_read_textfile,filename,n_lines,ignore_blank=ignore_blank,ignore_hash=ignore_hash
    ;/ignore_hash - ignore # lines
    ;/ignore_blank - ignore blank lines
    ;option n_lines - read n lines

  if (file_info(filename)).exists eq 0 then return,'no file'
  if n_elements(n_lines) eq 0 then n_lines=file_lines(filename)
  if n_lines eq 0 then return,''
  s=strarr(n_lines)
  openr,lun,filename,/get_lun
  readf,lun,s
  free_lun,lun
  
  if keyword_set(ignore_blank) then begin
                q=where(s ne '')
                if q[0] ne -1 then s=s[q] else s='blank file'
            end
            
  if keyword_set(ignore_hash) then begin
                q=where(strmid(s,0,1) ne '#')
                if q[0] ne -1 then s=s[q] else s='blank file'
            end
  
  return,s
end
