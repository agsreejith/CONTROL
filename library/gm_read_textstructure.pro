


; gm_read_textstructure.pro 
; Greg Michael 2015
; 
; purpose
;   To read a text file with 'keyword = value' pairs into a structure.
;   Should handle:
;     1. scalar values: keyword = value
;     2. array values: keyword = [val1,val2,...]
;     3. variable length arrays: keyword = *[val1,val2,...]
;     4. comments: # comment
;     5. tables:
;     
;             tablename = {column_name1,column_name2,column_name3 [,...]
;             val1    val2    val3  [... ]
;             val4    val5    val6
;             }
;     
;         result: tablename.column_namex=column
;     
;     6. text blocks:
;     
;            textblockname = {
;            the quick brown
;            fox jumps
;            over a lazy dog
;            }
;            
;            result: textblockname = pointer to string array
;     
;     
;     7. structure arrays: 
;     
;             object = name1
;                 tag1:val1
;                 tag2:val2
;             end_object = name1
;             
;             object = name1
;                 tag1:val3
;                 tag2:val4
;             end_object = name1         
;             
;        Note that items are not numbered, so may be easily added/removed. Sharing a common object
;        name indicates membership of the same array. The array structure will contain the union of
;        tags present in individual members                 
;             
;             
; 
; requirements
;   cmset_op (C. Markwardt - http://www.physics.wisc.edu/~craigm/idl/)
;   gm_read_textfile
; 
;
; example
;   a=gm_read_textstructure("sample_structured_textfile.txt")
;
; modifications
;   2010 written, GM
;   2011-04-19 changed table syntax to work with or without column tags
;   2015-11-05 change table syntax - no column tags => import as text block 
;              (may break old code! old version retained as gm_read_textstructure0.pro)
;-----------------------------------------------------------------------------------------------------------------------
;Copyright (c) 2010, Greg Michael
;All rights reserved.
;
;Redistribution and use in source and binary forms, with or without modification, are permitted 
;provided that the following conditions are met:
;
;       1. Redistributions of source code must retain the above copyright notice, this list of 
;       conditions and the following disclaimer.
;       
;       2. Redistributions in binary form must reproduce the above copyright notice, this list of 
;       conditions and the following disclaimer in the documentation and/or other materials 
;       provided with the distribution.
;       
;THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
;ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
;WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
;DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR 
;ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
;(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
;OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
;THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
;NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
;IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;-----------------------------------------------------------------------------------------------------------------------



function gmrts_strip_quotes,s
  q=where(strmid(s,0,1) eq "'" and strmid(s,0,1,/rev) eq "'")
  if q[0] ne -1 then s[q]=strmid(s[q],1,transpose(strlen(s[q]))-2)
  return,s
end

function gmrts_evaluate,s0
    q=where(strtrim(s0,2) ne '')
    if q[0] eq -1 then return,-1
    s=s0[q] 
    struct=-1 ;empty value
    
    i=0
    while i lt n_elements(s) do begin ;while still more '=' left...
    
      p=strpos([s[i:*],"="],"=") ;use lines from one '=' up to the next            
      q=where(p ne -1)
      if n_elements(q) lt 2 then break ;no more definitions to process
      
      q0=q[0]
      keyword=strtrim(strmid(s[i+q0],0,p[q0]),2)
      remainder=strtrim(strmid(s[i+q0],p[q0]+1),2)
      
      v0=strjoin(strtrim(s[i+q0:i+q[1]-1],2),' ') ;join up intermediate lines
      v=strtrim(strmid(v0,strpos(v0,'=')+1),2) ;remove keyword= to get clean value

        
      ptr=strmid(v,0,1) eq '*' ;store entries beginning * using pointer
      if ptr then v=strmid(v,1)
    
      case strmid(v,0,1) of
          '[':begin
                  w=stregex(v,'\[.*\]',/extract)
                  r=strsplit(strmid(w,1,strlen(w)-2),',',/extract)
                  r=gmrts_strip_quotes(r)
               end
          '(':begin
                  w=stregex(v,'\(.*\)',/extract)
                  r=strsplit(strmid(w,1,strlen(w)-2),',',/extract)
               end
          '"':begin ;double quote
                case remainder eq '"' of ;qutoed string, or start of text block?
                  0:begin
                      w=stregex(v,'".*"',/extract)
                      r=strmid(w,1,strlen(w)-2)                    
                    end
                  1:begin
                      qt=where(strtrim(s[i+q0:*],0) eq '"') ;lone quote as 1st char (trim trailing spaces)
                      if qt[0] ne -1 then begin
                        q[1]=qt[0]+1 ;reset end of entry position to point to next line
                        r=s[i+q0+1:i+q0+qt[0]-1]
                        ptr=1
                      endif               
                    end
                endcase                   
              end
          "'":begin ;single quote
                  w=stregex(v,"'.*'",/extract)
                  r=strmid(w,1,strlen(w)-2) 
               end
          '{':begin ;declare ascii table                    
                  tag0=strsplit(s[i+q0],'{',/extract)
                  tag=n_elements(tag0) eq 1?'':strsplit(tag0[1],',',/extract)  ;column tags provided?                    
                  
                  qt=where(strpos(s[i+q0:*],'}') ne -1)                  
                  if qt[0] ne -1 then begin
                    q[1]=qt[0]+1 ;reset end of entry position to point to next line
                    tab0=s[i+q0+1:i+q0+qt[0]-1]
                    rows=n_elements(tab0)                    
            
                    case tag[0] eq '' of
                      1:begin
                          r=tab0
                          ptr=1
                        end
                      0:begin
                          cols=n_elements(tag) ;how many columns? 
                          tab=strarr(cols,rows)
                          for j=0,rows-1 do tab[*,j]=(gm_quoted_strsplit(tab0[j]))[0:cols-1] ;truncate extra columns                              
                          for j=0,cols-1 do r=j eq 0?create_struct(tag[j],reform(tab[j,*])):create_struct(r,tag[j],reform(tab[j,*]))
                        end
                    endcase
                       
                  endif else r="table error"              
               end             
           else:r=v
      endcase    
      
      if strpos(keyword,'^') ne -1 then keyword="ptr_"+strmid(keyword[i],1)
      case i eq 0 of
        1:struct=create_struct(keyword,ptr?ptr_new(r):r)         ;first entry
        0:struct=create_struct(keyword,ptr?ptr_new(r):r,struct)  ;then concatenate
      endcase    
      
      i+=q[1]
    endwhile
    
    return,struct 
end


function gmrts_resolve,s
    q=(where(stregex(s,"^ *OBJECT *= ([^ ]*)",/boolean,/fold_case)))[0]
    
    if q eq -1 then return,gmrts_evaluate(s)
    
    name=(stregex(s[q],"= ([^ ]*)",/extract,/subexpr))[1]
    
    q_same=where(stregex(s,"^ *OBJECT *= "+name+" *$",/boolean,/fold_case))
    n=n_elements(q_same)
    struct=ptrarr(n)
    qs=lindgen(n_elements(s))
    q_remain=qs
    
    for i=0,n-1 do begin
      q_end=(where(stregex(s[q_same[i]:*],"^ *END_OBJECT *= "+name+" *$",/boolean,/fold_case)))[0]
      struct[i]=ptr_new(gmrts_resolve(s[q_same[i]+1:q_same[i]+q_end-1]))
             
    
      q_gone=qs[q_same[i]:q_same[i]+q_end]
      q_remain=cmset_op(q_remain, 'and', /not2, q_gone) 
      
      case i eq 0 of
          1:tmpstr=*struct[0]
          0:begin
                 s_tags=tag_names(*struct[i])
                 new=cmset_op(s_tags,'and',/not2,tag_names(tmpstr),/index)
                 case new[0] eq -1L of
                    1:
                    0:for j=0,n_elements(new)-1 do tmpstr=create_struct(tmpstr,s_tags[new[j]],(*struct[i]).(new[j]))
                 endcase
              end
       endcase
    endfor    
    
    case n eq 1 of
      1:objstr=tmpstr
      0:begin
            objstr=replicate(tmpstr,n)
            for i=0,n-1 do begin                
                struct_assign,*struct[i],tmpstr
                objstr[i]=tmpstr
            endfor
          end
    endcase

    rest=gmrts_resolve(s[q_remain])
    out=size(rest,/type) eq 8?create_struct(name,objstr,rest):create_struct(name,objstr)
    return,out      
end


function gm_read_textstructure,filename,strarr=strarr
    case keyword_set(strarr) of
      1:ss=filename
      0:begin
            ss=''
            for i=0,n_elements(filename)-1 do begin
              s=gm_read_textfile(filename[i],/ignore_hash)   
              ss=[ss,s]
            endfor
         end
    endcase
    return,gmrts_resolve(ss)
end

;would be nice to add nested structure... e.g. a={x:3,y:4} ...syntax probably as IDL
;
;change so that array is made by duplication of ANY entry, i.e.
;scalar=1
;scalar=2
;
;...yields: scalar=[1,2]
;...but what if types disagree? e.g. merge [1,2] with 1?
;
;then struct arr is made by
;chronology={
;  name='Mars2'
;  coefficients=[2.68e-14,6.93e00,0,4.13e-4]
;}
;
;maybe switch text block to
;text_block="
;Once upon a time, in a land far, far away,
;there was a beautiful young girl who had
;become lost while traveling. As the sun set
;and the moon rose, she happened upon a
;splendid castle.
;"
; - requiring the lone " to end the block?


;a=gm_read_textstructure("d:\mydocs\code\idl\file-format\gm_read_textstructure_sample.txt")
;




