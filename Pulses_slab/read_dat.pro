; 
; NAME: 
;	READ_DAT
;
; PURPOSE:
;	Reading from data file 
;
; EXPLANATION:
;	This procedure reads  from a ASCI data input file. The data are 
; 	stored in the matrix 'data'. 
;
; CALLING SEQUENCE:
; 	read_dat, datafile, data, quiet = quiet
; 
; INPUTS:
; 	datafile = input data file 
;	data : will be defined as multidimensional array
;
; OUTPUTS: 
;	data


PRO read_dat, datafile, data, quiet = quiet

; --- check file format  ---
readerror = 0 
openr, 1 , datafile 
line = 0 
rows = 0 
max_columns = 0 
; 
;
;
while (not eof(1)) do begin 
  b = ' ' 
  readf, 1,  b
  line = line + 1 
  ; -- check if comment line -- 
  comment = (strpos(b,'#')ge 0)or(strpos(b,'$')ge 0)or(strpos(b,'%')ge 0)$
          or(strpos(b,'&')ge 0)or(strpos(b,'*')ge 0)or(strpos(b,';')ge 0)$
          or(strpos(b,':')ge 0)or(strpos(b,'?')ge 0)or(strpos(b,'>')ge 0)
  if comment then begin
    if not keyword_set(quiet) then print, 'line', strcompress(string(line)), ' ' , b 
    endif  else begin 
    ; -- determine number of data per line --
    b = strtrim(b,2)
    b = strcompress(b) 
    len1 = strlen(b)
    if len1 gt 0 then begin 
      len2 = strlen(strcompress(b, /remove_all))
      columns = len1-len2+1
      if (columns ne max_columns) and (max_columns gt 0) $
          then begin 
            print, 'line', strcompress(string(line)), ' wrong number of data points' 
            readerror = 1 
      endif       
      max_columns = columns>max_columns
      rows = rows+1
    endif else begin 
      columns = 0 
      if not keyword_set(quiet) then print, 'line', strcompress(string(line)), ' empty'    
    endelse 
  endelse
endwhile 
close, 1 

; --- read file --- 

if readerror then begin 
  print, 'Data could not be read, because of inconsistent format.' 
endif else begin 
if not keyword_set(quiet) then print, 'data format:', strcompress(string(rows)) , ' rows,' $
		     , strcompress(string(max_columns)), ' columns'
  data = dblarr(rows, max_columns)
  openr, 1, datafile
  i = -1  
  while (not eof(1)) do begin   
    b = ' ' 
    readf, 1,  b
    ; -- check if comment line or empty -- 
    comment = (strpos(b,'#')ge 0)or(strpos(b,'$')ge 0)or(strpos(b,'%')ge 0)$
          or(strpos(b,'&')ge 0)or(strpos(b,'*')ge 0)or(strpos(b,';')ge 0)$
          or(strpos(b,':')ge 0)or(strpos(b,'?')ge 0)or(strpos(b,'>')ge 0)
    empty = strlen(strcompress(b, /remove_all)) eq 0
    if not(comment or empty) then begin 
      i = i + 1
      x0 = data(i,*)       
      reads, b, x0
      data(i,*) = x0 
    endif
  endwhile 
  close, 1 
endelse 

END
