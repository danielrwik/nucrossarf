pro nucrossarf__header,dir,evtfile,imhead

if file_test(dir+'/imhead.txt') then begin

l=''
undefine,imhead
openr,lun,dir+'/imhead.txt',/get_lun
while ~eof(lun) do begin
    readf,lun,l
    push,imhead,l
endwhile
free_lun,lun

endif else begin

h=headfits(evtfile,exten=1)
i=1
check=1
while i lt 100 and check ne 0 do begin
    if str(sxpar(h,'TTYPE'+str(i))) eq 'X' and $
          str(sxpar(h,'TTYPE'+str(i+1))) eq 'Y' then check=0
    i++
endwhile
if i ge 100 then $
      stop,'NUIMYLZE__HEADER: failed to find astrometry in evtfile '+evtfile
xx=str(i-1)
yy=str(i)
imhead=headfits(evtfile,exten=0)
sxaddpar,imhead,'BITPIX',32
sxaddpar,imhead,'NAXIS',2
sxaddpar,imhead,'NAXIS1',fix(sxpar(h,'TLMAX'+xx))-fix(sxpar(h,'TLMIN'+xx))+1
sxaddpar,imhead,'NAXIS2',fix(sxpar(h,'TLMAX'+yy))-fix(sxpar(h,'TLMIN'+yy))+1
sxaddpar,imhead,'CTYPE1',sxpar(h,'TCTYP'+xx)
sxaddpar,imhead,'CTYPE2',sxpar(h,'TCTYP'+yy)
sxaddpar,imhead,'CRPIX1',sxpar(h,'TCRPX'+xx)
sxaddpar,imhead,'CRPIX2',sxpar(h,'TCRPX'+yy)
sxaddpar,imhead,'CRVAL1',sxpar(h,'TCRVL'+xx)
sxaddpar,imhead,'CRVAL2',sxpar(h,'TCRVL'+yy)
sxaddpar,imhead,'CDELT1',sxpar(h,'TCDLT'+xx)
sxaddpar,imhead,'CDELT2',sxpar(h,'TCDLT'+yy)

openw,lun,dir+'/imhead.txt',/get_lun
for i=0,n_elements(imhead)-1 do printf,lun,imhead[i]
free_lun,lun

endelse

end
