pro nucrossarf__libraries,iobs,nucrossarf_val,nucrossarf_obs

cldir=nucrossarf_obs[iobs].obsdir+'/'+nucrossarf_obs[iobs].obsid+'/event_cl/'
obsid=nucrossarf_obs[iobs].obsid
ab=nucrossarf_obs[iobs].ab
xybox=nucrossarf_obs[iobs].xybox

nucrossarf__header,cldir,evtfile,header

bgddir=nucrossarf_val.bgddir
det=fltarr(1000,1000)
for i=0,3 do begin
    file=cldir+'/'+bgddir+'/det'+str(i)+ab+'im.fits'
    if file_test(file) then fits_read,file,im $
          else stop,'NUCROSSARF__LIBRARIES: Det image file '+file+' not found.'
    det+=im
endfor

pa=sxpar(header,'PA_PNT')+1.0
parad=pa*!pi/180.

;;390,300 --> 710,620
;xmin=390
;ymin=300
;dpsf=10
;npsf=320/dpsf
;darf=4
;narf=320/darf

dpsf=10
px=fltarr(1000/dpsf,1000/dpsf)
py=fltarr(1000/dpsf,1000/dpsf)
for i=0,n_elements(px[*,0])-1 do begin
    px[i,*]=i*dpsf+dpsf/2.
    py[*,i]=i*dpsf+dpsf/2.
endfor
pxnew=(px-500.)*cos(parad)-(py-500.)*sin(parad)+500.
pynew=(px-500.)*sin(parad)+(py-500.)*cos(parad)+500.
ii=where(pxnew ge xybox[0] and pxnew le xybox[2] and $
      pynew ge xybox[1] and pynew le xybox[3] and det[pxnew,pynew] gt 0.1)
pxnew=pxnew[ii]
pynew=pynew[ii]

darf=4
ax=fltarr(1000/darf,1000/darf)
ay=fltarr(1000/darf,1000/darf)
for i=0,n_elements(ax[*,0])-1 do begin
    ax[i,*]=i*darf+darf/2.
    ay[*,i]=i*darf+darf/2.
endfor
axnew=(ax-500.)*cos(parad)-(ay-500.)*sin(parad)+500.
aynew=(ax-500.)*sin(parad)+(ay-500.)*cos(parad)+500.
ii=where(axnew ge xybox[0] and axnew le xybox[2] and $
      aynew ge xybox[1] and aynew le xybox[3])
axnew=axnew[ii]
aynew=aynew[ii]


psfenstr=mrdfits(getcaldbfile('psfen',ab,refmjd=nucrossarf_obs[iobs].mjd),1,/silent)
epsf1=psfenstr.energ_lo
epsf2=psfenstr.energ_hi
energy=(epsf1+epsf2)/2.



; Make PSF library

if not file_test(nucrossarf_val.libdir+'/psf/',/directory) then $
      spawn,'mkdir '+nucrossarf_val.libdir+'/psf'
if not file_test(nucrossarf_val.libdir+'/psf/'+obsid,/directory) then $
      spawn,'mkdir '+nucrossarf_val.libdir+'/psf/'+obsid

if file_test(nucrossarf_val.libdir+'/psf/'+obsid+'/optax'+ab+'.sav') then begin
    restore,nucrossarf_val.libdir+'/psf/'+obsid+'/optax'+ab+'.sav'
    savopt=0
endif else begin
    undefine,optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop ;,pa
    savopt=1
endelse

undefine,xpsfs,ypsfs
psfsavs=file_search(nucrossarf_val.libdir+'/psf/'+obsid+'/psf'+ab+'*sav')
if psfsavs[0] ne '' then begin
    for i=0,n_elements(psfsavs)-1 do begin
        s=strsplit(psfsavs[i],'_',/extract)
        ns=n_elements(s)
        push,xpsfs,fix(s[ns-2])
        ss=strsplit(s[ns-1],'.',/extract)
        push,ypsfs,fix(ss[0])
    endfor
endif

print
print,'Creating PSF library for OBSID '+obsid+' '+ab

for i=0,n_elements(pxnew)-1 do begin
    x=pxnew[i]
    y=pynew[i]
    skip=0
    if size(xpsfs,/type) ne 0 then begin
        nearest=min(sqrt((xpsfs-x)^2+(ypsfs-y)^2))
        if nearest lt sqrt(2.)*dpsf/2. then skip=1
    endif
    if not skip then begin
;for i=0,npsf-1 do for j=0,npsf-1 do begin
;    xo=xmin+dpsf*i
;    yo=ymin+dpsf*j
;    x=(xo-500.)*cos(parad)-(yo-500.)*sin(parad)+500.
;    y=(xo-500.)*sin(parad)+(yo-500.)*cos(parad)+500.
;    if x ge xybox[0] and x le xybox[2] and y ge xybox[1] and y le xybox[3] and $
;          not file_test(nucrossarf_val.libdir+'/psf/'+obsid+'/psf'+ab+'_'+ $
;          str(round(x))+'_'+str(round(y))+'.sav') then begin
        psflib=dblarr(325,325,n_elements(energy))
        for e=0,n_elements(energy)-1 do begin
            psf=getnupsf(cldir,obsid,ab,round(x),round(y),energy[e],$
                  optax,xoptmin,yoptmin,pa)  ;,/nodev)
            psflib[*,*,e]=psf
        endfor
        save,psflib,filename=nucrossarf_val.libdir+'/psf/'+obsid+'/psf'+ $
              ab+'_'+str(round(x))+'_'+str(round(y))+'.sav'
    endif
;    if j eq npsf-1 then counter,i+1,npsf,/percent,'  PSF library: '
    counter,i+1,n_elements(pxnew),/percent,'  PSF library: '
endfor

print
print,'PSF library created for OBSID '+obsid+' '+ab
print



; Make ARF library

if not file_test(nucrossarf_val.libdir+'/arf',/directory) then $
      spawn,'mkdir '+nucrossarf_val.libdir+'/arf'
if not file_test(nucrossarf_val.libdir+'/arf/'+obsid,/directory) then $
      spawn,'mkdir '+nucrossarf_val.libdir+'/arf/'+obsid

if size(xoptax,/type) eq 0 then $
      undefine,optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop

undefine,xarfs,yarfs
arfsavs=file_search(nucrossarf_val.libdir+'/arf/'+obsid+'/arf'+ab+'*sav')
if arfsavs[0] ne '' then begin
    for i=0,n_elements(arfsavs)-1 do begin
        s=strsplit(arfsavs[i],'_',/extract)
        ns=n_elements(s)
        push,xarfs,fix(s[ns-2])
        ss=strsplit(s[ns-1],'.',/extract)
        push,yarfs,fix(ss[0])
    endfor
endif

print
print,'Creating ARF library for OBSID '+obsid+' '+ab

for i=0,n_elements(axnew)-1 do begin
    x=axnew[i]
    y=aynew[i]
    skip=0
    if size(xarfs,/type) ne 0 then begin
        nearest=min(sqrt((xarfs-x)^2+(yarfs-y)^2))
        if nearest lt sqrt(2.)*darf/2. then skip=1
    endif 
    if not skip then begin
;for i=0,narf-1 do for j=0,narf-1 do begin
;    xo=xmin+darf*i
;    yo=ymin+darf*j
;    x=(xo-500.)*cos(parad)-(yo-500.)*sin(parad)+500.
;    y=(xo-500.)*sin(parad)+(yo-500.)*cos(parad)+500.
;    if x ge xybox[0] and x le xybox[2] and y ge xybox[1] and y le xybox[3] and $
;          not file_test(nucrossarf_val.libdir+'/arf/'+obsid+'/arf'+ab+'_'+ $
;          str(round(x))+'_'+str(round(y))+'.sav') then begin
        arf=getnuarf(cldir,obsid,ab,round(x),round(y),nucrossarf_obs[iobs].mjd,$
              optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop,pa) ;,/apstopcorr)
        save,arf,filename=nucrossarf_val.libdir+'/arf/'+obsid+'/arf'+ab+'_'+ $
              str(round(x))+'_'+str(round(y))+'.sav'
    endif
;    if j eq narf-1 then counter,i+1,narf,/percent,'  ARF library: '
    counter,i+1,n_elements(axnew),/percent,'  ARF library: '
endfor

print
print,'ARF library created for OBSID '+obsid+' '+ab
print


if savopt then save,optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop, $
      filename=nucrossarf_val.libdir+'/psf/'+obsid+'/optax'+ab+'.sav'


end
