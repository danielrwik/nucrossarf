pro nucrossarf_arfs,nucrossarf_val,nucrossarf_obs,nucrossarf_src,clobber=clobber

if size(clobber,/type) eq 0 then clobber=0

print
print,'Acquiring ARFs'


nobs=n_elements(nucrossarf_obs)
regnames=nucrossarf_val.regnames
nreg=nucrossarf_val.nreg
if n_elements(regnames) ne nreg then $
      stop,'NUCROSSARF_ARFS: number of regions inconsistent somehow'
dim=float(get_screen_size())
fscl=(sqrt(dim[0]^2+dim[1]^2)/1686.56)^0.58
xyim=fix(600*fscl)
mkspec=nucrossarf_val.mkspec

; don't forget to multiply psf by exposure map so weighting is correct!

for iobs=0,nobs-1 do begin


outdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid
ab=nucrossarf_obs[iobs].ab
obsid=nucrossarf_obs[iobs].obsid
cldir=nucrossarf_obs[iobs].obsdir+'/'+obsid+'/event_cl/'
evtfile=cldir+'nu'+obsid+ab+'01_cl.evt'
nucrossarf__header,cldir,evtfile,head
nucrossarf__updateastrom,iobs,nucrossarf_obs,head,imhead
fits_write,cldir+'dummy.fits',intarr(1000,1000),imhead

regmask=intarr(nreg,1000,1000)
regpsfsave=fltarr(nreg,1000,1000)
for r=0,nreg-1 do begin
    mask=reg2mask(cldir+'dummy.fits',regnames[r])
    regmask[r,*,*]=mask
endfor
;spawn,'rm -f '+cldir+'dummy.fits'

;arffilearray=getcaldbfile('arf',ab)
;arfstr=mrdfits(arffilearray[0],1,harf,/silent)
arfstr=mrdfits(getcaldbfile('arf',ab,refmjd=nucrossarf_obs[iobs].mjd),1,harf,/silent)

bgddir=nucrossarf_val.bgddir
det=fltarr(1000,1000)
for i=0,3 do begin
    file=cldir+'/'+bgddir+'/det'+str(i)+ab+'im.fits'
    if file_test(file) then fits_read,file,im $
          else stop,'NUCROSSARF__ARFS: Det image file '+file+' not found.'
    det+=im
endfor

nmod=n_elements(nucrossarf_obs[iobs].model[0,*,0,0])
order=intarr(nmod)  ;-1)
ii=where(nucrossarf_src.ptorder ge 0)
if ii[0] ne -1 then order[ii]=1
ii=where(nucrossarf_src.extorder ge 0)
if ii[0] ne -1 then order[ii]=2

psfenstr=mrdfits(getcaldbfile('psfen',ab,refmjd=nucrossarf_obs[iobs].mjd),1,/silent)
epsf1=psfenstr.energ_lo
epsf2=psfenstr.energ_hi
epsf=(epsf1+epsf2)/2.
epsflog=10.^((alog10(epsf2)-alog10(epsf1))/2.)*epsf1
nepsf=n_elements(epsf)
eband=['3to4.5','4.5to6','6to8','8to12','12to20','20to79']

earf=findgen(4096)*0.04+1.6


ipt=0
iext=0

xpsf=325
psfsavs=file_search(nucrossarf_val.libdir+'/psf/'+obsid+'/psf'+ab+'*sav')
if psfsavs[0] eq '' then stop,'NUCROSSARF_ARFS: psf library not found'
undefine,xpsfs,ypsfs
for i=0,n_elements(psfsavs)-1 do begin
    s=strsplit(psfsavs[i],'_',/extract)
    ns=n_elements(s)
    push,xpsfs,fix(s[ns-2])
    ss=strsplit(s[ns-1],'.',/extract)
    push,ypsfs,fix(ss[0])
endfor

if min(xpsfs) lt 0 or max(xpsfs) gt 999 or $
      min(ypsfs) lt 0 or max(ypsfs) gt 999 then $
  stop,'NUCROSSARF_ARFS: psf library contains sav files outside of 1000x1000 area'

arfsavs=file_search(nucrossarf_val.libdir+'/arf/'+obsid+'/arf'+ab+'*sav')
if arfsavs[0] eq '' then stop,'NUCROSSARF_ARFS: arf library not found'
undefine,xarfs,yarfs
for i=0,n_elements(arfsavs)-1 do begin
    s=strsplit(arfsavs[i],'_',/extract)
    ns=n_elements(s)
    push,xarfs,fix(s[ns-2])
    ss=strsplit(s[ns-1],'.',/extract)
    push,yarfs,fix(ss[0])
endfor
if min(xarfs) lt 0 or max(xarfs) gt 999 or $
      min(yarfs) lt 0 or max(yarfs) gt 999 then $
  stop,'NUCROSSARF_ARFS: arf library contains sav files outside of 1000x1000 area'

;; meant to speed up arf generation, but limitation is shifting psf image
;;  pixel by pixel -- need instead to more coarsely bin PSF 
;deltaarf=nucrossarf_val.darf
;if deltaarf gt 0 then begin
;    xinit=min(xarfs)
;    yinit=min(yarfs)
;    nx=fix((max(xarfs)-min(xarfs))/deltaarf)+1
;    ny=fix((max(yarfs)-min(yarfs))/deltaarf)+1
;    xarfsorig=xarfs
;    yarfsorig=yarfs
;    undefine,xarfs,yarfs
;    xarfs=xarfsorig[0]
;    yarfs=yarfsorig[0]
;    for i=0,nx-1 do for j=0,ny-1 do begin
;        dist=min(sqrt((i*deltaarf+xinit-xarfsorig)^2+ $
;              (j*deltaarf+yinit-yarfsorig)^2),ii)
;        jj=where(xarfs eq xarfsorig[ii] and yarfs eq yarfsorig[ii])
;        if jj[0] eq -1 then begin
;            push,xarfs,xarfsorig[ii]
;            push,yarfs,yarfsorig[ii]
;        endif
;    endfor
;endif



print
print,'Observation '+str(iobs+1)+' out of '+str(nobs)
print,'  EVT: '+evtfile
print,'    Mapping Library to Pixels of this Observation'



allsrc=intarr(1000,1000)
allsrc[min(xarfs):max(xarfs),min(yarfs):max(yarfs)]=1
ii=where(allsrc gt 0.5)
ii2d=array_indices(allsrc,ii)
psfid=intarr(1000,1000)
arfarrall=intarr(2,n_elements(ii))
nlib=n_elements(ii)
pid=0
undefine,xpsfarr,ypsfarr,xarfarr,yarfarr,pidarr
for i=0L,nlib-1 do begin
    dist=min((ii2d[0,i]-xarfs)^2+(ii2d[1,i]-yarfs)^2,jj)
    arfarrall[*,i]=[xarfs[jj],yarfs[jj]]
    if dist lt 0.1 then begin
        push,xarfarr,xarfs[jj]
        push,yarfarr,yarfs[jj]
        blah=min((ii2d[0,i]-xpsfs)^2+(ii2d[1,i]-ypsfs)^2,kk)
        push,xpsfarr,xpsfs[kk]
        push,ypsfarr,ypsfs[kk]
        pid++
        push,pidarr,pid
    endif
    if i mod long(nlib/1000)*10 eq 0 then $
          counter,i+1,nlib,/percent,'      ARF/PSF binning: '
endfor
print
npid=n_elements(pidarr)
for i=0,npid-1 do begin
    ll=where(xarfarr[i] eq arfarrall[0,*] and yarfarr[i] eq arfarrall[1,*])
    if ll[0] eq -1 then stop,'NUCROSSARF_ARFS: dumb thing happened that shouldnt'
    psfid[ii[ll]]=pidarr[i]
    counter,i+1,n_elements(pidarr),/percent,'        Assigning ARF/PSF pixels: '
endfor
;pid=0
;undefine,xpsfarr,ypsfarr,xarfarr,yarfarr
;for i=min(xarfs),max(xarfs) do begin
;  for j=min(yarfs),max(yarfs) do begin
;    jj=where(arfarrall[0,*] eq i and arfarrall[1,*] eq j)
;    if jj[0] ne -1 then begin
;        pid++
;        psfid[ii2d[0,jj],ii2d[1,jj]]=pid
;        push,xarfarr,i
;        push,yarfarr,j
;        blah=min((mean(ii2d[0,jj])-xpsfs)^2+(mean(ii2d[1,jj])-ypsfs)^2,kk)
;        push,xpsfarr,xpsfs[kk]
;        push,ypsfarr,ypsfs[kk]
;    endif
;  endfor
;  counter,i+1-min(xarfs),max(xarfs)-min(xarfs),/percent,'      ARF/PSF binning: '
;endfor
print
print
if n_elements(xpsfarr) ne pid or n_elements(xarfarr) ne pid then $
      stop,'NUCROSSARF_ARFS: issue with mapping library to pixels'




for isrc=0,nmod-1 do begin

redo=0
for r=0,nreg-1 do begin
    arfname=outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf'
    if not file_test(arfname) then redo=1
endfor

if redo or clobber then begin

regpsfsave=fltarr(nreg,1000,1000)

undefine,offaxis
if order[isrc] eq 1 then begin
    ra=nucrossarf_src[isrc].ptra
    dec=nucrossarf_src[isrc].ptdec
    scl=1.
    ipt++
endif else if order[isrc] eq 2 then begin
    if nucrossarf_src[isrc].extsrcimg eq 'flat' then begin
        srcim=reg2mask(cldir+'dummy.fits',regnames[isrc])
        regim=srcim
        hreg=imhead
        hsrc=imhead
        ii=where(regim gt 0.5)
        if ii[0] eq -1 then $
              stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+' not found in region: '+$
                    nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg
    endif else begin
        fits_read,nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extsrcimg,srcim,hsrc
        fits_read,nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg,regim,hreg
        ii=where(regim eq nucrossarf_src[isrc].extregnum)
        if ii[0] eq -1 then $
              stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+ $
                    ' not found in region image '+$
                    nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg
        regim[*,*]=0
        regim[ii]=1
    endelse
    if nucrossarf_val.darf gt 1 then begin
        hcongrid,regim,hreg,newregim,newhreg, $
              n_elements(regim[*,0])/nucrossarf_val.darf,$
              n_elements(regim[0,*])/nucrossarf_val.darf,interp=0
        hcongrid,srcim,hsrc,newsrcim,newhsrc, $
              n_elements(srcim[*,0])/nucrossarf_val.darf,$
              n_elements(srcim[0,*])/nucrossarf_val.darf,interp=1
        regim=newregim
        hreg=newhreg
        srcim=newsrcim
        hsrc=newhsrc
        ii=where(regim gt 0.5)
    endif
    ii2d=array_indices(regim,ii)
    xyad,hreg,ii2d[0,*],ii2d[1,*],ra,dec
    xyad,hsrc,ii2d[0,0],ii2d[1,0],rach,decch
    if sphdist(rach,decch,ra[0],dec[0],/degrees)*3600. gt 1. then $
          stop,'NUCROSSARF_ARFS: extended source and region image astrometry '+$
          'is inconsistent at >1", should be identical'
    scl=srcim[ii]/total(srcim[ii])
    iext++
endif else stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+' neither pt nor ext???'

adxy,imhead,ra,dec,xo,yo
x=xo
y=yo
;x=(xo-499.5)*cos(-nucrossarf_obs[iobs].rotast*!pi/180.) - $
;      (yo-499.5)*sin(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
;      499.5 + nucrossarf_obs[iobs].xshast
;y=(xo-499.5)*sin(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
;      (yo-499.5)*cos(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
;      499.5 + nucrossarf_obs[iobs].yshast
xs=x
ys=y
scls=scl
if n_elements(x) gt 10. then plotprog=1 else plotprog=0



print,'  Source '+str(isrc+1)+' out of '+str(nmod)


rpsfid=psfid[x[0],y[0]]
idxy=intarr(n_elements(x))
for i=0,n_elements(x)-1 do idxy[i]=psfid[x[i],y[i]]
if n_elements(x) gt 1 then for is=1,n_elements(x)-1 do begin
    thisid=psfid[x[is],y[is]]
    ii=where(rpsfid eq thisid)
    if ii[0] eq -1 then push,rpsfid,thisid
endfor
rpsfid=rpsfid[sort(rpsfid)]



if plotprog then begin
    progim=contimg(det,0.01)
    ii=where(progim gt 0.5)
    ii2d=array_indices(progim,ii)
    proglim=[min(ii2d[0,*]),max(ii2d[0,*]),min(ii2d[1,*]),max(ii2d[1,*])]
;    ii=where(regim gt 0.5)
;    progim[ii]=2
    for i=0,n_elements(x)-1 do progim[x[i]-0.6:x[i]+0.6,y[i]-0.6:y[i]+0.6]=2
endif


undefine,pixperbin
for ilib=0,n_elements(rpsfid)-1 do begin
    ii=where(idxy eq rpsfid[ilib])
    push,pixperbin,n_elements(ii)
endfor


arfs=fltarr(nreg,n_elements(earf))

for ilib=0,n_elements(rpsfid)-1 do begin

    ii=where(idxy eq rpsfid[ilib])
    x=xs[ii]
    y=ys[ii]
    scl=scls[ii]
    if plotprog then t0=systime(/seconds)

;    blah=min((x[is]-xpsfs)^2+(y[is]-ypsfs)^2,ii)
    file=nucrossarf_val.libdir+'psf/'+obsid+'/psf'+ab+'_'+ $
          str(xpsfarr[rpsfid[ilib]])+'_'+str(ypsfarr[rpsfid[ilib]])+'.sav'
    if not file_test(file) then stop,'NUCROSSARF_ARFS: whaaaaaaaaaaaa? (psf)'
; psflib variable restored by below call
    restore,file


;    blah=min((x[is]-xarfs)^2+(y[is]-yarfs)^2,ii)
    file=nucrossarf_val.libdir+'arf/'+obsid+'/arf'+ab+'_'+ $
          str(xarfarr[rpsfid[ilib]])+'_'+str(yarfarr[rpsfid[ilib]])+'.sav'
    if not file_test(file) then stop,'NUCROSSARF_ARFS: whaaaaaaaaaaaa? (arf)'
; arf variable restored by below call
    restore,file


for is=0,n_elements(x)-1 do begin
    temppsf=dblarr(nepsf,1000,1000)
    pweight=fltarr(nepsf,nreg)
    for p=0,nepsf-1 do begin
        psf=reform(psflib[*,*,p])
        temp2psf=dblarr(1000,1000)
        temp2psf[xpsf/2-162:xpsf/2+162,xpsf/2-162:xpsf/2+162]=psf
        temppsf[p,*,*]=fshift(temp2psf,x[is]-162.,y[is]-162.)
        temppsf[p,*,*]*=det
        for r=0,nreg-1 do begin
            mask=reform(regmask[r,*,*])
            pweight[p,r]=total(temppsf[p,*,*]*mask)
            if p eq 0 then regpsfsave[r,*,*]+=temppsf[p,*,*]*mask*scl[is]
        endfor
    endfor

; for each region, fit pweight to linear function
    for r=0,nreg-1 do begin
        pw=pweight[*,r]
        pscl=interpol(pw,epsflog,earf)
        arfs[r,*]+=arf*pscl*scl[is]
    endfor

    if plotprog then progim[x[is]-0.6:x[is]+0.6,y[is]-0.6:y[is]+0.6]=3

; then scale arf by that and scl factor from dist. of emission
; add scaled arf to final arf for each region for this source

endfor   ; is loop

if plotprog then begin
    if ilib eq 0 then chan,0,xyim,xyim
    erase
    tv,congrid(progim[proglim[0]:proglim[1],proglim[2]:proglim[3]],xyim,xyim),$
          xsize=xyim,ysize=xyim
    xyouts,0.02,0.95,str(ilib+1)+' ARF bins out of '+str(n_elements(rpsfid))+$
          ' complete',/normal,charsi=3.,charth=2
    t1=systime(/seconds)
    totaltime=(t1-t0)*mean(pixperbin)/n_elements(x)*(n_elements(rpsfid)-(ilib+1))
    strtime=str(totaltime/3600./24.,format='(F10.1)')+' days'
    if totaltime/3600./24. lt 2. then $
          strtime=str(totaltime/3600.,format='(F10.1)')+' hr'
    if totaltime/3600. lt 2. then strtime=str(fix(totaltime/60.))+' min'
    if totaltime/60. lt 2. then strtime=str(fix(totaltime))+' sec'
    xyouts,0.02,0.02,'Est. time remaining: '+strtime,/normal,charsi=3.0,charth=2
    counter,ilib+1,n_elements(rpsfid),/percent,'    Status: '
endif

endfor   ; ilib loop

if plotprog then print,'    Completed: '+systime() $
      else print,'Point Source ARF made.'


; write arf & cross-arfs for that src (model)

for r=0,nreg-1 do begin
    sxaddpar,hspec,'LO_THRES',0.0
    arfname=outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf'
    arfstr.specresp=reform(arfs[r,*])
    mwrfits,arfstr,arfname,harf,/silent,/create
    if isrc eq 0 and mkspec then spawn,'fparkey '+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf '+outdir+'/'+ $
          nucrossarf_val.outbase+nucrossarf_obs[iobs].ab+str(r+1)+'.pha ANCRFILE'

; WARNING: if put in same dir, need to add obs # to filename above!!!

    fits_write,outdir+'/'+nucrossarf_val.outbase+nucrossarf_obs[iobs].ab+ $
          '_psfim_src'+str(isrc+1)+'_reg'+str(r+1)+'.fits',$
          reform(regpsfsave[r,*,*]),imhead
endfor

endif else print,'ARFs exist for Source '+str(isrc+1)+ $
      '.  To remake, set clobarf keyword.'

endfor   ; isrc loop


spawn,'rm -f '+cldir+'dummy.fits'


endfor   ; iobs loop


print
print,'ARFs and Cross-ARFs Created'
print



end
