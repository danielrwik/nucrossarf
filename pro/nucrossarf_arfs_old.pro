pro nucrossarf_arfs,nucrossarf_val,nucrossarf_obs,nucrossarf_src


print
print,'Acquiring ARFs'


nobs=n_elements(nucrossarf_obs)
regnames=nucrossarf_val.regnames
nreg=nucrossarf_val.nreg
if n_elements(regnames) ne nreg then $
      stop,'NUCROSSARF_ARFS: number of regions inconsistent somehow'



; don't forget to multiply psf by exposure map so weighting is correct!

for iobs=0,nobs-1 do begin


ab=nucrossarf_obs[iobs].ab
obsid=nucrossarf_obs[iobs].obsid
cldir=nucrossarf_obs[iobs].obsdir+'/'+obsid+'/event_cl/'
evtfile=cldir+'nu'+obsid+ab+'01_cl.evt'
nucrossarf__header,cldir,evtfile,imhead
fits_write,cldir+'dummy.fits',intarr(1000,1000),imhead

regmask=intarr(nreg,1000,1000)
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

nucrossarf__header,cldir,evtfile,header

tstart=systime(/seconds)
ipt=0
iext=0


xpsf=325
psfsavs=file_search(nucrossarf_val[iobs].libdir+'/psf/'+obsid+'/psf'+ab+'*sav')
if psfsavs[0] eq '' then stop,'NUCROSSARF_ARFS: psf library not found'
undefine,xpsfs,ypsfs
for i=0,n_elements(psfsavs)-1 do begin
    s=strsplit(psfsavs[i],'_',/extract)
    ns=n_elements(s)
    push,xpsfs,fix(s[ns-2])
    ss=strsplit(s[ns-1],'.',/extract)
    push,ypsfs,fix(ss[0])
endfor

arfsavs=file_search(nucrossarf_val[iobs].libdir+'/arf/'+obsid+'/arf'+ab+'*sav')
if arfsavs[0] eq '' then stop,'NUCROSSARF_ARFS: arf library not found'
undefine,xarfs,yarfs
for i=0,n_elements(arfsavs)-1 do begin
    s=strsplit(arfsavs[i],'_',/extract)
    ns=n_elements(s)
    push,xarfs,fix(s[ns-2])
    ss=strsplit(s[ns-1],'.',/extract)
    push,yarfs,fix(ss[0])
endfor


print
print,'Observation '+str(iobs+1)+' out of '+str(nobs)
print,'EVT: '+evtfile

for isrc=0,nmod-1 do begin

undefine,offaxis
if order[isrc] eq 1 then begin
    ra=nucrossarf_src[ipt].ptra
    dec=nucrossarf_src[ipt].ptdec
    scl=1.
    ipt++
endif else if order[isrc] eq 2 then begin
    if nucrossarf_src[isrc].extsrcimg eq 'flat' then begin
        srcim=reg2mask(cldir+'dummy.fits',nucrossarf_val.extdir+'/'+ $
              nucrossarf_src[isrc].extregimg)
        hreg=imhead
        hsrc=imhead
        ii=where(srcim gt 0.5)
        if ii[0] eq -1 then $
              stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+' not found in region: '+$
                    nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg
        ii2d=array_indices(srcim,ii)
    endif else begin
        fits_read,nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extsrcimg,srcim,hsrc
        fits_read,nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg,regim,hreg
        ii=where(regim eq nucrossarf_src[isrc].extregnum)
        if ii[0] eq -1 then $
              stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+ $
                    ' not found in region image '+$
                    nucrossarf_val.extdir+'/'+nucrossarf_src[isrc].extregimg
        ii2d=array_indices(regim,ii)
    endelse
    xyad,hreg,ii2d[0,*],ii2d[1,*],ra,dec
    xyad,hsrc,ii2d[0,0],ii2d[1,0],rach,decch
    if sphdist(rach,decch,ra[0],dec[0],/degrees)*3600. gt 1. then $
          stop,'NUCROSSARF_ARFS: extended source and region image astrometry '+$
          'is inconsistent at >1", should be identical'
    scl=srcim[ii]/total(srcim[ii])
    iext++
endif else stop,'NUCROSSARF_ARFS: src '+str(isrc+1)+' neither pt nor ext???'


adxy,header,ra,dec,xo,yo
x=(xo-499.5)*cos(-nucrossarf_obs[iobs].rotast*!pi/180.) - $
      (yo-499.5)*sin(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
      499.5 + nucrossarf_obs[iobs].xshast
y=(xo-499.5)*sin(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
      (yo-499.5)*cos(-nucrossarf_obs[iobs].rotast*!pi/180.) + $
      499.5 + nucrossarf_obs[iobs].yshast
xs=x
ys=y
scls=scl

print,x,y

;x=round(reform(xs))
;y=round(reform(ys))
;scl=scls


print,'  Source '+str(isrc+1)+' out of '+str(nmod)

arfs=fltarr(nreg,n_elements(earf))
for is=0,n_elements(x)-1 do begin
    blah=min((x[is]-xpsfs)^2+(y[is]-ypsfs)^2,ii)
;    blah=min(abs(y[is]-ypsfs),jj)
    file=nucrossarf_val[iobs].libdir+'psf/'+obsid+'/psf'+ab+'_'+ $
          str(xpsfs[ii])+'_'+str(ypsfs[ii])+'.sav'
print,file
    if not file_test(file) then stop,'NUCROSSARF_ARFS: whaaaaaaaaaaaa? (psf)'
; psflib variable restored by below call
    restore,file


    blah=min((x[is]-xarfs)^2+(y[is]-yarfs)^2,ii)
;    blah=min(abs(x[is]-xarfs),ii)
;    blah=min(abs(y[is]-yarfs),jj)
    file=nucrossarf_val[iobs].libdir+'arf/'+obsid+'/arf'+ab+'_'+ $
          str(xarfs[ii])+'_'+str(yarfs[ii])+'.sav'
    if not file_test(file) then stop,'NUCROSSARF_ARFS: whaaaaaaaaaaaa? (arf)'
; arf variable restored by below call
    restore,file
print,file

;    onearf=nucrossarf__arf(cldir,obsid,ab,x[icent],y[icent],$
;          optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop,pa,offang=oang)
;    arf+=onearf*total(scl[ii])
;    foldspec=refspec*onearf

    xp=fltarr(325,325)
    yp=xp
    for xx=0,324 do begin
        xp[xx,*]=xx
        yp[*,xx]=xx
    endfor
    rpsf=sqrt((xp-162.)^2+(yp-162.)^2)*2.46/60.
    ii=where(rpsf gt 8.0)

    temppsf=dblarr(nepsf,1000,1000)
    pweight=fltarr(nepsf,nreg)
    for p=0,nepsf-1 do begin
        psf=reform(psflib[*,*,p])
;apply psfradcoeff here: create radius image in arcsec, for loop to make polycorr image, multiply psf image by it
;        corr=fltarr(325,325)
;        for f=0,fix(polymax[2,p]) do corr+=psfradcoeff[f,p]*rpsf^f
;        ii=where(rpsf ge (polymax[0,p])*2.46/60.)
;        corr[ii]=polymax[1,p]
;        corr2=fltarr(325,325)
;        for f=0,fix(polymax2[2,p]) do corr2+=psfradcoeff2[f,p]*rpsf^f
;        ii=where(rpsf ge (polymax2[0,p])*2.46/60.)
;        corr2[ii]=polymax2[1,p]
;        jj=where(corr lt 0.5)
;        corr[jj]=1.0
;        jj=where(corr gt 1.5)
;        corr[jj]=1.0
;plot,rpsf,corr,psym=1
;if p eq 5 then stop
;        psf*=corr ;*corr2
        temp2psf=dblarr(1000,1000)
        temp2psf[xpsf/2-162:xpsf/2+162,xpsf/2-162:xpsf/2+162]=psf
        temppsf[p,*,*]=fshift(temp2psf,x[is]-162.,y[is]-162.)
;        for k=0,n_elements(ii)-1 do temppsf+=  $
;             fshift(temp2psf,x[ii[k]]-162.,y[ii[k]]-162.)*  $
;             scl[ii[k]]*total(foldspec[pchan1[p]:pchan2[p]])
;        if not finite(total(temppsf)) then $
;               stop,'NUIMYLZE_OPT: PSF image contains bad values for some reason.'
        temppsf[p,*,*]*=det
        for r=0,nreg-1 do begin
            mask=reform(regmask[r,*,*])
            pweight[p,r]=total(temppsf[p,*,*]*mask)
        endfor
        fits_write,cldir+'psf'+ab+eband[p]+'keV.fits',reform(temppsf[p,*,*]),imhead
    endfor

; for each region, fit pweight to linear function
col=[1,2,3]
    for r=0,nreg-1 do begin
        pw=pweight[*,r]
;        coeff=linfit(epsflog,pw)
;        pscl=coeff[0]+coeff[1]*earf
        pscl=interpol(pw,epsflog,earf)
if r eq 0 then        plot,epsflog,pw,psym=1,col=col[r] else oplot,epsflog,pw,psym=1,col=col[r]
        oplot,earf,pscl,col=col[r]
  ff=fltarr(n_elements(epsflog))
  for t=0,n_elements(epsflog)-1 do begin
   ii=where(earf ge epsf1[t] and earf lt epsf2[t])
   ff[t]=total(pscl[ii])/(epsf2[t]-epsf1[t])*0.04
  endfor
  oplot,epsflog,ff,psym=2,col=col[r]
;        if r eq 0 then arfs[r,*]+=arf*scl*mean(pscl) else arfs[r,*]+=arf*pscl*scl
;if r eq 1 then scl*=0.97
;        arfs[r,*]+=arf*(pscl+(1.-pw[n_elements(pw)-1]))*pw[n_elements(pw)-1]*scl
        arfs[r,*]+=arf*pscl*scl
    endfor

; then scale arf by that and scl factor from dist. of emission
; add scaled arf to final arf for each region for this source

endfor

; write arf & cross-arfs for that src (model)

for r=0,nreg-1 do begin
    sxaddpar,hspec,'LO_THRES',0.0
    arfname=nucrossarf_val.outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf'
    arfstr.specresp=reform(arfs[r,*])
    mwrfits,arfstr,arfname,harf,/silent,/create
    if isrc eq 0 then spawn,'fparkey '+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf '+ $
          nucrossarf_val.outdir+'/'+ $
          nucrossarf_val.outbase+nucrossarf_obs[iobs].ab+str(r+1)+'.pha ANCRFILE'
endfor

endfor

spawn,'rm -f '+cldir+'dummy.fits'

endfor

print
print,'ARFs and Cross-ARFs Created'
print

end
