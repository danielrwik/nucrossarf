pro nucrossarf_det,nucrossarf_val,nucrossarf_obs, $
      clobspec=clobspec,clobrmf=clobrmf,clobbgd=clobbgd

if not keyword_set(clobspec) then clobspec=0
if not keyword_set(clobrmf) then clobrmf=0
if not keyword_set(clobbgd) then clobbgd=0
nobs=n_elements(nucrossarf_obs)

;refcldir=nucrossarf_val.refdir+'/'+nucrossarf_val.refobsid+'/event_cl/'
;refevt=refcldir+'nu'+nucrossarf_val.refobsid+'A01_cl.evt'
;if not file_test(refevt) then $
;      stop,'NUIMYLZE_DET: reference event file '+refevt+' not found'
;nucrossarf__header,refcldir,refevt,refimhead
;refpix=[sxpar(refimhead,'CRPIX1'),sxpar(refimhead,'CRPIX2')]
;refpos=[sxpar(refimhead,'CRVAL1'),sxpar(refimhead,'CRVAL2')]
;
;adxy,refimhead,nucrossarf_val.boxcent[0],nucrossarf_val.boxcent[1],refxcent,refycent
;refxcent=round(refxcent)
;refycent=round(refycent)

regnames=nucrossarf_val.regnames
nreg=nucrossarf_val.nreg
if n_elements(regnames) ne nreg then $
      stop,'NUCROSSARF_DET: number of regions inconsistent somehow'

if getenv('NUSKYBGD_AUXIL') eq '' then $
      stop,'NUCROSSARF_DET: environment variable NUSKYBGD_AUXIL not defined.'

for iobs=0,nobs-1 do begin

cldir=nucrossarf_obs[iobs].obsdir+'/'+nucrossarf_obs[iobs].obsid+'/event_cl/'
evtfile=cldir+'nu'+nucrossarf_obs[iobs].obsid+nucrossarf_obs[iobs].ab+'01_cl.evt'

; acquire astrometry shifts for this obs

nucrossarf__header,cldir,evtfile,head
nucrossarf__updateastrom,iobs,nucrossarf_obs,head,imhead
fits_write,cldir+'dummy.fits',intarr(1000,1000),imhead
nucrossarf_obs[iobs].exposure=sxpar(imhead,'EXPOSURE')

;pix=[sxpar(imhead,'CRPIX1'),sxpar(imhead,'CRPIX2')]
;pos=[sxpar(imhead,'CRVAL1'),sxpar(imhead,'CRVAL2')]
;if pix[0] ne refpix[0] and pix[1] ne refpix[1] then $
;      stop,'Reference header issue? ',pos,refpos
;adxy,refimhead,pos[0],pos[1],xref,yref
;nucrossarf_obs[iobs].ximsh=round(refpix[0]-1-xref)+nucrossarf_obs[iobs].ximsh
;nucrossarf_obs[iobs].yimsh=round(refpix[1]-1-yref)+nucrossarf_obs[iobs].yimsh
xysh=[nucrossarf_obs[iobs].ximsh,nucrossarf_obs[iobs].yimsh]
pltscl=6.828076e-4
dx=nucrossarf_val.boxwidth[0]
dy=nucrossarf_val.boxwidth[1]
nx=round(dx/pltscl)
if nx mod 2 eq 0 then nx++
ny=round(dy/pltscl)
if ny mod 2 eq 0 then ny++

;nucrossarf_obs[iobs].xybox=[refxcent-xysh[0]-nx/2,refycent-xysh[1]-ny/2,$
;      refxcent-xysh[0]+nx/2,refycent-xysh[1]+ny/2]

adxy,imhead,nucrossarf_val.boxcent[0],nucrossarf_val.boxcent[1],xcent,ycent
nucrossarf_obs[iobs].xybox=[xcent-xysh[0]-nx/2,ycent-xysh[1]-ny/2,$
      xcent-xysh[0]+nx/2,ycent-xysh[1]+ny/2]

if nucrossarf_val.checkregions then begin
    regmask=intarr(1000,1000)
    for r=0,nreg-1 do begin
        mask=reg2mask(cldir+'dummy.fits',regnames[r])
        regmask[*,*]+=mask
    endfor
    if max(mask) gt 1 then $
     stop,'NUCROSSARF_DET: regions overlap, events will be included more than once'+$
            string(10B)+'        set param noregch to 1 to proceed anyway,'+$
            string(10B)+'        which I think would be a mistake, statistically'+$
            srting(10B)+'        speaking, but hey, it is your life.'
endif

; get regions and extract products

for r=0,nreg-1 do begin

mask=reg2mask(cldir+'dummy.fits',regnames[r])
;mask=shift(mask,nucrossarf_obs[iobs].ximsh,nucrossarf_obs[iobs].yimsh)
ii=where(mask gt 0.5)
ii2d=array_indices(intarr(1000,1000),ii)

outdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid
mkspec=nucrossarf_val.mkspec

; extract spectrum

if mkspec then begin 
    if file_test(outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.pha') and clobspec eq 0 then begin
        print,'Spectrum already exists, set clobspec=1 to overwrite:'
        print,'     '+outdir+'/'+nucrossarf_val.outbase+ $
              nucrossarf_obs[iobs].ab+str(r+1)+'.pha'
    endif else $
      nucrossarf__spec,iobs,r,mask,ii2d,nucrossarf_val,nucrossarf_obs
endif else print,'SPEC for region '+str(r+1)+' skipped, set mkspec to 1 create

; create rmf

if mkspec then begin
    if file_test(outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.rmf') $
          and clobrmf eq 0 then begin
        print,'Response already exists, set clobrmf=1 to overwrite:'
        print,'     '+outdir+'/'+nucrossarf_val.outbase+ $
              nucrossarf_obs[iobs].ab+str(r+1)+'.rmf'
    endif else $
      nucrossarf__rmf,iobs,r,mask,nucrossarf_val,nucrossarf_obs
endif else print,'RMF for region '+str(r+1)+' skipped, set mkspec to 1 create


; generate bgd

if mkspec then begin
    undefine,ratiofile
    if nucrossarf_obs[iobs].bgdratiofile ne '' then $
          ratiofile=nucrossarf_obs[iobs].bgdratiofile
    if file_test(outdir+'/bgd'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.pha') $
          and clobbgd eq 0 then begin
        print,'Background already exists, set clobbgd=1 to overwrite:'
        print,'     '+outdir+'/bgd'+nucrossarf_val.outbase+ $
              nucrossarf_obs[iobs].ab+str(r+1)+'.pha'
    endif else $
      nucrossarf__bgdspec,iobs,r,mask,nucrossarf_val,nucrossarf_obs,$
            ratiofile=ratiofile
endif else print,'BGD for region '+str(r+1)+' skipped, set mkspec to 1 create


endfor


; make psf and arf libraries

nucrossarf__libraries,iobs,nucrossarf_val,nucrossarf_obs


endfor

spawn,'rm -f '+cldir+'dummy.fits'

end
