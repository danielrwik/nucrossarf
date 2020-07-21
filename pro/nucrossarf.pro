pro nucrossarf,dir,base,savpsf=savpsf,savbin=savbin,quickshift=quickshift, $
      clobspec=clobspec,clobrmf=clobrmf,clobbgd=clobbgd,clobarf=clobarf, $
      checkast=checkast


; read in input files, load structures

nucrossarf_init,dir,base,nucrossarf_val,nucrossarf_obs,nucrossarf_src



; make output directories if they don't already exist

if not file_test(nucrossarf_val.outdir,/directory) then $
      spawn,'mkdir '+nucrossarf_val.outdir

nobs=n_elements(nucrossarf_obs)
for iobs=0,nobs-1 do if not file_test(nucrossarf_val.outdir+ $
      nucrossarf_obs[iobs].obsid,/directory) then $
            spawn,'mkdir '+nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid


; if checkast set, write image files with new astrometry and stop

if keyword_set(checkast) then begin
    if size(checkast,/type) eq 7 then outdir=checkast+'/' $
          else outdir=nucrossarf_val.outdir
    print,'Astrometry-corrected, 3-30 keV images written to '+outdir
    for iobs=0,nobs-1 do begin
        ab=nucrossarf_obs[iobs].ab
        obsid=nucrossarf_obs[iobs].obsid
        cldir=nucrossarf_obs[iobs].obsdir+'/'+obsid+'/event_cl/'
        outprint='checkast_obs'+str(iobs+1)+'.fits'
        outimg=outdir+outprint
        mkimgs,cldir,obsid,ab,3,30,outname=outimg
        fits_read,outimg,im,head
        nucrossarf_updateastrom,iobs,nucrossarf_obs,head,imhead
        fits_write,outimg,im,imhead
        print,'    '+outprint
    endfor
    stop,'Unset checkast keyword to proceed with spec/bgd/rmf/arf generation.'
endif



; read data, bgd, exp images and apply det. absorption/rmf 

nucrossarf_det,nucrossarf_val,nucrossarf_obs,$
      clobspec=clobspec,clobrmf=clobrmf,clobbgd=clobbgd



; create arf-modulated PSF images of each source in each energy band

nucrossarf_arfs,nucrossarf_val,nucrossarf_obs,nucrossarf_src,clobber=clobarf



;;;; would be a good place to save structures for later loading/fitting ;;;;

; write xcm file for xspec fitting

nucrossarf_xcm,nucrossarf_val,nucrossarf_obs,nucrossarf_src



end
