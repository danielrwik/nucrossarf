pro nucrossarf__rmf,iobs,r,mask,nucrossarf_val,nucrossarf_obs

cmprmf=1

dir=nucrossarf_obs[iobs].obsdir
obsid=nucrossarf_obs[iobs].obsid
ab=nucrossarf_obs[iobs].ab
bgddir=dir+'/'+obsid+'/event_cl/'+nucrossarf_val.bgddir
bgddirref=bgddir
specdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid+'/'
srcreg=nucrossarf_val.regnames[r]
rmfname=specdir+nucrossarf_val.outbase+ab+str(r+1)

bgddet=fltarr(1000,1000,4)
for i=0,3 do begin
    file=bgddir+'/det'+str(i)+ab+'im.fits'
    if file_test(file) then fits_read,file,im $
          else stop,'NUCROSSARF__RMF: Det image file '+file+' not found.'
    bgddet[*,*,i]=im
endfor

bgdval=fltarr(4)
for i=0,3 do bgdval[i]=total(mask*bgddet[*,*,i])
bgdval/=total(bgdval)
ii=where(bgdval lt 0.01)
bgdval[ii]=0.
ii=where(bgdval ge 0.01)
if ii[0] eq -1 then $
      stop,'NUCROSSARF__RMF: Region does not overlap with any detectors'
bgdval/=total(bgdval)
matrix=fltarr(4096,4096)
rmf1str=mrdfits(getcaldbfile('rmf',ab,0,refmjd=nucrossarf_obs[iobs].mjd),$
      2,r2h,/silent)
for i=0,n_elements(ii)-1 do begin
    absstr=mrdfits(getcaldbfile('detabs',ab,ii[i],refmjd=nucrossarf_obs[iobs].mjd),$
          ii[i]+1,dh,/silent)
    abs=absstr.detabs
    rmfstr=mrdfits(getcaldbfile('rmf',ab,ii[i],refmjd=nucrossarf_obs[iobs].mjd),$
          1,rh,/silent)
    thismatrix=rmfstr.matrix
    for e=0,4095 do thismatrix[*,e]=thismatrix[*,e]*abs[e]
    matrix+=thismatrix*bgdval[ii[i]]
endfor

mwrfits,rmf1str,rmfname+'_temp.rmf',r2h,/silent,/create
rmfstr.matrix=matrix
sxaddpar,rh,'LO_THRES',0.0,'modified by cmprmf'
mwrfits,rmfstr,rmfname+'_temp.rmf',rh,/silent

if cmprmf then spawn,'cmprmf '+rmfname+'_temp.rmf '+ $
      rmfname+'.rmf 1e-6 clobber=yes' $
    else spawn,'cp -f '+rmfname+'_temp.rmf '+rmfname+'.rmf'
print,'RMF '+rmfname+'.rmf created'
spawn,'rm -f '+rmfname+'_temp.rmf'

end
