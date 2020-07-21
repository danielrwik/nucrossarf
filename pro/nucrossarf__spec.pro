pro nucrossarf__spec,iobs,r,mask,ii2d,nucrossarf_val,nucrossarf_obs

cldir=nucrossarf_obs[iobs].obsdir+'/'+nucrossarf_obs[iobs].obsid+'/event_cl/'
evtfile=cldir+'nu'+nucrossarf_obs[iobs].obsid+nucrossarf_obs[iobs].ab+'01_cl.evt'
evts=mrdfits(evtfile,1,header,/silent)

refspechdr=getenv('NUSKYBGD_AUXIL')+'/spechdr.txt'
undefine,hdr
line=''
openr,lun,refspechdr,/get_lun
while ~eof(lun) do begin
    readf,lun,line
    push,hdr,line
endwhile
free_lun,lun

npix=0
specstr=replicate({CHANNEL:0, COUNTS:0, QUALITY:0, GROUPING:1},4096)
specstr.channel=indgen(4096)
for ix=min(ii2d[0,*]),max(ii2d[0,*]) do begin
    chcont=reform(mask[ix,*])
    shchcont=chcont-shift(chcont,1)
    i1=where(shchcont eq 1) 
    i2=where(shchcont eq -1)-1
    if i1[0] eq -1 or i2[0] eq -1 then stop,'NUCROSSARF_DET: invalid region'
    for iy=0,n_elements(i1)-1 do begin
        npix+=(i2[iy]-i1[iy])+1
        ii=where(evts.x-1 eq ix and evts.y-1 ge i1[iy] and evts.y-1 le i2[iy] and $
              evts.grade ge 0 and evts.grade le 26)
        if ii[0] ne -1 then begin
            if n_elements(ii) gt 4096 then begin
                for j=0,4095 do begin
                    jj=where(evts[ii].pi eq j)
                    if jj[0] ne -1 then specstr[j].counts+=n_elements(jj)
                endfor
            endif else for j=0,n_elements(ii)-1 do $
                  specstr[evts[ii[j]].pi].counts++
        endif
    endfor
endfor

if npix ne long(total(mask)) or npix ne n_elements(ii2d[0,*]) then $
      stop,'NUCROSSARF__SPEC: only filtered on '+str(npix)+ $
            ' pixels while the mask for region '+str(r)+' should have '+ $
            str(n_elements(ii2d[0,*]))+' pixels'


grp=3
bin=0
first=1
for i=0,n_elements(specstr)-1 do begin
    if not first then specstr[i].grouping=-1
    bin+=specstr[i].counts
    if bin lt grp then first=0 else begin
        first=1
        bin=0
    endelse
endfor
sxaddpar,hdr,'EXPOSURE',nucrossarf_obs[iobs].exposure
sxaddpar,hdr,'LIVETIME',nucrossarf_obs[iobs].exposure
sxaddpar,hdr,'INSTRUME','FPM'+nucrossarf_obs[iobs].ab
sxaddpar,hdr,'BACKSCAL',total(mask)/1000./1000.
sxaddpar,hdr,'RESPFILE',nucrossarf_val.outbase+nucrossarf_obs[iobs].ab+ $
      str(r+1)+'.rmf'

outdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid

mwrfits,specstr,outdir+'/'+nucrossarf_val.outbase+ $
      nucrossarf_obs[iobs].ab+str(r+1)+'.pha',hdr,/create,/silent

print,'Spectrum '+outdir+'/'+nucrossarf_val.outbase+ $
      nucrossarf_obs[iobs].ab+str(r+1)+'.pha created'

end
