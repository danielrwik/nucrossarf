pro nucrossarf_xcm,nucrossarf_val,nucrossarf_obs,nucrossarf_src

nobs=n_elements(nucrossarf_obs)
nmod=n_elements(nucrossarf_src)
regnames=nucrossarf_val.regnames
nreg=nucrossarf_val.nreg

openw,lun,nucrossarf_val.xcmfile,/get_lun

for iobs=0,nobs-1 do for r=0,nreg-1 do begin
    num=iobs*nreg+r+1
    outdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid
    phaname=outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.pha'
    bgdname=outdir+'/bgd'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.pha'
    printf,lun,'data '+str(num)+':'+str(num)+' '+phaname
    printf,lun,'back '+str(num)+' '+bgdname
endfor

for iobs=0,nobs-1 do for r=0,nreg-1 do for isrc=0,nmod-1 do begin
    num=iobs*nreg+r+1
    outdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid
    respfile=outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(r+1)+'.rmf'
    arffile=outdir+'/'+nucrossarf_val.outbase+ $
          nucrossarf_obs[iobs].ab+str(isrc+1)+'_'+str(r+1)+'.arf'
    printf,lun,'resp '+str(isrc+1)+':'+str(num)+' '+respfile
    printf,lun,'arf '+str(isrc+1)+':'+str(num)+' '+arffile
endfor

; write out models, and done!

for isrc=0,nmod-1 do begin
    printf,lun,'model '+str(isrc+1)+':s'+str(isrc+1)+' '+ $
          nucrossarf_src[isrc].modname+'*const'
    for p=0,nucrossarf_src[isrc].mpar-1 do begin
        if strmid(nucrossarf_src[isrc].p[p],0,1) eq '=' then $
              p0=str((float(nucrossarf_src[isrc].plo[p])+$
                    float(nucrossarf_src[isrc].phi[p]))/2.) $
        else p0=str(nucrossarf_src[isrc].p[p])
        printf,lun,p0+' '+ $
                str(nucrossarf_src[isrc].dp[p])+' '+ $
                str(nucrossarf_src[isrc].plo[p])+' '+ $
                str(nucrossarf_src[isrc].plo[p])+' '+ $
                str(nucrossarf_src[isrc].phi[p])+' '+ $
                str(nucrossarf_src[isrc].phi[p])
    endfor
    printf,lun,'1.0 -1 0. 0. 2. 2.'
endfor
for isrc=0,nmod-1 do for p=0,nucrossarf_src[isrc].mpar-1 do begin
    if strmid(nucrossarf_src[isrc].p[p],0,1) eq '=' then $
          printf,lun,'newpar s'+str(isrc+1)+':'+str(p+1)+ $
                nucrossarf_src[isrc].p[p]
endfor
const='s1:'+str(nucrossarf_src[0].mpar+1)
print
print,'XCM file:'
print,'  Obs 1 cross-calibration constant is '+const
if nobs gt 1 then begin
    for r=1,nreg-1 do printf,lun,'newpar s1:'+ $
          str((nucrossarf_src[0].mpar+1)*(nreg+r+1))+'=s1:'+ $
          str((nucrossarf_src[0].mpar+1)*(nreg+1))
    if nmod gt 1 then for isrc=1,nmod-1 do $
          printf,lun,'newpar s'+str(isrc+1)+':'+str(nucrossarf_src[isrc].mpar+1)+$
                '='+const
    for iobs=1,nobs-1 do begin
        const='s1:'+str((nucrossarf_src[0].mpar+1)*(nreg+1))
        print,'  Obs '+str(iobs+1)+' cross-calibration constant is '+const
        printf,lun,'untie '+const
        for isrc=1,nmod-1 do begin
            thisconst='s'+str(isrc+1)+':'+ $
                  str((nucrossarf_src[isrc].mpar+1)*(nreg+1))
            printf,lun,'newpar '+thisconst+'='+const
            for r=1,nreg-1 do printf,lun,'newpar s'+str(isrc+1)+':'+ $
                  str((nucrossarf_src[isrc].mpar+1)*(nreg+r+1))+'='+thisconst
        endfor
    endfor
endif


printf,lun,'statistic cstat'
printf,lun,'cpd /xs'
printf,lun,'setpl e'
printf,lun,'ignore **:**-3.,30.-**
printf,lun,'setpl reb 10 50'
printf,lun,'pl ld'




free_lun,lun

end
