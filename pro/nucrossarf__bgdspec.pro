; mostly cut-and-pasted from nuskybgd_spec, done to work better with crossarf setup

; NUSKYBGD_SPEC
;
; Read in a parameter file and generate a background spectrum (using fakeit in
; XSPEC) for any arbitrary region.
;
; OPTIONAL KEYWORDS
;
;   paramfile --> If paramdir is set and the parameter filename is the default
;                 (bgdfitparams[A/B].dat), this argument can be omitted.
;                 paramfile must include the full path, or at least the path
;                 relative to your IDL directory -- if set, paramdir is not used
;
;   specname --> Name and path (relative to the event_cl/ directory) of the
;                source spectrum if it does not follow the convention:
;                    event_cl/'+specdir+'/'+specdir+ab+'_sr_g30.pha'
;
;   expfactor --> Factor of the exposure time that the spectrum is simulated for.
;                 The default is 100x.  Smaller values can run into small number
;                 counts issues.  Setting expfactor=1 and /perfect can tell you
;                 the typical number of photons in a given energy range you expect
;                 in your source region.
;
;   perfect --> Suppresses the addition of shot noise to the bgd counts spectrum.
;
;   srcdir --> Relative to event_cl/, the output directory of the bgd spectra
;              if you want a location different from where your source spectrum is.
;
;   fakname --> If you don't want the bgd spectrum to have the same name as
;               the source spectrum (but with 'bgd' prefixed to it), you can
;               specify an alternative name to have 'bgd' prefixed on.
;
;   avgfcxb --> Instead of using the value for the fCXB normalization in the
;               parameter file, use the expected value.  In practice this doesn't
;               work very well, so I would only set this with care.
;
;pro nuskybgd_spec,indir,obsid,srcreg,specdir,bgddir,ab,bgddirref, $
;      paramfile=paramfile,specname=specname,expfactor=expfactor,grxe=grxe,$
;      fakname=fakname,perfect=perfect,srcdir=srcdir,avgfcxb=avgfcxb,$
;      forcermf=forcermf,savexcm=savexcm,rmfwithabs=rmfwithabs,ratio=ratio

pro nucrossarf__bgdspec,iobs,r,mask,nucrossarf_val,nucrossarf_obs,$
      ratiofile=ratiofile

dir=nucrossarf_obs[iobs].obsdir
obsid=nucrossarf_obs[iobs].obsid
ab=nucrossarf_obs[iobs].ab
bgddir=nucrossarf_val.bgddir
bgddirref=bgddir
specdir=nucrossarf_val.outdir+nucrossarf_obs[iobs].obsid+'/'
srcreg=nucrossarf_val.regnames[r]
specname=nucrossarf_val.outbase+ab+str(r+1)+'.pha'

auxildir=getenv('NUSKYBGD_AUXIL')+'/'
caldbdir=getenv('CALDB')+'/'
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
expfactor=100.
ctstat='y'
grxe=0

if size(bgddirref,/type) ne 0 then refdir=dir+obsid+'/event_cl/'+bgddirref+'/' $
      else begin
    if file_test(dir+obsid+'/event_cl/'+bgddir+'/bgdap0'+ab+'.fits') then $
          refdir=dir+obsid+'/event_cl/'+bgddir+'/' else $
        if file_test(dir+obsid+'/event_cl/bgdap0'+ab+'.fits') then $
          refdir=dir+obsid+'/event_cl/' else $
        stop,'Cannot find bgddirref: Please set correctly'
endelse
paramfile=dir+obsid+'/event_cl/'+bgddir+'/bgdfitparams'+ab+'.dat'

print,'Pulling BACKSCAL (assumed to be region area as percentage of image) from:'
print,'  '+specdir+specname
if not file_test(specdir+specname) then begin
    stop,'Spectrum does not exist.'
   ; write bit to automatically create, but a few extra steps to do that
endif

pha=mrdfits(specdir+specname,1,hh,/silent)

livetime=sxpar(hh,'LIVETIME')

;if size(srcreg,/type) eq 7 then $
;      mask=reg2mask(refdir+'bgdap0'+ab+'.fits',specdir+srcreg) $
;; addition of addabs2rmf below means we need the source region
;;  else if size(srcreg,/type) eq 2 then mask=srcreg $
;  else stop,'  NUSKYBGD_SPEC: Source region name/mask ill-defined.'
print, "Created mask"
backscl=total(mask)/1000.^2

detfrac=fltarr(4)
dettot=fltarr(4)
apreg=0.
grxereg=0.
for i=0,3 do begin
    fits_read,refdir+'det'+str(i)+ab+'im.fits',detim
    fits_read,refdir+'bgdap'+str(i)+ab+'.fits',apim
    if grxe then fits_read,refdir+'bgdgrxe'+str(i)+ab+'.fits',grxeim
    detfrac[i]=total(detim*mask)
    dettot[i]=total(detim)
    apreg+=total(apim*mask)
    if grxe then grxereg+=total(grxeim*mask)
endfor
dettotfrac=total(detfrac)/total(dettot)
detwt=detfrac/total(detfrac)
detfrac=detfrac/dettot

;readcol,auxildir+'ratios_lineE.dat',eline,width,/silent
;readcol,auxildir+'ratios_lineE.dat',blah,index1,ebreak,index2,$
;      format='(A,F,F,F)',/silent
;readcol,auxildir+'ratios'+ab+'.dat',f0,f1,f2,f3,/silent
if size(ratiofile,/type) ne 7 then ratiofile=auxildir+'ratios'+ab+'.dat'
readcol,ratiofile,eline,width,f0,f1,f2,f3,/silent
readcol,ratiofile,index1,index2,b0,b1,b2,b3,ebreak,/silent
;neut=[eline[n_elements(f0)-1],width[n_elements(f0)-1]]
eline=eline[0:n_elements(width)-3]
width=width[0:n_elements(width)-3]

readcol,paramfile,p,/silent
apnorm=p[0]*0.002353*apreg/32.
fcxbnorm=p[1]*dettotfrac
;neutnorm=p[2]*dettotfrac
readcol,paramfile,p0,p1,p2,p3,/silent,skipline=3
pinstr=fltarr(n_elements(p0))

for i=0,n_elements(pinstr)-1 do pinstr[i]=total([p0[i],p1[i],p2[i],p3[i]]*detfrac)
if grxe then begin
    readcol,paramfile,blah,gp,format='(A,F)',/silent
    ii=where(blah eq 'GRXE')
    if n_elements(ii) eq 2 and blah[ii[0]] eq 'GRXE' then $
          grxenorm=[gp[ii[0]],gp[ii[1]]]*grxereg $
    else stop,'NUSKYBGD_SPEC: Problem with GRXE values in parameter file'
endif

;namesplit=strsplit(specname,'.',/extract)
;if strmid(namesplit[0],strlen(namesplit[0])-4) eq '_g30' then $
;      rmfname=strmid(namesplit[0],0,strlen(namesplit[0])-4)+'.rmf' $
;      else rmfname=namesplit[0]+'.rmf'
h=headfits(specdir+'/'+specname,exten=1,/silent)
rmfname=sxpar(h,'RESPFILE')

if file_test(specdir+'/temp.rmf') then spawn,'rm -f '+specdir+'/temp.rmf'
;if not keyword_set(rmfwithabs) then $
;      addabs2rmf,cldir+specdir+'/'+rmfname,ab,refdir+'bgdap0'+ab+'.fits',$
;      cldir+srcreg,cldir+bgddir,cldir+specdir+'/temp.rmf',method=2 $
;else spawn,'cp '+cldir+specdir+'/'+rmfname+' '+cldir+specdir+'/temp.rmf'
spawn,'cp '+specdir+'/'+rmfname+' '+specdir+'/temp.rmf'
rmfname='temp.rmf'

fakname=specname
srcdir=specdir

openw,lun,'temp.xcm',/get_lun
printf,lun,'model cutoffpl'
printf,lun,'1.29 -1'
printf,lun,'41.13 -1'
printf,lun,str(apnorm)
spawn,'rm -f '+srcdir+'/bgdap'+fakname
printf,lun,'fakeit none & '+specdir+'/'+rmfname+' & '+auxildir+$
      'be.arf & '+ctstat+' &  & '+srcdir+'/bgdap'+fakname+' & '+$
      str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model cutoffpl'
printf,lun,'1.29 -1'
printf,lun,'41.13 -1'
if keyword_set(avgfcxb) then printf,lun,str(0.002353*(2.45810736/3600.*1000.)^2* $
      backscl) else printf,lun,str(fcxbnorm)
spawn,'rm -f '+srcdir+'/bgdfcxb'+fakname
printf,lun,'fakeit none & '+specdir+'/'+rmfname+' & '+$
;      caldbdir+'data/nustar/fpm/bcf/arf/nu'+ab+'20100101v004.arf & '+$
      auxildir+'fcxb'+ab+'.arf & '+$
      ctstat+' &  & '+srcdir+'/bgdfcxb'+fakname+' & '+$
      str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model ',format='($,A)'
for i=0,n_elements(eline)-1 do printf,lun,'lorentz+',format='($,A)'
printf,lun,'apec'
for i=0,n_elements(eline)-1 do begin
    printf,lun,str(eline[i])+' -1'
    printf,lun,str(width[i])+' -1'
    printf,lun,str(pinstr[i])
endfor
printf,lun,str(index1[0])+' -1'
printf,lun,str(index2[0])+' -1'
printf,lun,str(ebreak[0])+' -1'
printf,lun,str(pinstr[n_elements(eline)])
spawn,'rm -f '+srcdir+'/bgdinstr'+fakname
printf,lun,'fakeit none & '+specdir+'/'+rmfname+' &  & '+ctstat+' &  & '+$
      srcdir+'/bgdinstr'+fakname+' & '+str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model bknpo'
printf,lun,str(index1[1])+' -1'
printf,lun,str(ebreak[1])+' -1'
printf,lun,str(index2[1])+' -1'
printf,lun,str(pinstr[n_elements(eline)+1])
printf,lun,'fakeit none & '+auxildir+'diag.rmf &  & '+ctstat+' &  & '+$
      srcdir+'/bgdintcont'+fakname+' & '+str(livetime*expfactor)
printf,lun,'data none'

if grxe then begin

printf,lun,'model nuabs*(gauss+atable{'+auxildir+'polarmodel.fits})'
printf,lun,'6.7 -1'
printf,lun,'0.0 -1'
printf,lun,str(grxenorm[0])
printf,lun,'0.6 -1'
printf,lun,str(grxenorm[1])
spawn,'rm -f '+srcdir+'/bgdgrxe'+fakname
printf,lun,'fakeit none & '+specdir+'/'+rmfname+' & '+auxildir+$
      'be.arf & '+ctstat+' &  & '+srcdir+'/bgdgrxe'+fakname+' & '+$
      str(livetime*expfactor)

endif

printf,lun,'exit'
free_lun,lun

f = file_info(srcdir+'/bgdintcont'+fakname)
IF f.exists THEN spawn, 'rm -f '+srcdir+'/bgdintcont'+fakname

spawn,'xspec - temp.xcm'
if keyword_set(savexcm) then spawn,'cp -f temp.xcm '+srcdir+'/'+savexcm
spawn,'rm -f temp.xcm'

spawn,'rm -f '+specdir+'/temp.rmf'

spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdap'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdfcxb'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdinstr'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdintcont'+fakname+' BACKSCAL'
if grxe then spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdgrxe'+fakname+' BACKSCAL'

spawn,'pwd > '+srcdir+'/temp'
cd,srcdir
thisdir='' ;specdir+'/'
if not grxe then begin

      print,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
            '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'


      spawn,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
            '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'

      print, expfactor*livetime

      spawn,'mathpha '+$
            thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd_allint_'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'


   ENDIF ELSE  spawn,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
                     '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
                     '+'+thisdir+'bgdgrxe'+fakname+$
                     ' COUNTS '+thisdir+'bgd'+fakname+' '+$
                     str(livetime*expfactor)+' 1.0 0 clobber=yes'
spawn,'fparkey '+str(backscl)+' '+'bgd'+fakname+' BACKSCAL'
print,'fparkey '+str(backscl)+' '+'bgd'+fakname+' BACKSCAL'
spawn,'fparkey bgd'+fakname+' '+specname+' BACKFILE'
print,'fparkey bgd'+fakname+' '+specname+' BACKFILE'



spawn,'fparkey '+str(backscl)+' '+'bgd_allint_'+fakname+' BACKSCAL'
print,'fparkey '+str(backscl)+' '+'bgd_allint_'+fakname+' BACKSCAL'


line=''
openr,lun,'temp',/get_lun
readf,lun,line
free_lun,lun
spawn,'rm -f temp'
cd,line

if file_test(specdir+'temp.rmf') then spawn,'rm -f '+specdir+'temp.rmf'

print,str(expfactor)+'x faked background spectra created:'
print,'  '+srcdir+'bgd'+fakname


end
