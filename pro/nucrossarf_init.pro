pro nucrossarf_init,dir,base,$
      nucrossarf_val,nucrossarf_obs,nucrossarf_src,$
      altbase=altbase,datapath=datapath,bgddir=bgddir

if keyword_set(datapath) then dir=datapath+dir else datapath=''
if keyword_set(altbase) then inbase=altbase else inbase=base
if not keyword_set(bgddir) then bgddir='bgd'

obsfile=dir+'/'+inbase+'_obs.dat'
if not file_test(obsfile) then stop,'NUCROSSARF_INIT: File not found: '+obsfile
obsparfile=dir+'/'+inbase+'_par.dat'
if not file_test(obsparfile) then $
      stop,'NUCROSSARF_INIT: File not found: '+obsparfile
srcsfile=dir+'/'+inbase+'_srcs.dat'
if not file_test(srcsfile) then stop,'NUCROSSARF_INIT: File not found: '+srcsfile



l=''
nobs=0
undefine,obsdir,obsid,ab,ximsh,yimsh,rotang,rotx,roty
openr,lun,obsfile,/get_lun
while ~eof(lun) do begin
    readf,lun,l
    ent=strsplit(l,' ',/extract)
    if str(ent[0]) ne '#' then begin
        if n_elements(ent) ne 3 and n_elements(ent) ne 5 and $
              n_elements(ent) ne 8 then $
              stop,'NUCROSSARF_INIT: '+obsfile+' contains an invalid line'
        push,obsdir,ent[0]
        push,obsid,ent[1]
        push,ab,ent[2]
        if n_elements(ent) eq 5 or n_elements(ent) eq 8 then begin
            push,ximsh,float(ent[3])
            push,yimsh,float(ent[4])
            if n_elements(ent) eq 8 then begin
                push,rotang,float(ent[5])
                push,rotx,float(ent[6])
                push,roty,float(ent[7])
            endif else begin
                push,rotang,0.0
                push,rotx,0.0
                push,roty,0.0
            endelse
        endif else begin
            push,ximsh,0.0
            push,yimsh,0.0
            push,rotang,0.0
            push,rotx,0.0
            push,roty,0.0
        endelse
        nobs++
    endif
endwhile
free_lun,lun

if keyword_set(datapath) then obsdir=datapath+'/'+obsdir
cldir=obsdir+'/'+obsid+'/event_cl/'
mjuldate=fltarr(nobs)
for i=0,nobs-1 do begin
    if not file_test(cldir[i]+'nu'+obsid[i]+ab[i]+'01_cl.evt') then $
          stop,'NUCROSSARF_INIT: File not found: '+$
                cldir[i]+'nu'+obsid[i]+ab[i]+'01_cl.evt'
    head=headfits(cldir[i]+'nu'+obsid[i]+ab[i]+'01_cl.evt')
    obsstr=sxpar(head,'DATE-OBS')
    t=strsplit(obsstr,'T',/extract)
    tt=strsplit(t[0],'-',/extract)
    ttt=strsplit(t[1],':',/extract)
    juldate,[fix(tt),fix(ttt)],mjd
    mjuldate[i]=mjd
    if not file_test(cldir[i]+bgddir,/directory) then $
          stop,'NUCROSSARF_INIT: Directory not found: '+cldir[i]+bgddir
endfor

header=headfits(cldir[0]+'nu'+obsid[0]+'A01_cl.evt',exten=1,/silent)
radef=sxpar(header,'RA_OBJ')
decdef=sxpar(header,'DEC_OBJ')




pnames=['refobsid','obsincl','obsexcl',$
      'ftol','dftol','nmax','outdir','outbase',$
      'e1','e2','libbox','exclregimg','refdir',$
      'bknpo','extdir','combine',$
      'psh','fitflag','perr','deltac','regbase','nreg','libdir','outxcm',$
      'bgdratio','noregch','srcbin','mkspec']
pdefault=[obsid[0],'','',$  ;2
      '1e-7','1e-8','1000','nucrossarf','nucrossarf',$  ;7
      '4','25',str(radef)+','+str(decdef)+',60.0",60.0"','',obsdir[0],$   ;12
      '2.0,7.0,2.0','','obs0en0',$  ;15
      '','1','all','1.0','reg','2','nope','',$  ;23
      '','0','0','0']   ;27






pvals=pdefault
readcol,obsparfile,parname,value,format='(A,A)',/silent,delimiter=' '
for i=0,n_elements(pnames)-1 do begin
    n=strlen(pnames[i])
    ii=where(strcmp(strmid(str(parname),0,n),pnames[i],/fold_case))
    if ii[0] ne -1 then pvals[i]=value[ii[0]]
    if n_elements(ii) gt 1 then $
          print,'NUCROSSARF_INIT: WARNING: multiple entries found for '+$
                pnames[i]+', first entry used.'
endfor
if pvals[20] ne '' then begin
    regcore=pvals[20]
    ok=execute('nreg='+pvals[21])
    if not ok then stop,'NUCROSSARF_INIT: '+pnames[21]+' syntax invalid'
endif
if pvals[1] ne '' then begin
    ok=execute('iiobs='+pvals[1])
    if not ok then stop,'NUCROSSARF_INIT: '+pnames[1]+' syntax invalid'
endif
if size(iiobs,/type) eq 0 then iiobs=indgen(nobs)
if pvals[2] ne '' then begin
    ok=execute('iiexcl=['+pvals[2]+']')
    if not ok then stop,'NUCROSSARF_INIT: '+pnames[2]+' syntax invalid'
    iiok=replicate(1,nobs)
    iiok[iiexcl]=0
    ii=where(iiok eq 1)
    if ii[0] ne -1 then iiobs=iiobs[ii] else $
          stop,'NUCROSSARF_INIT: after include/exclude selections, '+$
                'no valid observations found'
endif
if pvals[17] eq '1' or strcmp(pvals[17],'yes',/fold_case) or $
      strcmp(pvals[17],'y',/fold_case) then fitflag=1 else fitflag=0
if fitflag eq 0 or pvals[18] eq '0' or strcmp(pvals[18],'no',/fold_case) or $
      strcmp(pvals[18],'n',/fold_case) then errflag=0 else errflag=1
if n_elements(iiobs) ne nobs then begin
    obsdir=obsdir[iiobs]
    obsid=obsid[iiobs]
    ab=ab[iiobs]
    ximsh=ximsh[iiobs]
    yimsh=yimsh[iiobs]
    rotang=rotang[iiobs]
    rotx=rotx[iiobs]
    roty=roty[iiobs]
    nobs=n_elements(iiobs)
endif

undefine,e1,e2
e1str=strsplit(pvals[8],',',/extract)
neng=n_elements(e1str)
e2str=strsplit(pvals[9],',',/extract)
if neng ne n_elements(e1str) then $
      stop,'NUCROSSARF_INIT: upper and lower energy ranges invalid'
e1=fltarr(neng)
e2=fltarr(neng)
for i=0,neng-1 do begin
    e1[i]=float(e1str[i])
    e2[i]=float(e2str[i])
endfor

obscomb=intarr(nobs,neng)
if pvals[15] ne 'obs0en0' and n_elements(obscomb) gt 1 then begin
    spl=strsplit(pvals[15],';',/extract)
    for j=0,n_elements(spl)-1 do begin
        splen=strsplit(spl[j],'en',/extract)
        if n_elements(splen) eq 2 then if splen[1] ne '0' then $
              stop,"NUCROSSARF_INIT: cannot yet deal with combining energy bands :'("
        if n_elements(splen) eq 1 then $
              stop,'NUCROSSARF_INIT: obscomb syntax incorrect'
        span=strmid(splen[0],3)
        if span eq '0' then for i=0,neng-1 do obscomb[*,i]=i else begin
            dash=strsplit(span,'-',/extract)
            comma=strsplit(span,',',/extract)
            if n_elements(dash) gt 1 then begin
                range=fix(dash)-1
                for i=0,neng-1 do obscomb[range[0]:range[1],i]=i+j*neng
            endif else if n_elements(comma) gt 1 then begin
                ids=fix(comma)-1
                for k=0,n_elements(ids)-1 do for i=0,neng-1 do $
                      obscomb[ids[k],i]=i+j*neng
            endif else for i=0,neng-1 do obscomb[fix(span)-1,i]=i+j*neng
        endelse
    endfor
endif else for i=0,neng-1 do obscomb[*,i]=i
nfit=max(obscomb)+1
;for i=0,nobs-1 do print,i,obscomb[i,0]

if pvals[16] ne '' then begin
    spl=strsplit(pvals[16],',',/extract)
    if n_elements(spl) mod 3 ne 0 or n_elements(spl)/3 ne nfit then $
          stop,'NUCROSSARF_INIT: param '+$
          pnames[16]+' must be a comma separated list, multiple of 3 values'
    psh=float(spl)
endif else psh=-999

boxstr=strsplit(pvals[10],',',/extract)
if n_elements(boxstr) ne 4 then $
      top,'NUCROSSARF_INIT: fitting box requires 4 inputs (ra,dec,dx,dy)'
boxcent=[float(boxstr[0]),float(boxstr[1])]
if strmid(str(boxstr[2]),strlen(boxstr[2])-1) eq "'" then $
      dx=float(strmid(str(boxstr[2]),0,strlen(boxstr[2])-1))/60. $
      else if strmid(str(boxstr[2]),strlen(boxstr[2])-1) eq '"' then $
      dx=float(strmid(str(boxstr[2]),0,strlen(boxstr[2])-1))/3600. $
      else dx=float(strmid(str(boxstr[2]),0,strlen(boxstr[2])))
if strmid(str(boxstr[3]),strlen(boxstr[3])-1) eq "'" then $
      dy=float(strmid(str(boxstr[3]),0,strlen(boxstr[3])-1))/60. $
      else if strmid(str(boxstr[3]),strlen(boxstr[3])-1) eq '"' then $
      dy=float(strmid(str(boxstr[3]),0,strlen(boxstr[3])-1))/3600. $
      else dy=float(strmid(str(boxstr[3]),0,strlen(boxstr[3])))
boxwidth=[dx,dy]
pltscl=6.828076e-4
nx=round(dx/pltscl)
if nx mod 2 eq 0 then nx++
ny=round(dy/pltscl)
if ny mod 2 eq 0 then ny++
splpo=strsplit(pvals[13],',',/extract)
if n_elements(splpo) ne 3 then $
      stop,'NUCROSSARF: wrong #params for broken power law model'
phindex=float(splpo)
elow=findgen(4096)*0.04+1.6
ehigh=findgen(4096)*0.04+1.64
ii=where(ehigh le phindex[1])
spec=fltarr(4096)
if phindex[0] eq 1.0 then spec[ii]=(alog(ehigh[ii])-alog(elow[ii])) $
      else spec[ii]=(ehigh[ii]^(1.-phindex[0])-elow[ii]^(1.-phindex[0]))/ $
            (1.-phindex[0])
ii=where(elow ge phindex[1])
if ii[0] ne -1 then if phindex[2] eq 1.0 then $
      spec[ii]=(alog(ehigh[ii])-alog(elow[ii])) $
  else spec[ii]=(ehigh[ii]^(1.-phindex[2])-elow[ii]^(1.-phindex[2]))/ $
        (1.-phindex[2])*phindex[1]^(phindex[2]-phindex[0])

if pvals[24] ne '' then begin
    files=strsplit(pvals[24],',',/extract)
    if n_elements(files) ne nobs then $
          stop,'NUCROSSARF_INIT: comma-delimited bgdratio list must equal # of obs'
    bgdratiofiles=files
endif else bgdratiofiles=replicate('',nobs)

if pvals[26] ne '0' then begin
    if fix(pvals[26]) le 1 then $
          stop,'NUCROSSARF_INIT: default ARF binning is 1 (equal to source image).'+$
   string(10B)+'                 Set srcbin > 1 (integer amounts) to speed up calc.'
    deltaarf=fix(pvals[26])
endif else deltaarf=0
if pvals[27] eq '1' or strcmp(pvals[27],'yes',/fold_case) or $
      strcmp(pvals[27],'y',/fold_case) then mkspec=1 else mkspec=0


; read in source positions / images (diffuse)

o=0
l=''
undefine,srcnum,ptra,ptdec,ptorder,extsrcimg,extregimg,extregnum,extorder
undefine,modname,mpar,p1,dp1,plo1,phi1,ptie1,pfix1,bgdorder
openr,lun,srcsfile,/get_lun
while ~eof(lun) do begin
    readf,lun,l
    spl=strsplit(l,' ',/extract)
    if strmid(str(spl[0]),0,1) ne '#' then begin
        push,srcnum,fix(spl[0])
        if strcmp(strmid(str(spl[1]),0,2),'pt',/fold_case) then begin
            push,ptra,float(spl[2])
            push,ptdec,float(spl[3])
            push,ptorder,o
            push,bgdorder,-999
            push,extorder,-999
            push,extsrcimg,''
            push,extregimg,''
            push,extregnum,0
        endif else if strcmp(strmid(str(spl[1]),0,3),'bgd',/fold_case) then begin
            push,ptra,0.0
            push,ptdec,0.0
            push,bgdorder,o
            push,ptorder,-999
            push,extorder,-999
            push,extsrcimg,''
            push,extregimg,''
            push,extregnum,0
        endif else begin
            push,extsrcimg,spl[1]
            push,extregimg,spl[2]
            push,extregnum,fix(spl[3])
            if spl[1] ne 'flat' then begin
              if not file_test(pvals[14]+'/'+extsrcimg[-1]) then $
                  stop,'NUCROSSARF_INIT: file not found '+pvals[14]+'/'+extsrcimg[-1]
              if not file_test(pvals[14]+'/'+extregimg[-1]) then $
                  stop,'NUCROSSARF_INIT: file not found '+pvals[14]+'/'+extregimg[-1]
            endif
            push,extorder,o
            push,ptorder,-999
            push,bgdorder,-999
            push,ptra,0.0
            push,ptdec,0.0
        endelse
      if n_elements(spl) gt 4 then begin
        splm=strsplit(spl[4],':',/extract)
        if n_elements(splm) ne 2 then $
              stop,'NUCROSSARF_INIT: model names must be of the form mname:npar'
        push,modname,splm[0]
        push,mpar,fix(splm[1])
        if n_elements(spl) ne mpar[o]+5 then $
              stop,'NUCROSSARF_INIT: # parameters ne specified # for model '+modname[o]
        for i=0,mpar[o]-1 do begin
            splp=strsplit(spl[5+i],',',/extract)
            if strmid(splp[0],0,1) eq 'f' then begin
                push,p1,float(strmid(splp[0],1))
                push,dp1,0.
                push,plo1,0.
                push,phi1,0.
                push,pfix1,1
                push,ptie1,''
;            endif else if strmid(splp[0],0,1) eq '=' then begin
;                push,p1,0.
;                push,dp1,0.
;                push,plo1,0.
;                push,phi1,0.
;                push,pfix1,2
;                push,ptie1,strmid(splp[0],1)
            endif else begin
                if n_elements(splp) ne 4 then $
                      stop,'NUCROSSARF_INIT: list parameters as val,delta_val,min,max'
                push,p1,splp[0]
                push,dp1,splp[1]
                push,plo1,splp[2]
                push,phi1,splp[3]
                push,pfix1,0
                push,ptie1,''
            endelse
        endfor
      endif
        o++
    endif
endwhile
free_lun,lun

;if max(bgdorder) lt 0 then for i=0,max(obscomb) do begin
;    push,ptra,0.0
;    push,ptdec,0.0
;    push,bgdorder,o
;    push,ptorder,-999
;    push,extorder,-999
;    push,extsrcimg,''
;    push,extregimg,''
;    push,extregnum,0
;    push,srcnum,max(srcnum)+1
;    push,modname,'bgd'
;    push,mpar,1
;    push,p1,1.0
;    push,dp1,0.
;    push,plo1,0.
;    push,phi1,0.
;    push,pfix1,1
;    push,ptie1,''
;    o++
;endfor
ii=where(bgdorder ge 0)
if n_elements(ii) ne max(obscomb)+1 then $
      stop,'NUCROSSARF_INIT: # bgd models must equal # non-combined obs'
if size(ptra,/type) ne 0 then nptsrc=n_elements(ptra) else nptsrc=0
if size(extsrcimg,/type) ne 0 then nextsrc=n_elements(extsrcimg) else nextsrc=0
nmod=o

maxpar=max(mpar)
parnames=strarr(nmod,maxpar)
ptie=strarr(nmod,maxpar)
pfix=intarr(nmod,maxpar)
p=strarr(nmod,maxpar)
dp=strarr(nmod,maxpar)
plo=strarr(nmod,maxpar)
phi=strarr(nmod,maxpar)
perr=intarr(nmod,maxpar)

np=0
for i=0,nmod-1 do begin
    for j=0,mpar[i]-1 do begin
        parnames[i,j]='s'+str(srcnum[i])+'p'+str(j+1)
        ptie[i,j]=ptie1[j+np]
        pfix[i,j]=pfix1[j+np]
        p[i,j]=p1[j+np]
        dp[i,j]=dp1[j+np]
        plo[i,j]=plo1[j+np]
        phi[i,j]=phi1[j+np]
    endfor
    np+=mpar[i]
endfor

;;if errflag and (pvals[18] eq '0' or strcmp(pvals[18],'all',/fold_case)) then $
;;      perr=pfix $
;;  else if errflag then begin
;;    perr[*,*]=1
;;    splt=strsplit(pvals[18],',',/extract)
;;    for i=0,n_elements(splt)-1 do begin
;;        if strmid(splt[i],0,1) ne 's' then $
;;              stop,'NUCROSSARF_INIT: invalid perr name: '+splt[i]
;;        splt2=strsplit(splt[i],'p',/extract)
;;        if n_elements(splt2) ne 2 then $
;;              stop,'NUCROSSARF_INIT: invalid perr name: '+splt[i]
;;        snum=fix(strmid(splt2[0],1))
;;        ii=where(snum eq srcnum)
;;        if ii[0] eq -1 then $
;;              stop,'NUCROSSARF_INIT: '+splt[i]+ $
;;                    ' references missing source #: '+str(snum)
;;        if n_elements(ii) gt 1 then $
;;              stop,'NUCROSSARF_INIT: '+splt[i]+' refers to multiple potential srcs'
;;        pnum=fix(splt2[1])
;;        if pnum ge maxpar then stop,'NUCROSSARF_INIT: invalid param #: '+splt[i]
;;        perr[ii[0],pnum-1]=0
;;    endfor
;;endif
;;
;;for i=0,nmod-1 do for j=0,mpar[i]-1 do if ptie[i,j] ne '' then begin
;;    splt=strsplit(ptie[i,j],'.',/extract)
;;    expr=''
;;    for k=0,n_elements(splt)-1 do if strmid(splt[k],0,1) eq 's' then begin
;;        ii=where(parnames eq splt[k])
;;        if ii[0] eq -1 or n_elements(ii) gt 1 then $
;;              stop,'NUCROSSARF_INIT: '+splt[k]+' is an invalid parameter name'
;;        pindex=strsplit(splt[k],'p')
;;        snum=fix(strmid(splt[k],1,pindex[1]-2))
;;        ii=where(snum eq srcnum)
;;        if ii[0] eq -1 then $
;;              stop,'NUCROSSARF_INIT: '+splt[k]+ $
;;                    ' references missing source #: '+str(snum)
;;        if n_elements(ii) gt 1 then $
;;              stop,'NUCROSSARF_INIT: '+splt[k]+' refers to multiple potential srcs'
;;        pnum=fix(strmid(splt[k],pindex[1]))-1
;;        expr+='p['+str(ii[0])+','+str(pnum)+']'
;;    endif else begin
;;        if k eq 0 and strmid(ptie[i,j],0,1) eq '.' then tstr='.' else tstr=''
;;        expr+=tstr+splt[k]
;;    endelse
;;    ptie[i,j]=expr
;;endif


; define structures and populate basic fields

nucrossarf_val={refdir:'',refobsid:'',e1:fltarr(neng),e2:fltarr(neng),bgddir:'',$
      outdir:'',outbase:'',ftol:0.0,dftol:0.0,nmax:0,$
      boxcent:fltarr(2),boxwidth:fltarr(2),exclregimg:'',$
      refspec:fltarr(4096),extdir:'',fitflag:1,errflag:0,deltac:0.0,$
      regnames:strarr(nreg),regcore:'',nreg:0,libdir:'',xcmfile:'',checkregions:1,$
      darf:0,mkspec:0}
nucrossarf_val.refobsid=pvals[0]
nucrossarf_val.refdir=pvals[12]
nucrossarf_val.bgddir=bgddir
nucrossarf_val.e1=e1
nucrossarf_val.e2=e2
nucrossarf_val.ftol=float(pvals[3])
nucrossarf_val.dftol=float(pvals[4])
nucrossarf_val.nmax=fix(pvals[5])
if keyword_set(datapath) then nucrossarf_val.outdir=datapath+'/'+pvals[6]+'/' $
      else nucrossarf_val.outdir=pvals[6]+'/'
if pvals[22] eq 'nope' then pvals[22]=nucrossarf_val.outdir else pvals[22]+='/'
nucrossarf_val.libdir=pvals[22]
nucrossarf_val.outbase=pvals[7]
nucrossarf_val.boxcent=boxcent
nucrossarf_val.boxwidth=boxwidth
nucrossarf_val.exclregimg=pvals[11]
nucrossarf_val.refspec=spec
nucrossarf_val.extdir=pvals[14]
nucrossarf_val.fitflag=fitflag
nucrossarf_val.errflag=errflag
nucrossarf_val.regcore=regcore
nucrossarf_val.nreg=nreg
nucrossarf_val.regnames=regcore+str(indgen(nreg)+1)+'.reg'
nucrossarf_val.deltac=float(pvals[19])
if pvals[23] eq '' then pvals[23]=nucrossarf_val.outdir+nucrossarf_val.outbase+'.xcm'
nucrossarf_val.xcmfile=pvals[23]
if pvals[25] eq '1' or pvals[25] eq 'y' or pvals[25] eq 'Y' or pvals[25] eq 'yes' $
      then nucrossarf_val.checkregions=0 else nucrossarf_val.checkregions=1
nucrossarf_val.darf=deltaarf
nucrossarf_val.mkspec=mkspec

nucrossarf_obs=replicate({dim:intarr(neng,nx,ny),eim:fltarr(neng,nx,ny), $
      model:dblarr(neng,nmod,1000,1000),obsdir:'',obsid:'',ab:'',$
      ximsh:0.0,yimsh:0.0,rotang:0.0,rotx:0.0,roty:0.0,xybox:fltarr(4),$
      exposure:0.0,combindex:intarr(neng),mjd:0.0,bgdratiofile:''},nobs)
nucrossarf_obs.obsdir=obsdir
nucrossarf_obs.obsid=obsid
nucrossarf_obs.ab=ab
nucrossarf_obs.ximsh=ximsh
nucrossarf_obs.yimsh=yimsh
nucrossarf_obs.rotang=rotang
nucrossarf_obs.rotx=rotx
nucrossarf_obs.roty=roty
nucrossarf_obs.mjd=mjuldate
nucrossarf_obs.bgdratiofile=bgdratiofiles
for i=0,nobs-1 do nucrossarf_obs[i].combindex=obscomb[i,*]

nucrossarf_src=replicate({ptra:0.0,ptdec:0.0,extsrcimg:'',extregimg:'',extregnum:0,$
      nptsrc:0,nextsrc:0,ptorder:0,extorder:0,srcnum:0,pnames:strarr(maxpar),$
      p:strarr(maxpar),dp:strarr(maxpar),plo:strarr(maxpar),phi:strarr(maxpar),$
      pfix:intarr(maxpar),ptie:strarr(maxpar),modname:'',mpar:0,bgdorder:0,$
      perr:intarr(maxpar)},nmod)
nucrossarf_src.ptra=ptra
nucrossarf_src.ptdec=ptdec
nucrossarf_src.extsrcimg=extsrcimg
nucrossarf_src.extregimg=extregimg
nucrossarf_src.extregnum=extregnum
nucrossarf_src.srcnum=srcnum
nucrossarf_src[*].nptsrc=nptsrc
nucrossarf_src[*].nextsrc=nextsrc
nucrossarf_src.ptorder=ptorder
nucrossarf_src.extorder=extorder
nucrossarf_src.bgdorder=bgdorder
nucrossarf_src.modname=modname
;;nucrossarf_src.modname=''
nucrossarf_src.mpar=mpar
;;nucrossarf_src.mpar=0
for i=0,nmod-1 do begin
    nucrossarf_src[i].pnames=parnames[i,*]
    nucrossarf_src[i].p=p[i,*]
    nucrossarf_src[i].dp=dp[i,*]
    nucrossarf_src[i].plo=plo[i,*]
    nucrossarf_src[i].phi=phi[i,*]
    nucrossarf_src[i].pfix=pfix[i,*]
    nucrossarf_src[i].perr=perr[i,*]
    nucrossarf_src[i].ptie=ptie[i,*]
endfor





end
