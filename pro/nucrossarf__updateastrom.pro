pro nucrossarf__updateastrom,iobs,nucrossarf_obs,oldheader,newheader

xsh=nucrossarf_obs[iobs].ximsh
ysh=nucrossarf_obs[iobs].yimsh
angle=nucrossarf_obs[iobs].rotang
xrot=nucrossarf_obs[iobs].rotx
yrot=nucrossarf_obs[iobs].roty

if angle gt 1e-3 then begin
    dummyim=intarr(1000,1000)
    hrot,dummyim,oldheader,blah,newheader,angle,xrot,yrot,0
endif else newheader=oldheader

sxaddpar,newheader,'CRPIX1',sxpar(newheader,'CRPIX1')-xsh
sxaddpar,newheader,'CRPIX2',sxpar(newheader,'CRPIX2')-ysh

end
