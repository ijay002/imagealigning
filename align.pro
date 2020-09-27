dir='D:\Analysis\align_sparkmass_github\'; directory path****
sr=read_tiff(dir+'2019_05_10_H_.tif'); file name of the DNA-PAINT image****
ca=rtiffstack(dir+'2019_05_02_BD_avg.tif'); read in averaged Ca image ****
	szc=size(ca,/dim)
	ca=float(ca)
	paddir=[10,10]
shrink=2 ; shrink factor for display****
	szc=size(ca,/dim)
cscale=[0.05,0.05]; xy scale for Ca data (in um) (0.053 rounded)
sscale=[0.005, 0.005]; xy scale for DNA-PAINT data (in um)
svec=[(-300),(-525)] ;shift vectors [x,y] in pixel units (can be -ve or +ve)*** To be entered manually

	exf=cscale[0]/sscale[0] ;scale up factor

;;;;;;;;;;;;; resizing images for display;;;;;;;;

	cal=rebin(reform(ca[*,*]),szc[0]*exf,szc[1]*exf) ; resizing ca to match pixel size of the super-resolution image
	cas=rebin(reform(ca[*,*]),szc[0]*exf/shrink,szc[1]*exf/shrink)
	szc2=size(cas,/dim)

	szs=size(sr,/dim)
	ssr=congrid(sr,round(szs[0]/shrink),round(szs[1]/shrink),cubic=-0.5) ; creating small version of DNA-PAINT image
	szs2=size(ssr,/dim)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


wfit,3,cas;
tvscl,(cas<(prank(cas,90)))-5>0,ch=1
tvscl,unpad(shift(pad(ssr<0.001,0),svec[0],svec[1]),0)>0,ch=2 ; displays shifted version of ssr

;end

;;;;;;;;;;;; saving out overlay;;;;;;;;;

rgb=tvrd(/true) ; captures the overlay image into variable called rgb (note that it is a 2D image variable with 3 dimensions; the first is the channels)
	write_tiff,dir+'overlay.tif',rgb ; saves out the overlay image to be opened externally for visualising


sr2=unpad(shift(pad(sr,paddir),svec[0]*shrink,svec[1]*shrink),0); new variable with super-resolution image aligned to the Ca image.
	write_tiff,dir+'sr2.tif',sr2,/float ; writes out the super-resolution image - to be used for puncta detection

tvscl,unpad(shift(pad(ssr<prank(ssr,99.9),0),svec[0],svec[1]),0),ch=2
tvscl,cas,ch=1
rgb=tvrd(/true)
	write_tiff,dir+'overlay_intensity_unscaled.tif',rgb; writes out the overlay with minimal intensity boosting

;;;;;;;;;;; correlation alignment;;;;;;;;

;ccorr=fltarr(szc2)
fshift=corr_sv(cas-prank(cas,5)>0,(sr2<prank(sr2,95))>0,500,corr=ccorr)

if abs(fshift[0]) le 40 then print,'cross-correlation successful in x, fine shift vectors in x (px)   :', fshift[0] $
else print, 'cross-correlation unreliable, fine shift vectors in x (px)   :', 0
if abs(fshift[1]) le 40 then print,'cross-correlation successful in y, fine shift vectors in y (px)   :', fshift[1] $
else print, 'cross-correlation unreliable, fine shift vectors in y (px)   :', 0


;;;;;;;;;;;;;outputs;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


write_tiff,dir+'cal.tif',cal
wncdf2,dir+'shift_vectors.nc',{arr:svec},/ov ;shift-vectors
wncdf2,dir+'paddir.nc',{arr:paddir},/ov;pad sizes
wncdf2,dir+'fshift.nc',{arr:fshift},/ov; fine shift vectors (from correlation)


end