dir='D:\Analysis\align_sparkmass_github\';****
sr2=read_tiff(dir+'sr2.tif') ; re-aligned DNA-PAINT image from align.pro
svec=rncdf(dir+'shift_vectors.nc') ; shift vectors defined in align.pro
fshift=rncdf(dir+'fshift.nc') ; fine shift vectors (from correlation step)
paddir=rncdf(dir+'paddir.nc');RyR image pad vectors
cscale=[0.05,0.05]; xy scale for Ca data (in um) (0.103 rounded)
sscale=[0.005, 0.005]; xy scale for DNA-PAINT data (in um)
shrink=2 ; shrink factor for display
exf=cscale[0]/sscale[0]


ca=rtiffstack(dir+'2019_05_02_BD_avg.tif'); read in original ca data****
	szc=size(ca,/dim)
	ca=float(ca)

for a=0,1 do if abs(fshift[a]) le 40 then svec[a]=svec[a]+fshift[a]

szc=size(ca,/dim)
szs=size(sr2,/dim)

	exf=cscale[0]/sscale[0] ;scale up factor
	cal=rebin(reform(ca[*,*]),szc[0]*exf,szc[1]*exf) ; resizing ca to match pixel sizes
	cas=rebin(reform(ca[*,*]),szc[0]*exf/shrink,szc[1]*exf/shrink);smaller version
	szc2=size(cas,/dim)
	szc=size(cal,/dim)

	szall=fltarr(2,2) ; array to compare dims of both scaled Ca and RyR image
		szall[0,*]=szc
		szall[1,*]=szs
	szfinal=fltarr(2) ; holds the sizes of the final image
		szfinal[0]=max(szall[*,0])
		szfinal[1]=max(szall[*,1])

;;;;;;;;;;;;;;; RyR coords from DNA-PAINT image;;;;
wfit,1,cas
tvscl,cas

dots=rd_tfile(dir+'2019_05_10_H_pos.txt',2,/convert);read in RyR xy coords****
vj=bytarr(szs)

for u=0, n_elements(dots[1,*])-1 do vj[(dots[1,u]+paddir[1]),(dots[0,u]+paddir[0])]=1
vj=shift(vj,svec[0]*shrink,svec[1]*shrink)
	vj2=fltarr(szfinal)
	vj2[0:szs[0]-1,0:szs[1]-1]=vj
	szs2=size(vj2,/dim)

psf=(bytscl(vgauss([30,30,3],[20.0,20.0,1.0],/fwh)) gt 116)[*,*,1]
vjd=dilate(vj2,psf)
vjs=congrid(vjd,szs2[0]/shrink,szs2[1]/shrink,cubic=-0.5)
sr2s=congrid(sr2,szs2[0]/shrink,szs2[1]/shrink,cubic=-0.5)
tvscl,vjs,ch=3



;;;;;;;;;;;;;; Ca sparks coords from XYspark;;;;;

dotc=rd_tfile(dir+'2019_05_02_BD.txt',10,/convert);read in xy coords ****

fsz=2.5
lst=fltarr(n_elements(dotc[0,*]),9)

for v=0, n_elements(dotc[1,*])-1 do begin

		vc0=fltarr(szc[0],szc[1])
		vc0[(round(dotc[2,v]*exf)),(round(dotc[3,v]*exf))]=dotc[0,v] ;assigning puncta coordinates to pixels - notice that x and y are flipped to correct for Baddeley rotation convention

		psf=(bytscl(vgauss([30,30,3],[25.0,25.0,1.0],/fwh)) gt 116)[*,*,1]
		vc=dilate(vc0 gt 0,psf)

			vc2=fltarr(szfinal)
			vc2[0:szc[0]-1,0:szc[1]-1]=vc0
			szc2=size(vc2,/dim) ; overwriting szc

		vcs=congrid(vc2,szc2[0]/shrink,szc2[1]/shrink,cubic=-0.5) ;

		tv,250*vcs,ch=2
		tvscl,vjs,ch=3

;;;;;;;;;;;;;;;;spark by spark analysis;;;;;;;;;

		lst[v,3]=dotc[4,v]; FWHM measured in microns with xyspark
		lst[v,4]=dotc[5,v] ; F/F0 estimated directly from xyspark
		lst[v,5]=dotc[6,v] ; r^2 estimated directly from xyspark
		lst[v,6]=dotc[9,v] ; spark mass estimated directly from xyspark

		fw=lst[v,3] ; looking up the fwhm
		psfl=(bytscl(vgauss([1200,1200,3],[fw*exf*5,fw*exf*5,1.0],/fwh)) gt 116)[*,*,1]
		field=unpad(dilate((pad(vc2,620)), psfl),620)


		if total(field gt 0) gt 0 then begin
			lst[v,1]=total(vj2[where(field gt 0)])

		;;;;;;;;
			coorde=transpose(getxyz(vj2[where(field gt 0)] gt 0))
			coorde2=transpose(getxyz(vj2[where(dfield gt 0)] gt 0))

		endif

	tvscl,congrid(field gt 0,szc2[0]/shrink,szc2[1]/shrink,cubic=-0.5),ch=2
;
	print,'analysing spark number: ',v+1, ' out of total sparks: ',n_elements(dotc[0,*])

endfor

window,3,xs=1000
plot,(lst[*,1])[where(lst[*,1] gt 0)],(lst[*,6])[where(lst[*,1] gt 0)],title='RyR numbers vs spark mass'$
,psym=5,/xstyle,/ylog, xr=[0,prank((lst[*,1])[where(lst[*,1] gt 0)],99.5)],yr=[1E-1,prank(lst[*,6],85)]$
,xtitle='RyRs', ytitle='Spark Mass (A.U.)'


lst2=fltarr(n_elements(dotc[0,*]),5)
lst2[*,0]=lst[*,1]
lst2[*,1:4]=lst[*,3:6]

write_csv,dir+'sparkmass.csv',transpose(lst2),header=transpose(['local RyR count','FWHM (microns)','F/F0','r^2','spark mass(AU)'])
wncdf2,dir+'sparkmass.nc',{arr:lst},/ov

end