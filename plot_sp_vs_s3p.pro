FUNCTION CIRCLE, xcenter, ycenter, radius
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + radius * COS(points )
   y = ycenter + radius * SIN(points )
   RETURN, TRANSPOSE([[x],[y]])
END

;-----------------------------------------------------------
pro read_output,filetype,spec,nt,nlng,nr,lng,val,rad
;-----------------------------------------------------------

lng = 0.0
val = 0.0
close,1
fil = '/Volumes/Scratch/2D_Model/run_vmass/'+filetype+spec+nt+'_3D.dat'
;fil = './plots/'+filetype+spec+nt+'_3D.dat'
;fil = './'+filetype+spec+nt+'_3D.dat'
print,fil
openr,1,fil
;readf,1,lng,val,rad
lng = 0.
val = 0.
rad = 0.
for j = 0,nr-1 do begin
;while not(eof(1)) do begin
   nval = 0.
   for i = 0,nlng-1 do begin
      readf,1,l,v,r
      ;print,l,r
      lng = [lng,l]
      nval = [nval,v]
      ;val = [val,v]
      rad = [rad,r]
   endfor
   nval = nval(1:*)
   minnv = min(nval)
   maxnv = max(nval)
   nval = -1 + 2.0*(nval-minnv)/(maxnv-minnv)
   val = [val,nval]
   if not(eof(1)) then readf,1,junk
;   if not(eof(1)) then readf,1,junk
;endwhile
 
endfor

end
;-----------------------------------------------------------


;-----------------------------------------------------------
pro read_data,nfil,img
;-----------------------------------------------------------

nlng = 12
nr = 10

;sp-------
filetype='MIXR'
spec = 'sp'

i=nfil
if (i lt 10) then tm = '000'+strtrim(string(i),2)
if ((i ge 10) and (i lt 100)) then tm = '00'+strtrim(string(i),2)
if (i ge 100) then tm = '0'+strtrim(string(i),2)
read_output,filetype,spec,tm,nlng,nr,lng,val,rad

pc = transpose([[lng],[rad]])
xy = cv_coord(FROM_POLAR=pc,/TO_RECT,/DEGREES)
;val = val(1:*)
x = reform(xy(0,*))
y = reform(xy(1,*))
triangulate,x,y,tr,b
arrsp = trigrid(x,y,val,tr,nx=401,ny=401)


;s3p-------
filetype='MIXR'
spec = 's3p'

read_output,filetype,spec,tm,nlng,nr,lng,val,rad

pc = transpose([[lng],[rad]])
xy = cv_coord(FROM_POLAR=pc,/TO_RECT,/DEGREES)
x = reform(xy(0,*))
y = reform(xy(1,*))
triangulate,x,y,tr,b
arrs3p = trigrid(x,y,val,tr,nx=401,ny=401)

arr = rotate(arrsp*arrs3p,7)


sz = size(arr)
xarr = min(x) + findgen(sz(1))*(max(x)-min(x))/(sz(1)-1)
yarr = min(y) + findgen(sz(1))*(max(y)-min(y))/(sz(1)-1)
xar = fltarr(sz(1),sz(2))
yar = fltarr(sz(1),sz(2))
for j = 0,sz(2)-1 do begin
   xar(*,j) = xarr
endfor
for i = 0,sz(1)-1 do begin
   yar(i,*) = yarr
endfor

wh = where(arr eq 0)
arr(wh) = max(arr)

wh = where((sqrt(xar^2 + yar^2) gt 6.0) and (sqrt(xar^2 + yar^2) lt 9.0))
minarr = min(arr(wh))

c = image(arr>minarr,xarr,yarr,/buffer,/current,rgb_table=33,axis_style=2,$
         xrange=[-11,11],yrange=[-11,11],/xsty,/ysty,$
         position=[0.1,0.15,0.8,0.85])

cb = colorbar(target=c,range=[min(arr(wh)),max(arr(wh))],orientation=1,$
             position=[0.82,0.15,0.85,0.85],textpos=1,font_size=14)
cb.title=filetype+' '+spec

print,max(arr),min(arr)

dr = 0.2
;for i = 0,29 do begin
;   f = circle(0,0,0.1+dr*i)
;   p = plot(f(0,*),f(1,*),'w7',/overplot,/buffer)
;endfor
;f = circle(0,0,5.9)
e = ellipse(0,0,major=5.9,/fill_background,/overplot,/buffer,/data)


;for i = 0,22 do begin
;   f = circle(0,0,10.0+dr*i)
;   p = plot(f(0,*),f(1,*),'w7',/overplot,/buffer)
;endfor

f1 = circle(0,0,9.9)
f2 = circle(0,0,20.0)

fx = [reform(f1(0,*)),reform(f2(0,*))]
fy = [reform(f1(1,*)),reform(f2(1,*))]

p = polygon(fx,fy,/fill_background,/data,linestyle=6)

t = text(5,0,'0',/data,font_size=14,alignment=0.5)
t = text(0,5,'270',/data,font_size=14,alignment=0.5)
t = text(-5,0,'180',/data,font_size=14,alignment=0.5)
t = text(0,-5,'90',/data,font_size=14,alignment=0.5)

img = c.CopyWindow()

end
;-----------------------------------------------------------


;main program
;@clr_win
xsz = round(1000/1.0)
ysz = round(1000/1.0)
nframe=41
nfrm0 = 0


XINTERANIMATE, SET=[xsz,ysz, nframe-nfrm0], /SHOWLOAD 

video_file = 'torus.mp4'
video = idlffvideowrite(video_file)
framerate = 7.5
;wdims = w.window.dimensions
stream = video.addvideostream(xsz, ysz, framerate)

cnt = 0
for i = 0,200,5 do begin
   w = window(window_title='torus',dimensions=[xsz,ysz],margin=0,$
              buffer=1)
   
   nfrm = i
   read_data,nfrm,img
   ;im = image(img)
   xinteranimate, frame = cnt, image = img
   print,'video:', video.put(stream, img)
   w.close
   cnt = cnt+1
endfor

video.cleanup
xinteranimate,/keep_pixmaps

end
