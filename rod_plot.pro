
window, xsize=800, ysize=800

FOR t=0,18 DO BEGIN
  
  idx=string(FORMAT='(I03)', t)
  data = READ_ASCII('D:\Desktop\jamming-master\1\'+idx+'.txt')
  data = data.FIELD1

  x_1 = data[2:3,*]
  x_2 = data[4:5,*]
  
  s = SIZE(data,/DIMENSIONS)

  plot, [x_1[0,0],x_2[0,0]], [x_1[1,0],x_2[1,0]],thick=2, xrange=[0,120],yrange=[0,120],/NODATA,/ISOTROPIC

  FOR i=0,s[1]-1 DO BEGIN
    oplot, [x_1[0,i],x_2[0,i]], [x_1[1,i],x_2[1,i]],color=255,thick=4
  ENDFOR
wait,0.1
ENDFOR



END