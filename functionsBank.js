
Math.abs(x)+.1*Math.abs(y)
//bruit
Math.random()*amplitude

x*x+y*y//sphere
Math.floor(12*x*x)+Math.floor(12*y*y)//stepshpere
.1*x*x+y*y//cigar
Math.pow(Math.abs(x),2)+Math.pow(Math.abs(y),12)//diffpow
-.01*x+Math.abs(y)//sharpR

//4tests
x*x+10*y*y//ellipsoid, x²+alpha*y²
Math.pow(100*x*x-y*10,2)+Math.pow(1-x*10,2)/100//rosenbrock, alpha is the last factor
4*x*x+Math.pow(2*y,2+10)//diffpow, y^(2+alpha)
.1*(20+100*x*x-10*Math.cos(2*Math.PI*10*x)+100*y*y-10*Math.cos(2*Math.PI*10*y))//Rastrigin
.1*(20+100*x*x-10*Math.cos(2*Math.PI*10*x)+100*y*y-10*Math.cos(2*Math.PI*10*y))+Math.random()*.2//Rastrigin bruitée
