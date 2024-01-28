using PyPlot

# ГОСТ 1497-84, образец плоский, тип I, образец номер №6
mm = 1e-3;
B = 40mm;
l0 = 140mm;
a0 = 20mm;
b0 = 30mm;
l = l0+2*sqrt(a0*b0);
h1 = 80mm;
r = (B-b0)/2;

points = [
    0             0;  
    h1            0; 
    # h1+r          0; 
    h1+r          r; 
    h1+r+l        r; 
    # h1+r+l        0; 
    h1+r+l+r      0; 
    h1+r+l+r+h1   0; 
    h1+r+l+r+h1   B; 
    h1+r+l+r      B; 
    # h1+r+l        B; 
    h1+r+l        B-r;
    h1+r          B-r;
    h1+r          B; 
    h1            B;  
    0             B;  
    0             0;
]

plot(points[:, 1], points[:, 2])
gca().set_aspect(1)