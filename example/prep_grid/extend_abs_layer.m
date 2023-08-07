function [bz1,bz2] = extend_abs_layer(bz1,bz2,dx,nx,nz,num_pml)

coef = 0.7;
for i=num_pml:-1:1
    dz2 = bz2(i+2,2)-bz2(i+1,2);
    slope2 = coef*dz2/dx;
    bz2(i,1) = bz2(i+1,1) - dx;
    bz2(i,2) = bz2(i+1,2) - slope2*dx;

    dz1 = bz1(i+2,2)-bz1(i+1,2);
    slope1 = coef*dz1/dx;
    bz1(i,1) = bz1(i+1,1) - dx;
    bz1(i,2) = bz1(i+1,2) - slope1*dx;
end

for i=nx-num_pml+1:nx
    dz2 = bz2(i-1,2)-bz2(i-2,2);
    slope2 = coef*dz2/dx;
    bz2(i,1) = bz2(i-1,1) + dx;
    bz2(i,2) = bz2(i-1,2) + slope2*dx;

    dz1 = bz1(i-1,2)-bz1(i-2,2);
    slope1 = coef*dz1/dx;
    bz1(i,1) = bz1(i-1,1) + dx;
    bz1(i,2) = bz1(i-1,2) + slope1*dx;
end
