function [output,est_p,edge1,edge2,Nr]  = phaseunwrap2d(input)
% usage: [output,est_p,edge1,edge2,Nr]  = phaseunwrap2d(input)
% input: a 2D complex image
% output: the unwrapped phase values
% note: this is good only for epi Nyquist artifact removal;
% author: nan-kuei chen

p = input;
L2 = find(p==0);
p(L2)=1;
p = p./abs(p);
xdim = size(p,1);
ydim = size(p,2);

p2 = zeros(xdim*2,ydim*2);
p2(1:xdim,1:ydim)=p;
p2(1:xdim,ydim+1:ydim*2)=flipdim(p,2);
p2(xdim+1:xdim*2,1:ydim)=flipdim(p2(1:xdim,1:ydim),1);
p2(xdim+1:xdim*2,ydim+1:ydim*2)=flipdim(p2(1:xdim,ydim+1:ydim*2),1);

gmatrix = zeros(xdim*2,ydim*2);
for cntx = 1:xdim*2
 for cnty = 1:ydim*2
  cx = cntx-(xdim+1);
  cy = cnty-(ydim+1);
  gmatrix(cntx,cnty)= (cx^2) + (cy^2);
 end
end
gmatrix(xdim+1,ydim+1)=0.01;

p2real = real(p2);
p2imag = imag(p2);
t1 =p2real.*qft(qift(p2imag).*gmatrix);
t2 =p2imag.*qft(qift(p2real).*gmatrix);
t3 = t1-t2;
t4 =qft(qift(angle(p2)).*gmatrix);
t5 = t3-t4;
t6 =qft(qift(t5)./gmatrix);
t7=2*pi*round(t6/2/pi); 
t8 =qft(qift(t3)./gmatrix);
output = angle(p)+t7(1:xdim,1:ydim);
edge1 = t3(1:xdim,1:ydim);
edge2 = t5(1:xdim,1:ydim);
Nr = t6(1:xdim,1:ydim);
est_p = t8(1:xdim,1:ydim);

return
