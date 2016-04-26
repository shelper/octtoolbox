close all;
dia = 3;
aline_num = 360;
cxs = dia * sin(linspace(0, 2*pi, aline_num));
cys = dia * cos(linspace(0, 2*pi, aline_num));
figure;plot(cxs, cys);
mov = [0.1, 0.1];
cxe = dia * sin(linspace(0, 2*pi, aline_num))+mov(1);
cye = dia * cos(linspace(0, 2*pi, aline_num))+mov(2);
hold on;plot(cxe, cye,'r');
 
x = [cxs; cxe];
y = [cys; cye];
hold on;plot(x, y);

cxe_shift = circshift(cxe, [0, 3]);
cye_shift = circshift(cye, [0, 3]);
x = [cxs; cxe_shift];
y = [cys; cye_shift];
figure;plot(x, y);

error = 0;
delta_s = sqrt((cxs-cxe).^2 + (cys-cye).^2);
figure;plot(delta_s);
for i = -10:10
    cxe_shift = circshift(cxe, [0, i]);
    cye_shift = circshift(cye, [0, i]);
    delta_c = sqrt((cxs-cxe_shift).^2 + (cys-cye_shift).^2);
    error = abs(delta_c - delta_s);
    hold on;plot(delta_c);
end


