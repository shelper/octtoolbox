for i=0:31
    if i<10
        movefile(['OCTImage000',num2str(i),'.dat'], ['.\OCTImage00',num2str(i),'.dat']);
    else
        movefile(['OCTImage00',num2str(i),'.dat'], ['.\OCTImage0',num2str(i),'.dat']);
    end    
end


fid = fopen('T:\8 Doppler OCT\Doppler RawData from THQ\RawData20110802161453\FlowImage000THQ.dat', 'rb');
Image = fread(fid, [885 512], 'single');
fclose(fid);

DisplayRange=GetDisplayRange(Image, [0.1 0]);
figure;imshow(Image,'DisplayRange',DisplayRange);

% figure;imshow(Image,[40 80]);
% figure;plot(mean(Image(:,1550:1560),2));