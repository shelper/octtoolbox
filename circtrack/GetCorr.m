function p = GetCorr(line1, line2, DepthSpan)
BandNum = floor(size(line1, 1)/DepthSpan);
Corr = zeros(BandNum, 1);
Weight = zeros(BandNum,1); 
for jDepth = 1:BandNum 
    SubLine1 = line1((jDepth-1)*DepthSpan+1:jDepth*DepthSpan);
    SubLine2 = line2((jDepth-1)*DepthSpan+1:jDepth*DepthSpan);
    p1 = sum((SubLine1-mean(SubLine1)) .*(SubLine2-mean(SubLine2)));
    p2 = sqrt(sum((SubLine1-mean(SubLine1)).^2) * sum((SubLine2-mean(SubLine2)).^2));
    Corr(jDepth) = p1/p2;
    Weight(jDepth) = mean(SubLine1+SubLine2);     
end
p = sum(Corr.* Weight)/sum(Weight);

