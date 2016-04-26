function Edge=FindEdge(InputImage,Method)
%default use Threshold Method
if nargin < 2
    Method='Threshold';
end
%preprocessing of Image 
FilterSize=5;
InputImage=imfilter(InputImage,fspecial('disk',FilterSize/2)); % TBF?-Find optimized FilterSize here to close the holes and smooth the edge
%Methods for initialzie the edges.
switch Method
    case 'Canny'
        BW=edge(InputImage,'canny');%canny edge detection
        PossibleEdgeNum=round(mean(sum(BW))/10);
    case 'DiffCanny'
        DiffImage=InputImage-circshift(InputImage,FilterSize);
        BW=edge(DiffImage,'canny');%canny edge detection with Diff
        PossibleEdgeNum=round(mean(sum(BW))/10);
    case 'Threshold'
        BW=imclose(im2bw(InputImage,graythresh(InputImage)),strel('disk',FilterSize));%edge withthresholding  
        BW=diff(BW);
        PossibleEdgeNum=round(mean(sum(BW>0))*2);
    otherwise
        disp('Unknown edge searching method.');
end
[foo,I]=sort(BW,'descend');
I=I(1:PossibleEdgeNum,:);
PossibleEdges=zeros(size(I));
PossibleEdges(I(:,1)>10,1)=I(I(:,1)>10,1);

%find edge by minimize position change
for j=1:size(I,2)-1
    %check the possible neighbour points with minimum position shift.
    SpotNum1=min(PossibleEdgeNum,sum(PossibleEdges(:,j)>10));%remove the fake edges detected by sort.
    SpotNum2=min(PossibleEdgeNum,sum(I(:,j+1)>10));
    if ~SpotNum1&&SpotNum2
        PossibleEdges(:,j+1)=I(:,j+1);
    else 
        %StepLineNum=min(PossibleEdgeNum,SpotNum1*SpotNum2);
        PosShift=abs(repmat(PossibleEdges(1:SpotNum1,j)',[SpotNum2,1])-repmat(I(1:SpotNum2,j+1),[1,SpotNum1]));
        %find the points with minumum position shift.
        [foo,i]=sort(PosShift(:));
%         i=i(1:min(sum(foo<10),PossibleEdgeNum));%prevent sudden change bigger than 10 this cause error with a bright line.
        i=i(1:min(end,PossibleEdgeNum));
        if ~isempty(i)
            i1=ceil(i/SpotNum2);
            i2=i+SpotNum2-i1*SpotNum2;
            PossibleEdges(:,j+1)=I(i2(1),j+1);
            PossibleEdges(i1,j+1)=I(i2,j+1);
        end
    end
end

%delete the zero sections
for j=1:PossibleEdgeNum
    [foo,col1]=find(diff(PossibleEdges(j,:)<10,1,2)==1);
    col1=col1+1;
    if PossibleEdges(j,end)==0
        PossibleEdges(j,col1(end):end)=PossibleEdges(j,col1(end)-1);
        col1(end)=[];
    end
    [foo,col2]=find(diff(PossibleEdges(j,:)<10,1,2)==-1);
    if PossibleEdges(j,1)==0
        PossibleEdges(j,1:col2(1))=PossibleEdges(j,col2(1)+1);
        col2(1)=[];
    end
    for i=1:size(col1,2)
        PossibleEdges(j,col1(i):col2(i))=interp1([col1(i)-1,col2(i)+1],[PossibleEdges(j,col1(i)-1),PossibleEdges(j,col2(i)+1)],col1(i):col2(i));
    end
end

%find the smoothest edge by combing smoothest sections from each edge
for j=1:size(PossibleEdges,2)-1
    iDup=find((diff(PossibleEdges(:,j))~=0&(diff(PossibleEdges(:,j+1)))==0)==1);
    if ~isempty(iDup)
        for i=1:size(iDup,1)
            if sum(diff(PossibleEdges(iDup(i), 1:j+1)).^2)>sum(diff(PossibleEdges(iDup(i)+1, 1:j+1)).^2)  
                PossibleEdges(iDup(i), 1:j+1)=PossibleEdges(iDup(i)+1, 1:j+1);
            else
                PossibleEdges(iDup(i)+1, 1:j+1)=PossibleEdges(iDup(i), 1:j+1);
            end
        end
    end
end

%choose the smoothest edge as result
[foo,i]=min(sum(diff(PossibleEdges,1,2).^2,2));
Edge=round(filtfilt(ones(1,20)/20,1,PossibleEdges(i,:)));


%% Debug scripts
% Margin=1;
% foo=I((I(:,1)>10),Margin);
% while isempty(foo)
%     PossibleEdges(:,1)=0;
%     Margin=Margin+1;
%     foo=I((I(:,Margin)>10),Margin);
% end
% PossibleEdges(1:size(foo,1),1:Margin)=repmat(foo,[1,Margin]);
% PossibleEdges(size(foo,1)+1:end,1:Margin)=repmat(PossibleEdges(size(foo,1),1:Margin),[PossibleEdgeNum-size(foo,1),1]);


%find the smoother edge when crossed
% for i=1:2000
%     InputImage(Edge(i),i)=65535;
% end
% figure;imshow(InputImage);
%     iDup=find((diff(PossibleEdges(:,j))~=0&(diff(PossibleEdges(:,j+1)))==0)==1);
%     if ~isempty(iDup)
%         for i=1:size(iDup,1)
%             if sum(diff(PossibleEdges(iDup, FilterSize:j+1)).^2)>sum(diff(PossibleEdges(iDup+1, FilterSize:j+1)).^2)   
%                 PossibleEdges(iDup, FilterSize:j+1)=PossibleEdges(iDup+1, FilterSize:j+1);
%             else
%                 PossibleEdges(iDup+1, FilterSize:j+1)=PossibleEdges(iDup, FilterSize:j+1);
%             end
%         end
%     end

%             PossibleEdges(:,j+1)=I(:,j+1);
%             [foo,i]=sort(i1);            
%             PossibleEdges(1:size(i2,1),j+1)=I(i2(i),j+1);  
%             PossibleEdges(size(i2,1):end,j+1)=PossibleEdges(size(i2,1),j+1);

%             PossibleEdges(1:size(i1,1),1:j)=PossibleEdges(i1,1:j);
%             PossibleEdges(size(i1,1):end,1:j)=repmat(PossibleEdges(size(i1,1),1:j),[PossibleEdgeNum-size(i1,1)+1,1]);
            
%             PossibleEdges(1:size(i2,1),j+1)=I(i2,j+1);  
%             PossibleEdges(size(i2,1):end,j+1)=PossibleEdges(size(i2,1),j+1);

%         else
%             foo=I((I(:,j+1)>10),j+1);
%             PossibleEdges(1:size(foo,1),j+1)=foo;
%             PossibleEdges(size(foo,1):end,j+1)=PossibleEdges(size(foo,1),j+1);
%             PossibleEdges(:,j+1)=I(:,j+1);


% Edge=PossibleEdges(1,:);
% NegShift=0;
% for j = 2: PossibleEdgeNum
%     Cross=diff(Edge~=PossibleEdges(j,:));
%     [foo, iCross]=sort(Cross,'descend');
%     iCrossPos=iCross(1:sum(foo>0));
%     iCrossNeg=iCross(end-sum(foo<0)+1:end)+1;
%     if  Edge(1)~=PossibleEdges(j,1)    
%         if sum(diff(Edge(1:iCrossNeg(1))).^2)>sum(diff(PossibleEdges(j,1:iCrossNeg(1))).^2)   
%             Edge(1:iCrossNeg(1))=PossibleEdges(j,1:iCrossNeg(1));
%         end
%         NegShift=1;
%     end
%     for i=1:size(iCrossNeg,2)-NegShift
%         if sum(diff(Edge(iCrossPos(i):iCrossNeg(i+NegShift))).^2)>sum(diff(PossibleEdges(j,iCrossPos(i):iCrossNeg(i+NegShift))).^2)           
%             Edge(iCrossPos(i):iCrossNeg(i+NegShift))=PossibleEdges(j,iCrossPos(i):iCrossNeg(i+NegShift));
%         end
%     end
%     if  Edge(end)~=PossibleEdges(j,end)    
%         if sum(diff(Edge(iCrossPos(end):end)).^2)>sum(diff(PossibleEdges(j,iCrossPos(end):end)).^2)           
%             Edge(iCrossPos(end):end)=PossibleEdges(j,iCrossPos(end):end);
%         end
%     end
% end
% 
% tmpImage=InputImage;
% for i=1:2000 
%     tmpImage(Edge(1,i)+1,i)=65535;
% end
% figure;imshow(tmpImage);
%     
% 
% [foo,i]=min(sum(abs(diff(PossibleEdges,1,2)),2));
% figure;plot(PossibleEdges(i,:));
%     
% tmpImage=InputImage;
% for i=1:2000
% tmpImage(PossibleEdges(1,i)+1,i)=65535;
% end
% figure;imshow(tmpImage)
% 
%     
% % Image1=imfilter(InputImage, fspecial('disk',5));
% % % BW = EdgeImage(Image1,'canny',0.3);figure;imshow(BW);
% % BW = EdgeImage(Image1,'canny',0.4)
% 
% 
% 
% 
% % Image2=Image1-circshift(Image1,10);
% % Image3=diff(im2bw(Image2,0.4));
% % [Image3,i]=sort(Image2,1,'descend');
% % 
% % i=i(1:5,:);
% % 
% % for j=1:2000
% %     m1=repmat(i(:,j),[1,5]);
% %     m2=repmat(i(:,j+1)',[5,1]);
% %     m3=abs(m1-m2);
% %     
% % 
% % Image3=imfilter(Image2,ones(20,1)/20);
% % 
% % for i = 1:sizeof(Image3,2)
% %     ALine=Image3(:,i);
% %     
% %     
% %     [max,I1]=max(Image1);[max,I2]=max(Image2);[max,I1]=max(Image1);[max,I1]=max(Image1);[max,I1]=max(Image1);
% %     Image1(I1,)=0;
% %     [max,I2]=max(ALine);ALine(I2)=0;
% %     [max,I3]=max(ALine);ALine(I3)=0;
% %     [max,I4]=max(ALine);ALine(I4)=0;
% %     [max,I5]=max(ALine);ALine(I5)=0;
% %     
% %     
% %     PossibleEdges(i)=
% % end
% % 
% % 
% % 
% % figure; imagesc(Image2);
% % figure;plot(Image3(:,100));
% % hold on; plot(Image1(:,100));