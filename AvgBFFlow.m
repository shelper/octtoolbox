%Combine Back forward scan results
function BFAvgFlow=AvgBFFlow(FlowB, FlowF)
if get(handles.BackForward,'Value')
        DivArray=ones(size(FlowB))+(FlowB~=0).*(FlowF~=0);
        BFAvgFlow=(FlowB-FlowF)./DivArray;
        BFAvgFlow=BFAvgFlow+abs(FlowB+FlowF)>1; %correct the flow area where B&F flow is relatively high and has same sign
        BFAvgFlow=BFAvgFlow-2*(BFAvgFlow>1);
end    

