function Non0Var=NonZeroVar(InputArray)

% if dim=2
%     InputArray=InputArray';
% end
Non0Num=sum(InputArray~=0);
InputArray(InputArray==0)=mean(InputArray(InputArray~=0));
Non0Var=var(InputArray).*size(InputArray,1)/Non0Num;

