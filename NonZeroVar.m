function Non0Var=NonZeroVar(InputArray)

% if dim=2
%     InputArray=InputArray';
% end
Non0Num=sum(InputArray~=0);
VarGain=size(InputArray,1)./Non0Num;
InputArray=InputArray+(InputArray==0)*sparse(diag(mean(InputArray).*VarGain));
Non0Var=var(InputArray).*VarGain;