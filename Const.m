function ConstValue = Const(ConstName)

switch ConstName
    case 'Null'
        ConstValue=0;
    case 'Recon'
    ConstValue = 100;
    case 'Thresh'
        ConstValue = 400;
    case 'Demotion'
        ConstValue = 300;
    case 'Filter'
        ConstValue = 200;
    case 'Display'
        ConstValue = 500;
    case 'LIMPRatio'
        ConstValue=1;
    case 'LIMPThreshP'
        ConstValue=0.3;
    case 'DFRMRes'
        ConstValue=2/128;
    case 'DFRMGPU'
        ConstValue=1;
    case 'DFRMThresh'
        ConstValue=0.1;
    case 'DFRMGain'
        ConstValue=0.1;
    case 'DEFRItNum'
        ConstValue=2;
    case 'THQDOCT'
        ConstValue=2;
    otherwise
        error('Unknown Constant Name');
end
