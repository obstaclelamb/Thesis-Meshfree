function [xI,BC_index] = sub_ReadNeutralInputFiles(NeutralFileName)
% Read-in the CAD/FEA neutral file in Patran neutral format
% e.g., NeutralFileName = 'FE_Neutral.dat';
%
filename = fullfile(NeutralFileName);  % open the full neutral file
T = readtable(filename); % read in table format in MATLAB
C = table2cell(T);       % convert table format to cell format
% Convert cell format to double precision to obtain coordinates
nLine = length(C); % number of lines in neutral file
C_double = cellfun(@str2num,C,'UniformOutput',false);
% Read discretization information
LineOfC = C_double{2};
NNODE = LineOfC(5); % number of nodes
NELEM = LineOfC(6); % number of elements

% 0429 test boundary
BC_index = [];

xI = zeros(NNODE,2); % nodal coordinates initialization
disp(['From ',filename,' file: #node is ',num2str(NNODE),', #element is ',num2str(NELEM)])
% Read the file 
iL = 2; idx_node = 1;
while (iL <= nLine) % read line by line from the input files
    LineMatrix = C_double{iL};
    if ~isempty(LineMatrix) % no empty line is read
    IDCARD = LineMatrix(1); 
    %% IDCARD = 01 is the coordinates list; Read coordinates
    if (IDCARD == 1) && length(LineMatrix) == 9
    iL = iL + 1;  % read next line 
    xI(idx_node,1:2) = C_double{iL}(1:2);
    idx_node = idx_node + 1; % next node
    iL = iL + 1;  % next line  
    elseif (IDCARD==21) && strcmp( C{iL+1},'BoundaryNode')
        number_boundary_node = C_double{iL}(3)/2;
        if mod(number_boundary_node,5) == 0
            number_of_line_for_BC = fix(number_boundary_node/5);  
        else
            number_of_line_for_BC = fix(number_boundary_node/5)+1;
        end
        for nL_idx = 1:number_of_line_for_BC
            BC_index = [BC_index;C_double{iL+1+nL_idx}(2:2:end)'];
        end
        iL = iL + (fix(number_boundary_node/5)+1);
    end
    end % end if ~isempty
    iL = iL + 1; % next line
end % end of reading each line
end % end of the function
