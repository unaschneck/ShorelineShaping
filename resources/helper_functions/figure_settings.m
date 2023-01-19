function figure_settings(fontsize_num,fontstyle,grid_yes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective:
%   Makes figures pretty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   fontsize = font size for text in figures 
%   fontstyle = style (e.g. bold)
%   grid_yes = Y (grid on) or N (grid off)
% outputs:
%   none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fh = findall(0,'Type','Figure');
set( findall(fh, '-property', 'fontsize'), 'fontsize', fontsize_num)
fh = findall(0,'Type','Figure');
set( findall(fh, '-property', 'FontWeight'), 'FontWeight',fontstyle)
box on;
if nargin == 3 && grid_yes == 1
    grid on;
elseif nargin <3
    grid off;
end

end