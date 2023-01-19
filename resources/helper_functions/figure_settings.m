function figure_settings(fontsize_num,fontstyle,grid_yes)



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