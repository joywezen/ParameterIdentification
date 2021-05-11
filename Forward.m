function [ sensor_x,sensor_y ] = Forward( varargin )

% This function is meant to be used for computing the predicted
% measurements.

%   This function takes 2 inputs and produces output of displacement field
%   Inputs are:
%          varargin{1} - x-position of the center of circular disc
%          varargin{2} - radius of the circular disc

% An example of a call to this function is:
%    [y] = Forward( 6, 0.3);
%       where, 'y' is a matrix of 2 columns in which the first column
%               corresponds to the x-displacements of nodes and second
%               column corresponds to the y-displacement of nodes.

    % position value
    p_true = varargin{1}; 
   % radius value
    r_true=varargin{2};
    % sample number
    i=varargin{3};
    
    dir=varargin{4};
    edge=varargin{5};
    xforce=varargin{6};
    yforce=varargin{7};
  % read input file in a string 
   str = fileread('ISTRIP');   
  % find line with position in the previous string and modify it with the sample value
   str = strrep(str, 'p=p', sprintf('p=%d',p_true)); 
  % find line with radius in the previous string and modify it with the sample value
   str = strrep(str, 'r=r', sprintf('r=%d',r_true));
 % dummy??!!
   str = strrep(str, '1 l x y', sprintf('%d %d %d %d',dir,edge,xforce,yforce)); 
 % open new file
  fid = fopen(['new_' sprintf('filename%d',i)], 'w');
  % write the modified string in new file
  fwrite(fid, str, '*char');         
  % close file 
  fclose(fid);
  % get matlab library path getenv: get environment variable  
   %MatlabPath = getenv('LD_LIBRARY_PATH');
  % set matlab library path
   %setenv('LD_LIBRARY_PATH',getenv('PATH'));
 % call  feap without interface
   brstr=strcat('feap', sprintf(' -inew_filename%d',i),sprintf(' -ooutputs%d',i),sprintf(' -rrout%d -ssout%d -ppout%d',i,i,i));
   %brstr=strcat('./feap',sprintf(' n -inew_filename%d -ooutputs%d -rrout%d -ssout%d -ppout%d',i,i,i,i,i),sprintf(' y'));
   system(brstr);
 % set back matlab path
   %setenv('LD_LIBRARY_PATH',MatlabPath);
  % delete copied input file (not to take space on main memory)
  delete(sprintf('new_filename%d',i));
  delete(sprintf('new_filename%d.rev%d',i));
  delete(sprintf('Lutputs%d',i));
  delete('feap.neu');
  delete('feapname');
  %delete(sprintf('Mew_filename%d',i));
  % open output file 
  fid=fopen(sprintf('outputs%d',i));
  % scan content 
  C = textscan(fid, '%s','delimiter', '\n');
  %###############the former guy forgot to close the file##################
  fclose(fid);
  delete(sprintf('outputs%d',i));
  %############################idiot!!#####################################
  %I HATE HIM: Muhammad Umar Anwar 
  % find the line where displacement begin 
  IndexCS = strfind(C{1}, 'Node     1 Coord     2 Coord     1 Displ     2 Displ');
  IndexS = find(not(cellfun('isempty', IndexCS)));
  IndexCE = strfind(C{1}, '*Command   6 * stre NODE');
  IndexE = find(not(cellfun('isempty', IndexCE)));
 % read values
   for j=1:(IndexE-IndexS-1)
    Rel_exp{j,1}=C{1}{IndexS+j};
   end
 % transform to numbers
  z_exp=char(Rel_exp);
  y_b=str2num(z_exp);
 % global position of sensors (does not depend on mesh)
  xp=[3 3.5 3.6 3.7 4 4.5 5 5.3 5.5 6 7 7.5 8 8.2 8.5 8.8 8.96 9.2 10 12 12.7 13.5 13.75 14 14.5 14.7 14.85 15 15.2 15.4 15.5 15.6 16 16.2 16.5 16.6 16.75 16.8 17 18];
  yp=[5 5 2.2 3.5 7.5 3 8.8 3.8 2.5 8 9 1.5 2.4 6.25 2.15 1.8 7.9 8.15 7.5 2.5  3  8.3  2.1 5 5 6.5 3.5 5 2.2 5 5 6.5 2.9 5 7.2 1.5 8 5 1.75 1.9];
 % interpolate Feap displacements to the sensor positions
  xx=y_b(:,2);
  yy=y_b(:,3);
  disp_x=y_b(:,4);
  disp_y=y_b(:,5);
   % for j=1:length(xp) % linear interpolation
  sensor_x = griddata(xx,yy,disp_x,xp,yp)'; %%#ok<GRIDD> % predicted measurement
  sensor_y = griddata(xx,yy,disp_y,xp,yp)';
    %end
    %zi = transpose(zi); %#ok<NASGU>
end

