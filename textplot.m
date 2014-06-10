classdef textplot < handle
    
    %% properties
    
    properties (Access = protected)
        
        title  = ''
        xlabel = ''
        ylabel = ''
        
        grid = true
        box  = true
        
        marker    = ''
        linestyle = '-'
        
        unicode = true
        
        xscale = 'linear'
        yscale = 'linear'
        
        xmin
        xmax
        ymin
        ymax
        
    end
    
    properties (Hidden, Access = protected)
        
        xdata
        ydata
        
        xspacing = 2
        yspacing = 2
        
        xsz
        ysz 
                
        valid_markers = {'.' 'o' 'x' '+' '*' 's' 'd' 'v' '^' '<' '>' 'p' 'h'}
        valid_linestyles = {'-.' ':' '--' '-'}
        
        plotstring
                
    end
    
    %% methods
    
    % public methods
    methods
        
        %% basic operation
        
        % constructor
        function obj = textplot(varargin)
             
            n = 0;
            
            % Call of the form 
            %    textplot(xdata, 'property', value, ...)
            %    textplot(xdata, '-.s', 'property', value, ...)            
            if nargin == 1 || ischar(varargin{2})
                obj.ydata = varargin{1};
                if ~isnumeric(obj.ydata)
                    error('textplot:nonnumeric_ydata',...
                        'Non-numeric ydata received.'); 
                end
                obj.xdata = 1:numel(obj.ydata);
                                
                if nargin > 1                   
                    n = obj.parse_linestyle(varargin{2}); end
             
            % Call of the form 
            %    textplot(xdata, ydata, 'property', value, ...)
            %    textplot(xdata, ydata, '-.s', 'property', value, ...)                
            elseif (nargin >= 2 && isnumeric(varargin{2})) || ischar(varargin{3})
                obj.xdata = varargin{1};
                obj.ydata = varargin{2};
                if ~isnumeric(obj.xdata) || ~isnumeric(obj.ydata)
                    error('textplot:nonnumeric_xydata',...
                        'Non-numeric x- or y-data received.'); 
                end
                if numel(obj.xdata) ~= numel(obj.ydata)
                    error('textplot:dimension_mismatch',...
                        'X and Y data must have the same number of elements.'); 
                end
                
                n = 1;
                if nargin > 2                    
                    n = n+obj.parse_linestyle(varargin{3}); end    
                
            end
            
            % default minima/maxima
            obj.xmin = min(obj.xdata(:));   obj.ymin = min(obj.ydata(:));
            obj.xmax = max(obj.xdata(:));   obj.ymax = max(obj.ydata(:));           
                        
            % parse any additional property/value pairs 
            % NOTE: indexing with 'end' makes the array empty 
            % for nargin < 3
            pvpairs = varargin(2+n:end);            
            if ~isempty(pvpairs)
                obj.set(pvpairs{:}); end
            
        end
        
        % display
        function disp(obj)
            obj.plotstring = obj.draw;
            disp(obj.plotstring)
        end
                
               
        %% get/set
        
        function set(obj, varargin)
            
            % get all valid properties
            class_props = ?textplot;
            class_props = class_props.Properties;
            valid_props = class_props(cellfun(@(x)~x.Hidden, class_props));
            valid_props = cellfun(@(x) x.Name, valid_props, 'UniformOutput', false);
                        
            % display ALL 
            if (nargin == 1)
                % TODO                
                return
            end
            
            % display valid and default values for given option
            if (nargin == 2)
                % TODO                
                return
            end
            
            % these might become inconsistent after parsing; save 
            % current values
            xmin_P = obj.xmin;  ymin_P = obj.ymin;
            xmax_P = obj.xmax;  ymax_P = obj.ymax;
            
            % parse property/value pairs            
            props = varargin(1:2:end);
            vals  = varargin(2:2:end);
            
            if numel(props) ~= numel(vals)
                error('textplot:invalid_property_value_pairs',...
                    'Textplot expects options in pairs of (''property'', value, ...)');
            end
            
            for ii = 1:numel(props)
                
                property = props{ii};
                value    = vals {ii};
                
                switch lower(property)
                    
                    case {'xlabel' 'ylabel' 'title'}
                        if ~ischar(value) && ~iscellstr(value)
                            error('textplot:invalid_argument',...
                                ['Property ''' property ''' must have value of type char.']);
                        end
                        obj.(property) = value;
                        
                        
                    case {'grid' 'unicode' 'box'}
                        if ischar(value)
                            if any(strcmpi(value, {'on' 'off'}))
                                value = strcmpi(value, 'on');
                            else
                                error('textplot:invalid_argument',...
                                    ['Property ''' property ''' must equal ''on'' or ''off''.']);
                            end
                            
                        elseif isscalar(value)
                            if isnumeric(value)
                                value = value~=0;
                                
                            elseif ~islogical(value)
                                error('textplot:invalid_argument',...
                                    ['Property ''' property ''' must equal true or false.']);
                            end
                            
                        else
                            error('textplot:invalid_argument',...
                                ['Property ''' property ''' must equal true or false, ''on'' or ''off'', or a scalar numeric value.']);
                            
                        end
                        obj.(property) = value;
                        
                        
                    case 'marker'
                        if ~ischar(value) || ~isscalar(value)
                            error('textplot:invalid_argument',...
                                'Property ''marker'' must be specified with a single character.');
                        end
                        obj.marker = value;
                        
                    case 'linestyle'
                        if ~ischar(value)
                            error('textplot:invalid_argument',...
                                'Property ''linestyle'' must have a value of type char.');
                        end
                        obj.linestyle = value;
                    
                        
                    case {'xscale' 'yscale'}
                        if ~ischar(value) || ~any(strcmpi(value, {'linear', 'log'}))
                            error('textplot:invalid_argument',...
                                ['Property ''' property ''' must equal ''linear'' or ''log''.']);
                        end
                        obj.(property) = lower(value);
                        
                        
                    case {'xmin' 'xmax' 'ymin' 'ymax'}
                        if ~isnumeric(value) || ~isscalar(value)
                            error('textplot:invalid_argument',...
                                ['Property ''' property ''' must be a numeric scalar.']);
                        end
                        if ~isfinite(value)
                            error('textplot:invalid_axis_limits',...
                                'Axes limits must be non-NaN, strictly increasing values.');
                        end
                        obj.(property) = value;
                        
                        
                        
                    otherwise
                        warning('textplot:invalid_property',...
                            'Unknown property: ''%s''. Ignoring...', property);
                        continue;
                        
                end
                
            end
            
            
            % check for inconsistent options
            if obj.xmin > obj.xmax || obj.ymin > obj.ymax
                
                obj.xmin = xmin_P;  obj.ymin = ymin_P;
                obj.xmax = xmax_P;  obj.ymax = ymax_P;
                
                error('textplot:invalid_axis_limits',...
                    ['Axes limits must be non-NaN, strictly increasing values. ',...
                    'Axes limits have been reset.']);
            end
            
        end
        
        function get(obj, varargin)
            props = properties(obj);
            
            % display ALL current values
            if nargin == 1
                % TODO                
                return
                
            % incorrect usage
            elseif nargin > 3
                error('',...
                    '');
            end
            
            % check & return requested property
            % TODO            
            
        end
              
        
    end
    
    % private methods
    methods (Access = private)
        
        %% initialization
        
        function bool = parse_linestyle(obj, str)
            
            % str could indeed be compact linestyle string
            if ischar(str) && ~isempty(str)
                     
%                 % find linespec
%                 strcmp(obj.valid_linestyles, '.')
%                 
%                 inds = regexp(str, obj.valid_linestyles);
%                 
%                 
%                 % find marker
%                 obj.marker = regexp(str, obj.valid_markers);
     
bool = false;
                
            
            % str is something else
            else 
                bool = false;
                return
            end            
        end
        
        %% draw methods
        
        function str = draw(obj)
            
            % update command window sizes
            sz  = get(0, 'CommandWindowSize');
            obj.xsz = sz(1) - 2*obj.xspacing;            
            obj.ysz = min(sz) - 2*obj.yspacing;
            
            % generate all components            
            ttl = obj.generate_title;                        
            [xax,yax, xtop,yright, xtick,ytick] = obj.generate_axes(...
                obj.xsz, ...
                obj.ysz-size(ttl,1));
            plt = obj.generate_plot(...
                size(xax,2), ...
                size(yax,1) - size(xax,1),...
                xtick,ytick);
                        
            hspace = repmat(' ', size(yax,1), obj.xspacing);
            vspace = repmat(' ', obj.yspacing, 1);
           
            % and combine this into a single string
            if obj.box
                plt(1,:) = xtop;
                str = char(vspace, ttl, [hspace, yax, char(plt, xax), yright], vspace);
                
            else
                str = char(vspace, ttl, [hspace, yax, char(plt, xax)], vspace);
                
            end
            
        end
        
        % generate title area
        function ttl = generate_title(obj)
            if ~isempty(obj.title)                
                center = @(str) [repmat(' ', 1, floor((obj.xsz-size(str,2))/2)) str];
                if iscellstr(obj.title)
                    ttl = char( cellfun(center, obj.title, 'Uniformoutput', false) );
                else
                    ttl = center(obj.title);
                end
                
            else
                ttl = '';
            end
        end
        
        % create plot area
        function plt = generate_plot(obj, width,height, xtick,ytick)
            
            % initial plot is empty (possibly with gridlines)
            plt = repmat(' ', height,width);
            if obj.grid
                
                if obj.unicode
                    xgrid = char(8231);
                    ygrid = char(8231);
                else
                    xgrid = ':';
                    ygrid = '.';                    
                end               
                
                plt(:, xtick) = xgrid;
                plt(ytick,1:2:end) = ygrid;
                                 
                % possibly introduces some doubles; remove them
                dbl = plt(ytick,xtick)==xgrid & plt(ytick,xtick-1)==ygrid;                
                for ii = 1:numel(ytick)
                    plt(ytick,xtick(dbl(ii,:))) = ' '; end                
            end
            
            % increase in x/y value per single character
            xincr = (obj.xmax-obj.xmin)/width;
            yincr = (obj.ymax-obj.ymin)/height;
            
            % bin the data into character-wide cells
            xN = histc(obj.xdata, 0:width);
            yN = histc(obj.ydata, 0:height);
            
            
            return
            
            
            if obj.unicode
            end
            obj.linestyle
            obj.marker            
        end
        
        % create axes
        function [xax, yax, xtop, yright, xtick_locations, ytick_locations] = ...
                generate_axes(obj, width, height)
            
            % generate labels
            
            xlen =  width;
            ylen = height;  
            
            for ii = 1:2 % yes, iteratively...think about it. 
            
                % x-label  
                xlbl = [];
                if ~isempty(obj.xlabel) 

                    center = @(str) [repmat(' ', 1,floor((xlen-size(str,2))/2)) str];                

                    if iscellstr(obj.xlabel)
                        xlbl = char( cellfun(center, obj.xlabel, 'Uniformoutput', false) );
                    else
                        xlbl = center(obj.xlabel);
                    end  

                    xlbl = char(' ', xlbl);
                                        
                end

                % y-label
                ylbl = [];
                if ~isempty(obj.ylabel)

                    center = @(str) [...
                        repmat(' ', 1,floor((ylen-size(str,2))/2))...
                        str ...
                        repmat(' ', 1,ceil((ylen-size(str,2))/2))];  

                    if iscellstr(obj.ylabel)
                        ylbl = char( cellfun(center, obj.ylabel, 'Uniformoutput', false) );
                    else
                        ylbl = center(obj.ylabel);
                    end

                    ylbl = char(ylbl, ' ').';                              
                end
                
                % y-tick
                if strcmp(obj.yscale, 'linear')
                    
                    
                else
                    % TODO
                end
                  
                % adjust lenghts with label sizes, and tick lengths
                % TODO: tick lengths are variable.....
                xlen =  width - size(ylbl,2) - 11;
                ylen = height - size(xlbl,1) - 2;

            end
            
            
            % add axes lines/ticks
                            
            if obj.unicode
                % regular axes
                xnorm   = char(9472); % ─
                xtick   = char(9524); % ┴
                                
                ynorm   = char(9474); % │
                ybottom = char(9492); % └                
                ytick   = char(9500); % ├
                
                % box axes
                xtoptick     = char(9516); % ┬
                
                ytop         = char(9484); % ┌
                yrighttop    = char(9488); % ┐
                yrightbottom = char(9496); % ┘
                yrighttick   = char(9508); % ┤                

            else
                % regular axes
                xnorm = '-';
                xtick = '|';

                ynorm   = '|';
                ytick   = '-';
                ybottom = 'L';
                
                % box axes
                xtoptick     = '|';
                
                ytop         = 'T'; 
                yrighttop    = 'T';
                yrightbottom = 'J';
                yrighttick   = '-';                

            end

            % x-axis
            % ------------------
            xax  = repmat(xnorm, 1,xlen);
            xtop = xax;

            % range in x-values per tick
            xincr = (obj.xmax-obj.xmin)/xlen;

xtick_locations = 2:13:xlen;

            xax (xtick_locations) = xtick;
            xtop(xtick_locations) = xtoptick;
            

            % y-axis
            % ------------------                
            vspace = repmat(' ', size(xlbl,1)+1,1);                
            yax    = [repmat(ynorm, ylen-1,1); ybottom; vspace];            
            yright = yax;
            if obj.box
                yax(1) = ytop; end            
            
            % range in y-values per tick
            yincr = (obj.ymax-obj.ymin)/ylen;
            
ytick_locations = 2:3:ylen-1;
            
            yax   (ytick_locations) = ytick;
            yright(ytick_locations) = yrighttick;
            
            yright(1)   = yrighttop;
            yright(yright == ybottom) = yrightbottom;

            % both axes
            % ------------------
            if ~isempty(xlbl)
                xax = char(xax, xlbl); end
            if ~isempty(ylbl)
                ylbl = char(ylbl,vspace);
                yax = [ylbl yax]; 
            end 
            
        end
        
        %% class management
        
        % construct display string similar to how structs are displayed
        function str = construct_struct_display(obj)
            % TODO
        end
        
    end
    
end

    
    