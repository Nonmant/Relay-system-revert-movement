classdef DataPoint
    %Point of data to solve
    %   Initial coordinate (x,y) with params, to find previous value(s)
    
    properties
        x
        y
        f    %1,-1,0 - f той траектории, которая привела в эту точку!
        lNum %1,2,3,4
        id@uint32
        prevId@uint32
        
        comm@char=''%additional info, would be added to final plot within plotDataPoints function
    end
    
    methods
        function obj = DataPoint(xIn,yIn,fIn,lNumIn, idIn, prevIdIn, varargin)
            %Construct an instance of this class
            %   xIN,yIn,fIn,lNum, comment
            obj.x = xIn;
            obj.y = yIn;
            obj.f = fIn;
            obj.lNum = lNumIn;
            obj.id=idIn;
            obj.prevId=prevIdIn;
            if(nargin>6)
                obj.comm=varargin{1};
            end
        end
    end
end

