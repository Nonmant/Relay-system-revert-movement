function output_txt = cursorUpdate(~, event_obj, varargin)
%dataPointClickedCallback prints summary of data point at data tip
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors
% third input is a cell formatted as below:
% (:,1)|(:,2)
% [x,y]|"text"
% [0,0]|"It's an origin point"
%
% initialise it like: info={[0,0],"text"} or
% info(2,:) = {[1,1],"text1"};
% info(3,:) = {[1,2],"text2"}
%
% example:
% dataY=rand(5,1);
% dataX=1:5;
% for i=1:5
% addInfo(i,:)={[dataX(i) dataY(i)],strcat('¹',num2str(i))}
% end

pos = event_obj.Position;


% Display the x and y values:
output_txt = {['X:',num2str(pos(1),'%.4g'),...
    ' Y:',num2str(pos(2),'%.4g')]};

if(nargin>2)%if there is extra data to show
    addInfo=varargin{1};
    
    for i=1:length(addInfo)
        if(isequal(addInfo{i,1},pos))
            output_txt{end+1} = ["id="+num2str(addInfo{i,3})];
            if(~isempty(addInfo{i,2}))
                output_txt{end+1} = [addInfo{i,2}];
            end
            break;
        end
    end
    
end
end

%%
%example:
% fig=figure;
% plot(rand(6,1))
% dcm_obj = datacursormode(fig);
% datacursormode on
% set(dcm_obj,'UpdateFcn',{@cursorUpdate,addInfo})

