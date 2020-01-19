function [] = plotDataPoints(dataPointsIn, tol, a, m, k, alpha, h, varargin)
%Fuction 2 plot data points
%   %plotDataPoints(points2End, 100*tol, a, m, k, alpha, h)
%
%uses cursorUpdate function to add annotations, stored in DataPoint.comm
%on plot, with data cursor mode

%%��������� �������
fig=[];%figure('Name','�������� ��������','Color', 'white');

%���� ���� ������� � �������������� ������� ����������
for i=1:length(varargin)
    if isa(varargin{i},'matlab.ui.Figure')
        fig=varargin{i};
    end
end

if isempty(fig)
    fig=figure;
end

fig.set('Name', '�������� ��������')
fig.set('Color', 'white')

axes1=findall(fig,'type','axes');

if isempty(axes1)
    axes1 = axes('Parent',fig);
end

hold(axes1,'on');
box(axes1,'on');

xlabel('x, ���')
ylabel('y, ���/�')
set(axes1,'FontName','Times New Roman');
hold on
grid minor

%��������� �������������� ���������� �� ������
addInfo=cell(length(dataPointsIn),3);
for i=1:length(dataPointsIn)
    point=dataPointsIn(1,i);
    addInfo(i,:)={[point.x point.y],point.comm, point.id};
end

dcm_obj = datacursormode(fig);
datacursormode on
set(dcm_obj,'UpdateFcn',{@cursorUpdate,addInfo})

%%����������

L1x=[];L1y=[];
L2x=[];L2y=[];
L3x=[];L3y=[];
L4x=[];L4y=[];

%���������� ������� ����������

%����� ������� ������ ������
xLims=[dataPointsIn(1,end).x dataPointsIn(1,end).x];
yLims=[dataPointsIn(1,end).y dataPointsIn(1,end).y];

for i=1:length(dataPointsIn)
    point=dataPointsIn(i);
    if(point.x>xLims(end))
        xLims(end)=point.x;
    else
        if(point.x<xLims(1))
            xLims(1)=point.x;
        end
    end
    
    if(point.y>yLims(end))
        yLims(end)=point.y;
    else
        if(point.y<yLims(1))
            yLims(1)=point.y;
        end
    end
end

L1x=xLims;
L1y=(alpha-L1x)/k;

L2x=xLims;
L2y=(alpha-h-L2x)/k;

L4x=xLims;
L4y=(-(alpha-h)-L4x)/k;

L3x=xLims;
L3y=(-alpha-L3x)/k;

plot(L1x,L1y, 'Color', [0.8500 0.3250 0.0980]);%��������� ��� ����� ���������
plot(L2x,L2y, 'Color', [0 0.4470 0.7410]);%����� ��� ����� ����������
plot(L4x,L4y, 'Color', [0 0.4470 0.7410]);
plot(L3x,L3y, 'Color', [0.8500 0.3250 0.0980]);

count=0;


while(count<length(dataPointsIn))%(~isempty(dataPointsIn))
    point=dataPointsIn(1,end-count);%����������� ����� � �����
    %dataPointsIn=dataPointsIn(:,1:end-1);%������� � �� ������� �� ����������
    count=count+1;
    
    prevPoint=point;%���������� ������ ����������� �������
    
    found=0;%���� �������
    %���� ������ ���������� ����� �� id
    for i=(length(dataPointsIn)-count):-1:1%���� �� ���������� ������, � �������
        if(point.prevId==dataPointsIn(i).id)%���� ����� �������� ����������
            prevPoint=dataPointsIn(i);%����������� � ������� �� ����� ������
            found=1;
            break;
        end
    end
    
    if(~found)%���� � ���� ������� �� ����� - ���� � ����������
        for i=length(dataPointsIn):-1:(length(dataPointsIn)-count)%���� �� ���������� ������, � �������
            if(point.prevId==dataPointsIn(i).id)%���� ����� �������� ����������
                prevPoint=dataPointsIn(i);%����������� � ������� �� ����� ������
                break;
            end
        end
    end
    
    %���� f � ���������� ����� ���������� ��� ����������
    if(prevPoint.f~=0)%���� f=1,-1 - �������� �������
        
        %��� �������� �������� �������� ����������� x(y), �.�. ����������
        %����� ���: (��������� y^2�ax=c), + ��� f=1, - ��� f=-1
        %     ^y
        %    *|-_   -yEnd
        %     |  \  
        %-----|---|->x     
        %     |  /  
        %     | *    -yStart
        
        c=point.y^2+prevPoint.f*a*point.x;%���������� � ���������
        
        yStart=min(point.y, prevPoint.y);
        yEnd=max(point.y, prevPoint.y);
        %stepsNum=ceil(abs(yStart-yEnd)/tol);
        
        y=[yStart:tol:yEnd yEnd];
        x=(prevPoint.f*(c-y.^2))/a;
        
        plot(x,y,'r');
        
        plot(point.x,point.y,'xr');%������� ������� ����� ������� �������
        
        continue
    end
    %f=0 - ��������� �������
    %��� ��������� �������� �������� ����������� y(x), �.�. ����������
    %����� ���: (��������� y^2 + m*x^2=�)
    %     ^y
    %   _-|-_   
    %  /  |  *
    %-*---|----->x
    %     |     
    % |      |
    %xStart  xEnd
    c=point.y^2+m*point.x^2;
    
    %� ������, ����� ��� ����� ����� �� ���� ������� �� ��� x
    if sign(atan2(point.y, point.x))==sign(atan2(prevPoint.y, prevPoint.x))
        % ������ ����������� y(x) ��� ������ ��������� x
        xStart=min(point.x, prevPoint.x);
        xEnd=max(point.x, prevPoint.x);
        
        x=[xStart:tol:xEnd xEnd];
        y1=[];
        %�����, ������� ��� ������ �������������
        if(point.y>0)%������ ����� � �������
            y=sqrt(c-m*x.^2);
        else
            if(point.y<0)%������ ����� � ������
                y=-sqrt(c-m*x.^2);
            else%������ ����� �� ����
                if(prevPoint.y>0)%������ ����� � �������
                    y=sqrt(c-m*x.^2);
                else
                    if(prevPoint.y<0)%������ ����� � ������
                        y=-sqrt(c-m*x.^2);
                    else%��� ����� �� ����, ������ ��� ��������
                        y=sqrt(c-m*x.^2);
                        y1=-y;
                    end
                end
            end
        end
        
        if(~isempty(y1))
            plot(x,y, '--b')
            plot(x,y1, '--b')
        else
            plot(x,y, 'b')
        end
    else
        %���� ���� ����������� ��� � - �������� 2 �������
        a1=sqrt(c/m);%����� ������� �
        
        if sign(atan2(prevPoint.y, prevPoint.x))>0
            %��������:
            %     ^y
            %   _-|-_
            %  /  |  *    -prevPoint
            %-|---|----->x
            %  \  |
            %   *         -point
            xPoz=[-a1:tol:prevPoint.x prevPoint.x];%������� ���� ��� �
            xNeg=[-a1:tol:point.x point.x];%������� ���� ��� �
        else
            %��������:
            %     ^y
            %     | *     -point
            %     |  \   
            %-----|---|->x
            %     | _/
            %   *-|-      -prevPoint
            xPoz=[point.x:tol:a1 a1];%������� ���� ��� �
            xNeg=[prevPoint.x:tol:a1 a1];%������� ���� ��� �
        end
        yPoz=sqrt(c-m*xPoz.^2);
        yNeg=-sqrt(c-m*xNeg.^2);
        
        plot(xPoz,yPoz, 'b')
        plot(xNeg,yNeg, 'b')
    end

    plot(point.x,point.y,'+b');%������� ������� ����� ������� ������
    
end
