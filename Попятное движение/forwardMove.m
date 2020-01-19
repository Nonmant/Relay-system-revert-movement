function [points2End, points2Start, count] = forwardMove(dataPointsIn, tol, a, m, k, alpha, h, maxSteps, yLims, varargin)
%������� ��� ������� �������� �������
%   dataPointsIn -  ������ ����� ������
%   tol - ����������� ���������� ����� ����������������� ���������
%   a - ������������� ������������ �����������, ���/�2
%   m - ������������� ��������������� �������, ���/�2
%   k - ����. ��� �������� (x+ky)
%   alpha - �������� ������� �����������, ���
%   h - ������ ����� �����������, ���
%   maxSteps - ������������ ����� �����
%
%   func2end - ������� (x, y), ��� �������� �������� <tol ������������
%   �����
%
%   points2End - ������������ �����
%   points2Start - �������������� �����
%   count - ����������� ����� �����

func2end=@(x,y)10*tol;

quiteMode=false(1);

%���� ������� ��������� � �������������� ������� ����������
for i=1:length(varargin)
    if isa(varargin{i},'function_handle')
        if nargin(varargin{i})==2 %���� 2 ������� ���������
            func2end=varargin{i}; %�������� ����������� �������
        end
    end
    
    if isa(varargin{i},'string')||isa(varargin{i},'char')
        switch varargin{i}
            case 'quiteMode'
                quiteMode=true(1);
        end
    end
    
end

points2Start=dataPointsIn;
points2End=[];
count=0;

id=uint32(length(points2Start));

yLims=sort(yLims);

%��������� ����������
syms x y
solFound=false(1);

%��������� ����� ������������
LEqns=[...
    x+k*y==alpha,...%L1
    x+k*y==alpha-h,...%L2
    x+k*y==-alpha,...%L3
    x+k*y==-(alpha-h)...%L4
    ];
while(count<maxSteps)&&(~isempty(points2Start))
    count=count+1;
    point = points2Start(1,1);%���� ������ �����
    points2Start=removeElementByIndex(points2Start,1);%������� � �� ������� �� ���������
    
    if func2end(point.x, point.y) < tol %������� ������
        if ~quiteMode
            disp("����� �� ������� " + func2str(func2end) + newline)
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        continue
    end
    
    %���� ����� �� ����� � �������� �������� y
    if (point.y>yLims(end))||(point.y<yLims(1))
        %���������� �����
        if ~quiteMode
            warning(strcat("����� [",num2str(point.x,4), " ", num2str(point.y,4),...
                "] ����� ��� �������� ���������� y:[",...
                num2str(yLims(1),4), ":", num2str(yLims(end),4),...
                "]\n�� ������������ ��� �����"),...
                "")
        end
        if isempty(point.comm)
            point.comm='����� ��� �������� y';
        else
            point.comm=char(string(point.comm) + newline +"����� ��� �������� y");
        end
        
        points2End=[points2End point];
        continue
    end
    
    if (point.lNum == 1)%������ �� L1
        %��� �������� ����������� �����������
        
        solFound=false(1);%���� ���������� ����������� � L2
        
        % ���� ����������� � L2
        eqns = [LEqns(2), y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        while(~isempty(sol_x))%���� ���� �������
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ��������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %���� ����� ���� ������ �������
            if yd<yLims(1)
                %���������� �����
                continue
            end
            
            %��������� �����
            %� �������� ������������, F=1, 2-� �����
            p2add=DataPoint(xd,yd,0,2, id, point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            solFound=true(1);
        end
        
        % ���� ���� ������� ����������� � L2 - �� ����� ������ �����������
        % � ������ ��������
        if(solFound)
            %������� ����� ������������
            points2End=[points2End point];
            continue
        end
        
        % ���� ����������� � ������ ��������
        eqns = [y==yLims(1), y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ����� - ���������� �����
        if(isempty(sol_x))
            if ~quiteMode
                warning(strcat("�� ���������� ������� � L1 �� ������ �������\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
            end
            continue
        end
        
        
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ��������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %���� ����� ����� ����� L1
            if xd+k*yd<alpha
                %���������� �����
                continue
            end
            
            %��������� �������� �����
            %� �������� ������������, F=1, 0-� ����� (�� �� �����), id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,1,0,id,point.id, "����������� � y_{min}");
            id=id+1;%����������� ������� id
            points2End=[points2End p2add];%��������� ����� ����� � ������������
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
        
    end
    
    if (point.lNum ==2)%������ �� L2
        %��� ��������� ����������� �����������
        
        solFound=false(1);%���� ���������� ����������� � L1
        
        % ���� ����������� � ������ ��������� L1
        eqns = [LEqns(1), y^2 + m*x^2==point.y^2 + m*point.x^2];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ���� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ���������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %��������� �����
            
            %� �������� ������������, F=0, 1-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            solFound=true(1);
        end
        
        % ���� ���� ������� ����������� � L1 - �� ����� ������ �����������
        % � L3
        if(solFound)
            %������� ����� ������������
            points2End=[points2End point];
            continue
        end
        
        % ���� ����������� � ������ ��������� L3
        eqns = [LEqns(3), y^2 + m*x^2==point.y^2 + m*point.x^2];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ���� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ���������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %��������� �����
            
            %� �������� ������������, F=-1, 3-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum == 3)%������ �� L3
        %��� �������� ����������� �����������
        
        solFound=false(1);%���� ���������� ����������� � L4
        
        % ���� ����������� � L4
        eqns = [LEqns(4), y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        while(~isempty(sol_x))%���� ���� �������
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ��������
            if yd<point.y
                %���������� �����
                continue
            end
            
            %���� ����� ���� ������� �������
            if yd>yLims(end)
                %���������� �����
                continue
            end
            
            %��������� �����
            %� �������� ������������, F=-1, 4-� �����
            p2add=DataPoint(xd,yd,0,4, id, point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            solFound=true(1);
        end
        
        % ���� ���� ������� ����������� � L4 - �� ����� ������ �����������
        % � ������� ��������
        if(solFound)
            %������� ����� ������������
            points2End=[points2End point];
            continue
        end
        
        % ���� ����������� � ������� ��������
        eqns = [y==yLims(end), y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ����� - ���������� �����
        if(isempty(sol_x))
            if ~quiteMode
                warning(strcat("�� ���������� ������� � L3 �� ������� �������\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
            end
            continue
        end
        
        
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ��������
            if yd<point.y
                %���������� �����
                continue
            end
            
            %���� ����� ������ ����� L4
            if xd+k*yd>-(alpha-h)
                %���������� �����
                continue
            end
            
            %��������� �������� �����
            %� �������� ������������, F=-1, 0-� ����� (�� �� �����), id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,-1,0,id,point.id, "����������� � y_{max}");
            id=id+1;%����������� ������� id
            points2End=[points2End p2add];%��������� ����� ����� � ������������
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
        
    end
    
    if (point.lNum ==4)%������ �� L4
        %��� ��������� ����������� �����������
        
        solFound=false(1);%���� ���������� ����������� � L3
        
        % ���� ����������� � ������ ��������� L3
        eqns = [LEqns(3), y^2 + m*x^2==point.y^2 + m*point.x^2];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ���� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ���������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %��������� �����
            
            %� �������� ������������, F=0, 3-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            solFound=true(1);
        end
        
        % ���� ���� ������� ����������� � L3 - �� ����� ������ �����������
        % � L1
        if(solFound)
            %������� ����� ������������
            points2End=[points2End point];
            continue
        end
        
        % ���� ����������� � ������ ��������� L1
        eqns = [LEqns(1), y^2 + m*x^2==point.y^2 + m*point.x^2];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ���� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %���� ������� �����������
            if(~isreal(xd))||(~isreal(yd))
                %���������� �����
                continue
            end
            
            %���� ����� ���� ���������
            if yd>point.y
                %���������� �����
                continue
            end
            
            %��������� �����
            
            %� �������� ������������, F=0, 1-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
    end
    
end



