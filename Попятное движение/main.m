%��������� ��� ����������� ���������� ��������� �������� ������� (���������
%�������� ������� � ������������).
%�� ����� ������� ������ ��������� �������� points2Start, ��������� ��
%����� ������ DataPoint. �� ������ - ������ points2End ��� �� ���������, �
%������� �������� ��� ����� ������������, ������� ������� ����� ���������.
%
%������ ����� ����� ���������� ������������� (DataPoint.id), � �����
%������������� ���������� ����� (DataPoint.prevId), � ������� ������� �����
%���� ������������� ������������������ �����. ��� ��������� �����
%�������������� ���������� � 0 � ���� �� �������.
%
%��� ������������ ������������� ����������� ������������ �������
%plotDataPoints
%
%����� ������������:
%\ \    \ \
% \ \    \ \
%  \ \    \ \
%  L3 L4  L2 L1

xLims=[-pi pi];%���
yLims=[-1 1];%���/�

%��������� �� �����������
xLims=sort(xLims);
yLims=sort(yLims);

points2End=[];%�������� ������
maxSteps=250;%����������� �� ����� ������������ �����
count=0;

a=1;%������������� ������������ �����������, ���/�2
m=0.1;%������������� ��������������� �������, ���/�2
k=0.5;%����. ��� �������� (x+ky)
alpha=deg2rad(1);%�������� ������� �����������, ���
h=deg2rad(0.1);%������ ����� �����������, ���
tol=deg2rad(0.005);%����������� �������

%��������� �����:

%{
%������ ��������� �����
p1=DataPoint(alpha-h-k*(0.1),...%x
    0.1,...%y
    1,...%f
    2,...%l num
    uint32(0),...%id
    uint32(0),...%prev id
    '��������� ����� 1');fig=figure;
%������ ��������� �����
%p2=DataPoint(alpha-h-k*(-0.1),-0.1,1,2,uint32(1),uint32(1), '��������� ����� 2');
points2Start=[p1];%points2Start=[p1 p2];

%points2Start=[p1];%������ ��� ������������ 
%}

[points2Start,fig]=getPointsPassiveMove(tol, m, k, alpha, h);

id=uint32(length(points2Start));

%���������� �� ������� ����� �� �������
simplePlot=false(1);%true(1);
%���� ��� � ��������������� �������, �������� ��������� �������� �������� plotDataPoints

%��������� ����������
syms x y

if(simplePlot) plot([],[]); end
if(simplePlot) hold('on'); end

while(count<maxSteps)&&(~isempty(points2Start))
          
    count=count+1;
    point = points2Start(1,1);%���� ������ �����
    points2Start=removeElementByIndex(points2Start,1);%������� � �� ������� �� ���������
        
    %���� ����� �� ����� � �������� �������� y
    if (point.y>yLims(end))||(point.y<yLims(1))
        %���������� �����
        warning(strcat("����� [",num2str(point.x,4), " ", num2str(point.y,4),...
            "] ����� ��� �������� ���������� y:[",...
            num2str(yLims(1),4), ":", num2str(yLims(end),4),...
            "]\n�� ������������ ��� �����"),...
            "")
        
        if isempty(point.comm)
            point.comm='����� ��� �������� y';
        else
            point.comm=char(string(point.comm) + newline +"����� ��� �������� y");
        end
        
        points2End=[points2End point];
        continue
    end
    
    %�������������� ���������� �����
    minLenSq=tol^2;
    removePoint=false(1);
    for i=1:length(points2End)
        %���������� ����� � ���������� ��������� f
        if point.f~=points2End(i).f
            continue
        end
        
        len=(k*(point.y-points2End(i).y))^2 + (point.x - points2End(i).x)^2;
        if len<=minLenSq
            %���������� �����
            removePoint=true(1);
                break
        end
    end
    
    if removePoint
        warning(strcat("����� ������� ������ � ��� ������������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
        continue
    end
    
    if (point.lNum == 1)%������ �� L1
        
        if (point.f == -1)
            %���������� �����
            warning(strcat("����� ������ �� L1 �� ���������� F=-1\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
        end
        
        
        %��� �������� ����������� �����������
        if (point.f == 1)
            
            % ���� ����������� � ������� ��������
            eqns = [y==yLims(end), y^2 + a*x==point.y^2+a*point.x];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %���� �� ����� - ���������� �����
            if(isempty(sol_x))
                warning(strcat("�� ���������� ������� � ������� ������� �� L1\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
            end
            
            
            while(~isempty(sol_x))%���� ���� �������
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %���� ����� ���� ��� y
                if( yd<0 )
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
                p2add=DataPoint(xd,yd,1,0,id,point.id, "����������� � y_{max}");
                id=id+1;%����������� ������� id
                points2End=[points2End p2add];%��������� ����� ����� � ������������
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %������� ����� ������������
            points2End=[points2End point];
            
            continue
        end
        
        %��� ��������� ����������� �����������
        if(point.f==0)
            
            %���� ����������� � ������ ���������� L4
            eqns = [x+k*y==-(alpha-h), y^2 + m*x^2==point.y^2 + m*point.x^2];
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
                
                %���� ������� ������ � ��������� �����
                if( (xd-point.x)^2 + (yd-point.y)^2<tol^2)
                    %���������� �����
                    continue
                end
                
                %���� ���� ������� �������
                if(yd>yLims(end))
                    %���������� �����
                    continue
                end
                
                %���� ����� ����� ������� - ���� �� ����������
%                 if(xd<xLims(1))
%                     %���������� �����
%                     continue
%                 end
                
                %���� ���� ��� y
                if(yd<0)
                    %���������� �����
                    continue
                    %������, ��������, ���� ��������
                end
                
                %���� ���� y ��������
                if sign(yd)~=sign(point.y)
                    %���������� �����
                    continue
                    %��. ��������
                end
                
                %��������� �����
                %� �������� ������������, F=-1, 4-� �����
                p2add=DataPoint(xd,yd,-1,4, id, point.id);
                id=id+1;%����������� ������� id
                points2Start=[points2Start p2add];%��������� ����� � ������� ���������
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %������� ����� ������������
            points2End=[points2End point];
            
            continue
        end
    end
    
    if (point.lNum ==2)%������ �� L2
        %��� f!=1 ����������� �����������
        if (point.f ~= 1)
            %���������� �����
            warning(strcat("����� ������ �� L2 �� �� ���������� F=1\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        % ���� ����������� � ������ ��������� L1
        eqns = [x+k*y==alpha, y^2 + a*x==point.y^2+a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ����� - ���������� �����
        if(isempty(sol_x))
            warning(strcat("�� ���������� ������� � �� L1 �� L2\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %������� ������������ �����
        
        iCount=length(sol_x);
        while iCount>0
            
            xd=double(sol_x(iCount));
            yd=double(sol_y(iCount));
            
            %���� ����� ���� ��������� ��� ���� �������� �������
            if (yd<point.y)||(yd>yLims(end))
                %������� �����
                sol_x=removeElementByIndex(sol_x,iCount);
                sol_y=removeElementByIndex(sol_y,iCount);
            end
            iCount=iCount-1;
        end
        
        %���� ��� ����� ��������� �������������
        if(isempty(sol_x))
            warning(strcat("��� �������� � �� L1 �� L2 ������������\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %���� �� ���� ���������� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %��������� �����
            
            %� �������� ������������, F=1, 1-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,1,1,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            %� �������� ������������, F=0, 1-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,0,1,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum ==4)%������ �� L4
        %��� f!=-1 ����������� �����������
        if (point.f ~= -1)
            %���������� �����
            warning(strcat("����� ������ �� L4 �� �� ���������� F=-1\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        % ���� ����������� � ������ ��������� L3
        eqns = [x+k*y==-alpha, y^2 - a*x==point.y^2 - a*point.x];
        [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
        
        %���� �� ����� - ���������� �����
        if(isempty(sol_x))
            warning(strcat("�� ���������� ������� � L3 �� L4\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %������� ������������ �����
        
        iCount=length(sol_x);
        while iCount>0
            
            xd=double(sol_x(iCount));
            yd=double(sol_y(iCount));
            
            %���� ����� ���� ��������� ��� ���� ������� �������
            if (yd>point.y)||(yd<yLims(1))
                %������� �����
                sol_x=removeElementByIndex(sol_x,iCount);
                sol_y=removeElementByIndex(sol_y,iCount);
            end
            iCount=iCount-1;
        end
        
        %���� ��� ����� ��������� �������������
        if(isempty(sol_x))
            warning(strcat("��� �������� � �� L4 �� L3 ������������\n��� ������������ ������\n���������� ����� [",...
                num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                "")
            continue
        end
        
        %���� �� ���� ���������� ��������
        while(~isempty(sol_x))%���� ���� �������
            
            xd=double(sol_x(end));
            yd=double(sol_y(end));
            
            sol_x=sol_x(1:end-1);
            sol_y=sol_y(1:end-1);
            
            %��������� �����
            
            %� �������� ������������, F=-1, 3-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,-1,3,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            %� �������� ������������, F=0, 3-� �����, id
            %- ������� ��������, prevId - id ������� �����
            p2add=DataPoint(xd,yd,0,3,id,point.id);
            id=id+1;%����������� ������� id
            points2Start=[points2Start p2add];%��������� ����� � ������� ���������
            
            if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            
        end
        
        %������� ����� ������������
        points2End=[points2End point];
        
        continue
    end
    
    if (point.lNum ==3)%������ �� L3
        if (point.f == 1)
            %���������� �����
            warning(strcat("����� ������ �� L3 �� ���������� F=1\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
        end
        
        
        %��� �������� ����������� �����������
        if (point.f == -1)
            
            % ���� ����������� � ������ ��������
            eqns = [y==yLims(1), y^2 - a*x==point.y^2 - a*point.x];
            [sol_x, sol_y]=vpasolve(eqns, [x y], [point.x point.y]);
            
            %���� �� ����� - ���������� �����
            if(isempty(sol_x))
                warning(strcat("�� ���������� ������� � ������ ������� �� L3\n��� ������������ ������\n���������� ����� [",...
                    num2str(point.x,4), " ", num2str(point.y,4), "]\n"),...
                    "")
                continue
            end
            
            
            while(~isempty(sol_x))%���� ���� �������
                
                xd=double(sol_x(end));
                yd=double(sol_y(end));
                
                sol_x=sol_x(1:end-1);
                sol_y=sol_y(1:end-1);
                
                %���� ����� ���� ��� y
                if( yd>0 )
                    %���������� �����
                    continue
                end
                
                %���� ����� ������ ����� L3
                if xd+k*yd>-alpha
                    %���������� �����
                    continue
                end
                
                %��������� �������� �����
                %� �������� ������������, F=-1, 0-� ����� (�� �� �����), id
                %- ������� ��������, prevId - id ������� �����
                p2add=DataPoint(xd,yd,-1,0,id,point.id,"����������� � y_{min}");
                id=id+1;%����������� ������� id
                points2End=[points2End p2add];%��������� ����� ����� � ������������
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %������� ����� ������������
            points2End=[points2End point];
            
            continue
        end
        
        %��� ��������� ����������� �����������
        if(point.f==0)
            
            %���� ����������� � ������ ���������� L2
            eqns = [x+k*y==alpha-h, y^2 + m*x^2==point.y^2 + m*point.x^2];
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
                
                %���� ������� ������ � ��������� �����
                if( (xd-point.x)^2 + (yd-point.y)^2<tol^2)
                    %���������� �����
                    continue
                end
                
                %���� ���� ������ �������
                if(yd<yLims(1))
                    %���������� �����
                    continue
                end
                
                %���� ������ ������ ������� - ���� �� ����������
%                 if(xd>xLims(end))
%                     %���������� �����
%                     continue
%                 end
                
                %���� ���� ��� y
                if(yd>0)
                    %���������� �����
                    continue
                    %������, ��������, ���� ��������
                end
                
                %���� ���� y ��������
                if sign(yd)~=sign(point.y)
                    %���������� �����
                    continue
                    %��. ��������
                end
                
                %��������� �����
                %� �������� ������������, F=1, 2-� �����
                p2add=DataPoint(xd,yd,1,2, id, point.id);
                id=id+1;%����������� ������� id
                points2Start=[points2Start p2add];%��������� ����� � ������� ���������
                
                if(simplePlot) plot([point.x xd], [point.y yd], '-o'); end
            end
            
            %������� ����� ������������
            points2End=[points2End point];
            
            continue
        end
    end
end
if count>=maxSteps
    warning(strcat("��������� ��������� �� ���������� ",num2str(count),...
        " �����\n������� ������� ����� ���� ��������\n"),...
                "")
end
plotDataPoints([points2End points2Start], 10*tol, a, m, k, alpha, h, fig)
