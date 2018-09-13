function drawlattice(Offset, Scaling, hAxes)
%DRAWLATTICE - Draws the AT lattice to a figure
%  drawlattice(Offset {0}, Scaling {1}, hAxes {gca})
%
%  Programmers Note: The AT index is stored in the Userdata of each symbol.
%
%  Written by Greg Portmann


global THERING


% Minimum icon width in meters
MinIconWidth = .09;

if nargin < 1
    Offset = 0;
end
Offset = Offset(1);
if nargin < 2
    Scaling = 1;
end
Scaling = Scaling(1);

if nargin < 3
    hAxes = gca;
end

SPositions = findspos(THERING, 1:length(THERING)+1);
L = SPositions(end);
plot(hAxes, [0 L], [0 0]+Offset, 'k');

% Remember the hold state then turn hold on
HoldState = ishold(hAxes);
hold(hAxes, 'on');


ATIndexHCM = family2atindex(gethcmfamily);
ATIndexVCM = family2atindex(getvcmfamily);


% Make default icons for elements of different physical types
for i = 1:length(THERING)
    NumberOfFinds = 0;
    
    SPos = SPositions(i);
    if isfield(THERING{i},'BendingAngle') && THERING{i}.BendingAngle
        % make icons for bending magnets
        NumberOfFinds = NumberOfFinds + 1;
        IconHeight = .3;
        IconColor = [1 1 0];
        IconWidth = THERING{i}.Length;
        if IconWidth < MinIconWidth % meters
            IconWidth = MinIconWidth;
            SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
        end
        vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
        vy = [IconHeight IconHeight -IconHeight -IconHeight];
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-', 'Parent',hAxes);
        set(h, 'UserData', i);

        %if IconWidth < .1 % meters
        %    set(h, 'EdgeColor', IconColor);
        %end

    elseif isfield(THERING{i},'K') && THERING{i}.K
        % Quadrupole
        NumberOfFinds = NumberOfFinds + 1;
       if THERING{i}.K >= 0
            % Focusing quadrupole           
            IconHeight = .8;
            IconColor = [1 0 0];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos SPos+IconWidth/2  SPos+IconWidth SPos+IconWidth/2 SPos];
            vy = [0          IconHeight               0      -IconHeight    0];
        else
            % Defocusing quadrupole
            IconHeight = .7;
            IconColor = [0 0 1];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                % Center about starting point
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos+.4*IconWidth    SPos    SPos+IconWidth  SPos+.6*IconWidth  SPos+IconWidth    SPos      SPos+.4*IconWidth];
            vy = [     0            IconHeight   IconHeight          0              -IconHeight  -IconHeight    0];
        end
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-', 'Parent',hAxes);
        set(h, 'UserData', i);
        %if IconWidth < .1 % meters
        %    set(h, 'EdgeColor', IconColor);
        %end

    elseif isfield(THERING{i},'PolynomB') && length(THERING{i}.PolynomB)>2 && (THERING{i}.PolynomB(3) || any(strcmpi(THERING{i}.FamName,{'SF','SFF','SD','SDD'})))
        % Sextupole
        NumberOfFinds = NumberOfFinds + 1;
        if THERING{i}.PolynomB(3)>0 || any(strcmpi(THERING{i}.FamName,{'SF','SFF'}))
            % Focusing sextupole
            IconHeight = .6;
            IconColor = [1 0 1];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos          SPos+.33*IconWidth  SPos+.66*IconWidth  SPos+IconWidth   SPos+IconWidth   SPos+.66*IconWidth  SPos+.33*IconWidth      SPos          SPos];
            vy = [IconHeight/3      IconHeight          IconHeight        IconHeight/3    -IconHeight/3      -IconHeight          -IconHeight     -IconHeight/3  IconHeight/3];
        else
            % Defocusing sextupole
            IconHeight = .6;
            IconColor = [0 1 0];
            IconWidth = THERING{i}.Length;
            if IconWidth < MinIconWidth % meters
                IconWidth = MinIconWidth;
                SPos = SPos - IconWidth/2 + THERING{i}.Length/2;
            end
            vx = [SPos          SPos+.33*IconWidth  SPos+.66*IconWidth  SPos+IconWidth   SPos+IconWidth   SPos+.66*IconWidth  SPos+.33*IconWidth      SPos          SPos];
            vy = [IconHeight/3      IconHeight          IconHeight        IconHeight/3    -IconHeight/3      -IconHeight          -IconHeight     -IconHeight/3  IconHeight/3];
        end
        h = patch(vx, Scaling*vy+Offset, IconColor,'LineStyle','-', 'Parent',hAxes);
        set(h, 'UserData', i);
        %if IconWidth < .1 % meters
        %    set(h, 'EdgeColor', IconColor);
        %end

    elseif isfield(THERING{i},'Frequency') && isfield(THERING{i},'Voltage')
        % RF cavity
        NumberOfFinds = NumberOfFinds + 1;
        IconColor = [1 0.5 0];
        h = plot(hAxes, SPos, 0+Offset, 'o', 'MarkerFaceColor', IconColor, 'Color', IconColor, 'MarkerSize', 4);
        set(h, 'UserData', i);

    elseif strcmpi(THERING{i}.FamName,'BPM')
        % BPM
        NumberOfFinds = NumberOfFinds + 1;
        IconColor = 'k';
        h = plot(hAxes, SPos, 0+Offset, '.', 'Color', IconColor);
        %h = plot(hAxes, SPos, 0, 'o', 'MarkerFaceColor', IconColor, 'Color', IconColor, 'MarkerSize', 1.5)
        set(h, 'UserData', i);
        
    elseif strcmpi(THERING{i}.FamName,'TV')
        % TV screen
        NumberOfFinds = NumberOfFinds + 1;
        IconHeight = .7;
        IconColor = [.5 0 0];  %'k';
        %h = plot(hAxes, SPos, 0+Offset, 'x', 'Color', IconColor);
        h = plot(hAxes, SPos, Scaling*IconHeight+Offset, 'Marker','Square', 'MarkerFaceColor', IconColor, 'Color', IconColor, 'MarkerSize', 3.5);
        set(h, 'UserData', i);
    end
    
    % Since correctors could be a combined function magnet, test separately
    if any(strcmpi(THERING{i}.FamName,{'COR','XCOR','YCOR','HCOR','VCOR'})) || isfield(THERING{i},'KickAngle')
        % Corrector
        NumberOfFinds = NumberOfFinds + 1;
        
        if NumberOfFinds > 1
            IconWidth = 0;
        else
            IconWidth = THERING{i}.Length;
        end
        IconHeight = 1.0;  % was .8
        vx = [SPos   SPos];

        % Draw a line above for a HCM and below for a VCM
        % If it's not in the ML, then draw a line above and below
        CMFound = 1;
        if any(i == ATIndexVCM)
            IconColor = [0 0 0];
            vy = [-IconHeight 0];
            if IconWidth < MinIconWidth
                h = plot(hAxes, vx, Scaling*vy+Offset, 'Color', IconColor);
            else
                IconWidth = THERING{i}.Length;
                vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
                vy = [0 0 -IconHeight -IconHeight];
                %vy = [IconHeight IconHeight -IconHeight -IconHeight];
                h = patch(vx, Scaling*vy+Offset, IconColor, 'LineStyle', '-', 'Parent',hAxes);
                if IconWidth < MinIconWidth*2 % meters
                    set(h, 'EdgeColor', IconColor);
                end
            end
            set(h, 'UserData', i);
            CMFound = 0;
        end

        if any(i == ATIndexHCM)
            IconColor = [0 0 0];
            vy = [0 IconHeight];
            if IconWidth < MinIconWidth
                h = plot(hAxes, vx, Scaling*vy+Offset, 'Color', IconColor);
            else
                IconWidth = THERING{i}.Length;
                vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
                vy = [IconHeight IconHeight 0 0];
                %vy = [IconHeight IconHeight -IconHeight -IconHeight];
                h = patch(vx, Scaling*vy+Offset, IconColor, 'LineStyle', '-', 'Parent',hAxes);
                if IconWidth < MinIconWidth*2 % meters
                    set(h, 'EdgeColor', IconColor);
                end
            end
            set(h, 'UserData', i);
            CMFound = 0;
        end
        
        if CMFound
            IconColor = [0 0 0];
            vy = [-IconHeight IconHeight];
            if IconWidth < MinIconWidth
                h = plot(hAxes, vx, Scaling*vy+Offset, 'Color', IconColor);
            else
                IconWidth = THERING{i}.Length;
                vx = [SPos SPos+IconWidth SPos+IconWidth SPos];
                vy = [IconHeight IconHeight -IconHeight -IconHeight];
                h = patch(vx, Scaling*vy+Offset, IconColor, 'LineStyle', '-', 'Parent',hAxes);
                if IconWidth < MinIconWidth*2 % meters
                    set(h, 'EdgeColor', IconColor);
                end
            end
            set(h, 'UserData', i);
            CMFound = 0;
        end

    end
end


% Leave the hold state as it was at the start
if ~HoldState
    hold(hAxes, 'off');
end

a = axis(hAxes);
axis(hAxes, [0 L a(3:4)]);

