%% PLOT THE TMSMultiLab / TMS@40 network map...
hubs=NaN(2,6);                                                               % main test sites: lat, long
hubs(:,1)= [52.449446, -1.930496];                                           % University of Birmingham = Command & Control
hubs(:,2)= [51.296819,  1.063497];                                           % University of Kent = South
hubs(:,3)= [52.938871, -1.198138];                                           % University of Nottingham = Midlands
hubs(:,4)= [54.010541, -2.786361];                                           % Lancaster University = North
hubs(:,5)= [53.229616, -4.129268];                                           % Bangor University = Wales
%hubs(:,5)=[50.000940,  8.257245];                                           % Mainz = Europe

%% UK & IRELAND______________________________________________________________
sats=NaN(4,60);                                                             % satellite TMS centres (lat, long, N, Hub to connect to)
sats(:,1)= [52.938871, -1.198138,9,NaN];                                    % Nottingham (9)
sats(:,2)= [52.449446, -1.930496,7,NaN];                                    % birmingham (7)
sats(:,3)= [52.205474,  0.114251,6,NaN];                                    % Cambridge (6)
sats(:,4)= [51.296819,  1.063497,3,NaN];                                    % Kent (3)
sats(:,5)= [53.466913, -2.233304,3,NaN];                                    % Manchester, MMU (3)
sats(:,6)= [51.501176, -2.545384,3,NaN];                                    % UWE, Bristol (3)
sats(:,7)= [55.872564, -4.290076,2,NaN];                                    % Glasgow (2)
sats(:,8)= [54.010541, -2.786361,2,NaN];                                    % Lancaster (2)
sats(:,9)= [51.754909, -1.253605,2,NaN];                                    % Oxford (2)
sats(:,10)=[53.405139, -2.963729,2,NaN];                                    % Liverpool (2)
sats(:,11)=[53.381503, -1.487746,1,NaN];                                    % Sheffield (1)
sats(:,12)=[57.164847, -2.100603,1,NaN];                                    % Aberdeen (1)
sats(:,13)=[53.229616, -4.129268,1,NaN];                                    % Bangor (1)
sats(:,14)=[50.737218, -3.534460,1,NaN];                                    % Exeter (1)
sats(:,15)=[51.523395, -0.132872,1,NaN];                                    % London: Brunel, Middlesex, LMU, QMUL, UEL, Roehampton, Goldsmiths, South Bank, Greenwich, City, KCL, ICL, UCL/ICN/IoN (1)
sats(:,16)=[54.979297, -1.614167,1,NaN];                                    % Northumbria, Newcastle (1)
sats(:,17)=[51.441500, -0.941075,1,NaN];                                    % Reading (1)
sats(:,18)=[56.146169, -3.917071,1,NaN];                                    % Stirling (1)
sats(:,19)=[53.946184, -1.051063,1,NaN];                                    % York (1)
sats(:,20)=[53.306514, -6.218243,1,NaN];                                    % Dublin (1)
sats(:,21)=[56.458309, -2.981037,1,NaN];                                    % Dundee
sats(:,22)=[51.782344, -0.061678,1,NaN];                                    % Herts
sats(:,23)= [51.378303,-2.325991,1,NaN];                                    % Bath

% EUROPE
sats(:,24)=[51.335591, 12.391673,2,NaN];                                    % MPI (2)
sats(:,25)=[45.523622, 10.202295,1,NaN];                                    % Brescia (1)
sats(:,26)=[49.992727,  8.261164,1,NaN];                                    % Mainz
sats(:,27)=[46.067068, 11.123677,1,NaN];                                    % Trento (1)
sats(:,28)=[51.046746,  3.728124,1,NaN];                                    % Ghent (1)

% WORLD
sats(:,29)=[37.427779,-122.16960,1,NaN];                                    % Stanford (1)

%sats(:,5)= [52.621243, -1.124117,1,NaN];                                    % Leicester
%sats(:,6)= [52.622123,  1.241748,1,3];                                      % East Anglia
%sats(:,8)= [50.933738, -1.393510,1,NaN];                                    % Southampton
%sats(:,10)=[51.487348, -3.176352,1,NaN];                                    % Cardiff
%sats(:,11)=[51.877884,  0.947957,1,NaN];                                    % Essex
%sats(:,14)=[51.425773, -0.562375,1,NaN];                                    % Royal Holloway
%sats(:,16)=[53.806782, -1.554002,1,NaN];                                    % Leeds
%sats(:,20)=[54.608469, -5.934983,1,NaN];                                    % Ulster
%sats(:,23)=[50.867800, -0.086894,1,NaN];                                    % Sussex
%sats(:,25)=[56.341796, -2.793580,1,NaN];                                    % St Andrews 

%sats(:,30)=[51.338998, 12.380758,1,5];                                      % Leipzig
%sats(:,31)=[47.993694,  7.846489,1,5];                                      % Freiburg
%sats(:,32)=[47.797274, 13.048677,2,5];                                      % Salzburg
%sats(:,33)=[50.669833,  4.616744,1,5];                                      % Louvain
%sats(:,34)=[55.680327, 12.573203,1,5];                                      % Copenhagen

%sats(:,46)=[-38.18823,144.363968,5,6];                                      % Deakin
%sats(:,47)=[-33.88823,151.188058,1,6];                                      % Sydney
%sats(:,48)=[ 31.23075,121.387320,1,6];                                      % East China
%sats(:,49)=[ 32.57217,-83.255147,1,6];                                      % Georgia
%sats(:,50)=[ 42.37703,-71.116149,1,6];                                      % Harvard
%sats(:,51)=[ 44.04492,-123.07206,1,6];                                      % Oregon
%sats(:,52)=[ 39.16828,-86.522949,1,6];                                      % Indiana

UK=1:23;
EU=24:28;
WO=29;

%% PLOT FIGURES WITH ALL THIS INFORMATION ON_________________________________
geoaxes;                                                                    % UK
hold on;
geoplot(hubs(1,1:5),hubs(2,1:5),'bo','MarkerSize',6,'MarkerFaceColor','b');
geoplot(sats(1,UK),sats(2,UK),'bo','MarkerSize',6,'MarkerFaceColor','b');
for s=UK                                                                    % for each satellite
    dist=sqrt((hubs(1,:)-sats(1,s)).^2 + (hubs(2,:)-sats(2,s)).^2);         % find the distance to all hubs
    if isfinite(sats(4,s))
        h=sats(4,s);                                                        % pre-selected hub?
    else
        [~,h]=min(dist);                                                    % find the closest hub
    end
    geoplot([hubs(1,h),sats(1,s)],[hubs(2,h),sats(2,s)],'b-','LineWidth',2);% plot a line between them
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_UK.png','-dpng');
close;

% REPEAT, WITHOUT HUB/SAT distinction_______________________________________
geoaxes;                                                                    % UK
hold on;
for s=UK
    geoplot(sats(1,s),sats(2,s),'bo','MarkerSize',sats(3,s).*4,'MarkerFaceColor','b');
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_UK2.png','-dpng');
close;

%% EUROPE___________________________________________________________________
geoaxes;
hold on;
geoplot(hubs(1,5),hubs(2,5),'bo','MarkerSize',6,'MarkerFaceColor','b');
geoplot(sats(1,EU),sats(2,EU),'bo','MarkerSize',6,'MarkerFaceColor','b');
for s=EU                                                                    % for each satellite
    dist=sqrt((hubs(1,:)-sats(1,s)).^2 + (hubs(2,:)-sats(2,s)).^2);         % find the distance to all hubs
    if isfinite(sats(4,s))
        h=sats(4,s);                                                        % pre-selected hub?
    else
        [~,h]=min(dist);                                                    % find the closest hub
    end
    geoplot([hubs(1,h),sats(1,s)],[hubs(2,h),sats(2,s)],'b-','LineWidth',2);% plot a line between them
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_Europe.png','-dpng');
close;

% REPEAT, WITHOUT HUB/SAT distinction_______________________________________
geoaxes;
hold on;
for s=EU
    geoplot(sats(1,s),sats(2,s),'bo','MarkerSize',sats(3,s).*4,'MarkerFaceColor','b');
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_Europe2.png','-dpng');
close;

%% WORLD____________________________________________________________________
geoaxes;
hold on;
geoplot(hubs(1,6),hubs(2,6),'bo','MarkerSize',6,'MarkerFaceColor','b');
geoplot(sats(1,WO),sats(2,WO),'bo','MarkerSize',6,'MarkerFaceColor','b');
for s=WO                                                                    % for each satellite
    dist=sqrt((hubs(1,:)-sats(1,s)).^2 + (hubs(2,:)-sats(2,s)).^2);         % find the distance to all hubs
    if isfinite(sats(4,s))
        h=sats(4,s);                                                        % pre-selected hub?
    else
        [~,h]=min(dist);                                                    % find the closest hub
    end
    geoplot([hubs(1,h),sats(1,s)],[hubs(2,h),sats(2,s)],'b-','LineWidth',2);% plot a line between them
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_World.png','-dpng');
close;

% REPEAT, WITHOUT HUB/SAT distinction_______________________________________
geoaxes;
hold on;
for s=WO
    geoplot(sats(1,s),sats(2,s),'bo','MarkerSize',sats(3,s).*4,'MarkerFaceColor','b');
end
set(gcf,'Position',[0,0,800,600]);
print('TMSMultiLab_Network_World2.png','-dpng');
close all;