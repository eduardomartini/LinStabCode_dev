function [D1,D2,D1s,D1a,D2s,D2a] = Dmats_SBP(Nz,dz,order)
% Aaron Towne's original file with modified outputs
% Mohseni_pole function added as subfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- Define finite difference scheme with SBP closure ----------- %

if order == 2

    % Interior discretization
    a_i1 = [-1/2,0,1/2];                     % First derivative
    a_i2 = [1,-2,1];                         % Second derivative
    
    % Summation-by-parts closure
    a_b1 = [-1,1,0];                         % First derivative
    a_b2 = [1,-2,1];                         % Second derivative

elseif order == 4
    
    % Interior discretization
    a_i1 = [1/12,-2/3,0,2/3,-1/12];          % First derivative
    a_i2 = [-1/12,4/3,-5/2,4/3,-1/12];       % Second derivative
    
    % Summation-by-parts closure
    a_b1 = [-24/17,59/34,-4/17,-3/34,0,0;...
        -1/2,0,1/2,0,0,0;...
        4/43,-59/86,0,59/86,-4/43,0;...
        3/98,0,-59/98,0,32/49,-4/49];       % First derivative
    
    a_b2 = [2,-5,4,-1,0,0;
        1,-2,1,0,0,0;
        -4/43,59/43,-110/43,59/43,-4/43,0;
        -1/49,0,59/49,-118/49,64/49,-4/49]; % Second derivative
    
    
elseif order == 6    
    
    % Interior discretization
    a_i1 = [-1/60,3/20,-3/4,0,3/4,-3/20,1/60];
    a_i2 = [1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90];
    
    % Summation-by-parts closure
    a_b1 = [-21600/13649, 104009/54596, 30443/81894, -33311/27298, 16863/27298, -15025/163788, 0,0,0;
        -104009/240260,0,-311/72078,20229/24026,-24337/48052,36661/360390,0,0,0;
        -30443/162660,311/32532,0,-11155/16266,41287/32532,-21999/54220,0,0,0;
        33311/107180, -20229/21436,485/1398,0,4147/21436, 25427/321540,72/5359,0,0;
        -16863/78770,24337/31508,-41287/47262,-4147/15754,0,342523/472620,-1296/7877,144/7877,0;
        15025/525612,-36661/262806,21999/87602,-25427/262806,-342523/525612,0,32400/43801,-6480/43801,720/43801];
        
    a_b2 = [ 114170/40947, -438107/54596, 336409/40947, -276997/81894, 3747/13649, 21035/163788,0,0,0;
        6173/5860, -2066/879, 3283/1758, -303/293, 2111/3516, -601/4395,0,0,0;
        -52391/81330, 134603/32532, -21982/2711, 112915/16266, -46969/16266, 30409/54220, 0, 0, 0;
        68603/321540, -12423/10718, 112915/32154, -75934/16077, 53369/21436, -54899/160770, 48/5359, 0, 0;
        -7053/39385, 86551/94524, -46969/23631, 53369/15754, -87904/23631, 820271/472620, -1296/7877, 96/7877, 0;
        21035/525612, -24641/131403, 30409/87602, -54899/131403, 820271/525612, -117600/43801, 64800/43801, -6480/43801, 480/43801]; 


elseif order == 8
    
    % Interior discretization
    a_i1 = [1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280]; % First derivative
    a_i2 = [-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]; % Second derivative
    
    % Summation-by-parts closure
    a_b1 = [-2540160/1498139, 5544277/5992556, 198794991/29962780, -256916579/17977668, 20708767/1498139, -41004357/5992556, 27390659/17977668, -2323531/29962780,0,0,0,0;
        -5544277/31004596, 0, -85002381/22146140, 49607267/4429228, -165990199/13287684, 7655859/1107307, -7568311/4429228, 48319961/465068940, 0, 0,0,0;
        -66264997/8719620, 9444709/415220, 0, -20335981/249132, 32320879/249132, -35518713/415220, 2502774/103805, -3177073/1743924, 0,0,0,0;
        256916579/109619916, -49607267/5219996, 61007943/5219996, 0, -68748371/5219996, 65088123/5219996, -66558305/15659988, 3870214/9134993, 0,0,0,0;
        -20708767/2096689, 165990199/3594324, -96962637/1198108, 68748371/1198108, 0, -27294549/1198108, 14054993/1198108, -42678199/25160268, -2592/299527,0,0,0;
        13668119/8660148, -850651/103097, 35518713/2061940, -21696041/1237164, 9098183/1237164, 0, -231661/412388, 7120007/43300740, 3072/103097, -288/103097, 0,0;
        -27390659/56287644, 7568311/2680364, -22524966/3350455, 66558305/8041092, -14054993/2680364, 2084949/2680364, 0, 70710683/93812740, -145152/670091, 27648/670091, -2592/670091,0;
        2323531/102554780, -48319961/307664340, 9531219/20510956, -3870214/5127739, 2246221/3238572, -21360021/102554780, -70710683/102554780, 0, 4064256/5127739, -1016064/5127739, 193536/5127739, -18144/5127739];
    
    
    a_b2 = [4870382994799/1358976868290, -893640087518/75498714905, 926594825119/60398971924, -1315109406200/135897686829, 39126983272/15099742981, 12344491342/75498714905, -451560522577/2717953736580, 0,0,0,0,0;
        333806012194/390619153855, -154646272029/111605472530, 1168338040/33481641759, 82699112501/133926567036, -171562838/11160547253, -28244698346/167408208795, 11904122576/167408208795, -2598164715/312495323084,0,0,0,0;
        7838984095/52731029988, 1168338040/5649753213, -88747895/144865467, 423587231/627750357, -43205598281/22599012852, 4876378562/1883251071, -5124426509/3766502142, 10496900965/39548272491,0,0,0,0;
        -94978241528/828644350023, 82699112501/157837019052, 1270761693/13153084921, -167389605005/118377764289, 48242560214/39459254763, -31673996013/52612339684, 43556319241/118377764289, -44430275135/552429566682,0,0,0,0;
        1455067816/21132528431, -171562838/3018932633, -43205598281/36227191596, 48242560214/9056797899, -52276055645/6037865266, 57521587238/9056797899, -80321706377/36227191596, 8078087158/21132528431, -1296/299527,0,0,0;
        10881504334/327321118845, -28244698346/140280479505, 4876378562/9352031967, -10557998671/12469375956, 57521587238/28056095901, -278531401019/93520319670, 73790130002/46760159835, -137529995233/785570685228, 2048/103097,-144/103097,0,0;
        -135555328849/8509847458140, 11904122576/101307707835, -5124426509/13507694378, 43556319241/60784624701, -80321706377/81046166268, 73790130002/33769235945, -950494905688/303923123505, 239073018673/141830790969, -145152/670091, 18432/670091, -1296/670091,0;
        0, -2598164715/206729925524, 10496900965/155047444143, -44430275135/310094888286, 425162482/2720130599, -137529995233/620189776572, 239073018673/155047444143, -144648000000/51682481381, 8128512/5127739, -1016064/5127739, 129024/5127739, -9072/5127739];
    
    
else
    
    error('Only second, fourth, sixth, and eigth order FD are currently implemented');
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------ Get modified Mohseni & Colonius pole stencils ------------ %

% First derivative, symmetric
a_p1s = Mohseni_pole(a_i1,1);

% First derivative, anti-symmetric
a_p1a = Mohseni_pole(a_i1,-1);

% Second derivative, symmetric
a_p2s = Mohseni_pole(a_i2,1);

% Second derivative, anti-symmetric
a_p2a = Mohseni_pole(a_i2,-1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------- Construct differentiation matrixes ------------------ %


% Stencile sizes
[r,m] = size(a_b1);
w = (length(a_i1)-1)/2;


% Initialize matrixes
D1s = zeros(Nz);
D1a = zeros(Nz);
D2s = zeros(Nz);
D2a = zeros(Nz);

% Interior points
for jz = 1+w:Nz-r
    D1s(jz,jz-w:jz+w) = a_i1./dz;           % First derivative, symmetric
    D1a(jz,jz-w:jz+w) = a_i1./dz;           % First derivative, anti-sym
    D2s(jz,jz-w:jz+w) = a_i2./dz^2;         % Second derivative, symmetric
    D2a(jz,jz-w:jz+w) = a_i2./dz^2;         % Second derivative, anti-sym
end


% Plus boundary
for jz = Nz-r+1:Nz
    for jm = 1:m
        D1s(jz,Nz+1-jm) = (-a_b1(1+Nz-jz,jm))./dz;
        D1a(jz,Nz+1-jm) = (-a_b1(1+Nz-jz,jm))./dz;
        D2s(jz,Nz+1-jm) = (a_b2(1+Nz-jz,jm))./dz^2;
        D2a(jz,Nz+1-jm) = (a_b2(1+Nz-jz,jm))./dz^2;
    end
end

% matrices without pole condition
D1  = D1s;
D2  = D2s;
% Plus boundary
for jz = 1:r
    for jm = 1:m
        D1(jz,jm) = (a_b1(jz,jm))./dz;
        D2(jz,jm) = (a_b2(jz,jm))./dz^2;
    end
end

% matrices with pole condition
% Pole conditions
for j = 1:w
    D1s(j,1:w+j) = a_p1s(j,w+2-j:2*w+1)./dz;
    D1a(j,1:w+j) = a_p1a(j,w+2-j:2*w+1)./dz;
    D2s(j,1:w+j) = a_p2s(j,w+2-j:2*w+1)./dz^2;
    D2a(j,1:w+j) = a_p2a(j,w+2-j:2*w+1)./dz^2;
end

% Make sparse
D1  = sparse(D1);
D2  = sparse(D2);
D1s = sparse(D1s);
D1a = sparse(D1a);
D2s = sparse(D2s);
D2a = sparse(D2a);


end

function a_mod = Mohseni_pole(a_interior,h)


% Half-width of finite difference coefficients 
w = (length(a_interior(1,:))-1)/2;

% Storage for modified FD coefficients
a_mod = zeros(w,2*w+1);

% Find modified coefficeints near pole
for j = 0:w-1
    
    % Unaffected nodes
    a_mod(j+1,2*w+1-2*j:2*w+1) = a_interior(1,2*w+1-2*j:2*w+1);
    
    % Affected nodes: out of bounds nodes are reflected over axis
    a_mod(j+1,w+1-j:2*w-2*j) = a_interior(1,w+1-j:2*w-2*j) + h*fliplr(a_interior(1,1:w-j));

end




end