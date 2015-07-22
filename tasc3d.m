function [icode,ics]=tasc3d(model,vpar,s)
%%
% Function to test the admissibility of a non-lmc model.
% syntax: [icode,ics]=tasc3d(model,vpar,s)
%
% model: nvar x nvar cell array with model{i,j}=model{j,i};
%        Each cell as n_ij structures
%        Each structure is specified by 4 (isotropic case) or 9
%        (anisotropic case) parameters
%            iso: model type, range, shape, sill   (shape is alpha for the
%            Cauchy and nu for the K-Bessel)
%            aniso: model type, range 1, range 2, range 3, rot x, rot y rot z, shape, sill
%                   the rotations are all done conterclockwise, z first,
%                   then rotated y, then rotated x.
% 
% One or two parameters can be specifed as nan. In this case, the program finds the domain of admissible values.
% The program gives the Cauchy-Schawrtz admissible domain (necessary
% condition) and the real admissible domain.
% more than two nan parameters: a message is issued and the execution stops
%
% vpar: a cell(2 x 1). Each entry is a vector giving the values to consider for
% the nan paramater (in order where they appear).
%
% s: a vector with frequencies where to test admissibility (in the anisotropic case, this
% vector should be rather short as it is used to define a 3D grid; if the vector is of size 100, then the grid of frequencies is 1 million nodes)
% typical example: s=[0.005:0.005:0.01 0.02:0.02:0.1 0.2:0.2:7]';
%
% Output: icode= 1 the model is admissible, 0 it is not.
%         ics=1, the model respect the C-S condition, 0, it does not. One
%         always has icode<= ics
%         icode and ics are vectors the same size as vpar{1} if a single
%         nan is provided, they are matrices when 2 nans are found
%
% Author: D. Marcotte April 2015

nvar=size(model,1);
an=-1;cn=-1;
nt=0;
vni=[];  % holds the position of free parameters
ncol=size(model{1,1},2);  % number of parameters in anisotropic model
iso=ncol<9;  % isotropic case

if isempty(s)
   s=10.^[-3:0.1:2]';  % set of default frequencies where to test admissibility
end

if iso
    ir=2;
    inu=3;
    is=4;
else
    ir=4;
    inu=8;
    is=9;
end

for i=1:nvar;
    for j=i:nvar
        an=max(an,nanmax(nanmax(model{i,j}(:,2:ir))));  % fing the largest range
        cn=max(cn,nanmax(model{i,j}(:,end)));
    end
end

icode=[];ics=[];

% some normalization
for i=1:nvar;
    for j=i:nvar
        model{i,j}(:,2:ir)= model{i,j}(:,2:ir)/an;
        model{i,j}(:,is)=model{i,j}(:,is)/cn;
        model{j,i}=model{i,j};
    end
end

icode=check_nugget(model);  % separate treatment for the nugget effect
if icode==0                 % stop if nugget is not admissible
   ics=nan;
    return
end

for i=1:nvar;
    for j=1:nvar
        id=model{i,j}(:,1)>1; % identify components different from nugget
        model{i,j}=model{i,j}(id,:); % keep only those components for the rest
        model{j,i}=model{i,j};
    end
end

for i=1:nvar;
    for j=i:nvar
        ni=find(isnan(model{i,j}(:)));
        ntt=length(ni);
        [subi,subj]=ind2sub(size(model{i,j}),ni);
        vni=[vni;i*ones(ntt,1) j*ones(ntt,1) subi(:) subj(:)];
        nt=nt+ntt;
    end
end

for i=1:nt
    if vni(i,4)>1 && vni(i,4)<=ir
        vpar{i}=vpar{i}/an;
    elseif vni(i,4)==is
        vpar{i}=vpar{i}/cn;
    end
end

if nt==0  % there is no free parameters
    ics=test_cs(model,iso);
    if ics==1
        icode=check_admiss(model,iso,s);
    else
        icode=0;
    end
elseif nt<=2
    ics=draw_cs(model,vni,vpar,iso);              % compute the CS domain for the free parameters
    icode=find_admiss(model,vni,vpar,ics,iso,s);  % compute the admissibility domain
    showgraph(ics,icode,model,vni,vpar,an,cn,iso,ir,is); % draw figure if one or two free parameters
else
    disp('Number of free parameters >2; reduce this number')
    return
end

function admiss=draw_cs(model,vni,vpar,iso)
%% This function calls test_cs at all parameter values specified in vpar
global wait
nvar=size(model,1);
wait=waitbar(0,'Check of the CS condition');
n1=length(vpar{1});

if size(vni,1)==1;
    vpar=vpar{1};
    for i=1:n1
        waitbar(i/(2*n1),wait);
        model{vni(1,1),vni(1,2)}(vni(1,3),vni(1,4))=vpar(i);
        admiss(i)=test_cs(model,iso);
    end
    
elseif size(vni,1)==2
    for i=1:n1
        waitbar(i/(2*n1),wait);
        model{vni(1,1),vni(1,2)}(vni(1,3),vni(1,4))=vpar{1}(i);
        for j=1:length(vpar{2})
            model{vni(2,1),vni(2,2)}(vni(2,3),vni(2,4))=vpar{2}(j);
            admiss(i,j)=test_cs(model,iso);
        end
    end
    
end


function ics=test_cs(model,iso);
%% This function checks the Cauchy-Schwarz condition
h=[0 0.001:0.001:0.01 0.02:0.01:0.1 0.2:0.1:1, 2:1:10 15:5:50]'; % vector of normalized distances where CS is checked
if ~iso
  [h1,h2,h3]=ndgrid(h,h,h);
  h=[h1(:) h2(:), h3(:)];
end
nvar=size(model,1);
ics=1; % the model meets CS condition

for i=1:nvar-1;
    g11=sum(model{i,i}(:,end))-covar(h,zeros(1,size(h,2)),model{i,i}(:,1:end-2),model{i,i}(:,end),model{i,i}(:,end-1));
    for j=i+1:nvar
        g22=sum(model{j,j}(:,end))-covar(h,zeros(1,size(h,2)),model{j,j}(:,1:end-2),model{j,j}(:,end),model{j,j}(:,end-1));
        g12=sum(model{i,j}(:,end))-covar(h,zeros(1,size(h,2)),model{i,j}(:,1:end-2),model{i,j}(:,end),model{i,j}(:,end-1));
        d=sqrt(g11.*g22)-abs(g12);
        if min(d)<-1e-09
            ics=0;  % one of the models does not meet CS condition
            return
        end
    end
end


function showgraph(ics,icode,model,vni,vpar,an,cn,iso,ir,is);
%% This function draws the figure showing CS condition and admissibility
figure(100);

nvar=size(model,1);
vpar1=vpar{1};
if length(vpar)>1;
    vpar2=vpar{2};
end

if iso
    tit=char('Type','Range','Shape','Sill');
else
    tit=char('Type','Range 1','Range 2','Range 3','Rot x','Rot y','Rot z','Shape','Sill');
end

if size(vni,1)==1
    if vni(1,4)>1 && vni(1,4)<=ir     % range
        vpar1=vpar1*an;
    elseif vni(1,4)==is              % sill
        vpar1=vpar1*cn;
    end
    
    plot(vpar1,ics,'-k',vpar1,icode,'--r','linewidth',2)
    title(['Model (',num2str(vni(1,1)),',',num2str(vni(1,2)),' ; ',num2str(vni(1,3)),') ', tit(vni(1,4),:)],'fontsize',14)
    grid on
    ylim([0 2]);
    t=legend('CS','Admiss');
    set(t,'fontsize',12);
    set(gca,'ytick',1)
else
    if vni(1,4)>1 && vni(1,4)<=ir     % range
        vpar1=vpar1*an;
    elseif vni(1,4)==is              % sill
        vpar1=vpar1*cn;
    end
     if vni(2,4)>1 && vni(2,4)<=ir    % range
        vpar2=vpar2*an;
    elseif vni(2,4)==is              % sill
        vpar2=vpar2*cn;
    end
    
    [X,Y]=ndgrid(vpar1,vpar2);
    contourf(X,Y,ics+icode,[-1:2]);
    xlabel(['Model (',num2str(vni(1,1)),',',num2str(vni(1,2)),' ; ',num2str(vni(1,3)),') ', tit(vni(1,4),:)],'fontsize',14)
    ylabel(['Model (',num2str(vni(2,1)),',',num2str(vni(2,2)),' ; ',num2str(vni(2,3)),') ', tit(vni(2,4),:)],'fontsize',14)
    colormap gray
    title('2: admissible, 1: Cauchy-Schwartz ok but not admissible')
    colorbar
end

function icode=find_admiss(model,vni,vpar,ics,iso,s);
%% Caller for check_admiss at all values shown in vpar
global wait
[n,m]=size(ics);
icode=zeros(n,m);

nvar=size(model,1);
n1=length(vpar{1});
waitbar(0.5,wait,'Checking the admissibility')
if size(vni,1)==1;
    for i=1:n1
        waitbar(0.5+i/(2*n1),wait);
        if ics(i)==1
            model{vni(1),vni(2)}(vni(3),vni(4))=vpar{1}(i);
            model{vni(2),vni(1)}(vni(3),vni(4))=vpar{1}(i);
            icode(i)=check_admiss(model,iso,s);
        end
    end
    
    
elseif size(vni,1)==2
    for i=1:n1
        waitbar(0.5+i/(2*n1),wait);
        model{vni(1,1),vni(1,2)}(vni(1,3),vni(1,4))=vpar{1}(i);
        model{vni(1,2),vni(1,1)}(vni(1,3),vni(1,4))=vpar{1}(i);
        for j=1:length(vpar{2})
            if ics(i,j)==1
                model{vni(2,1),vni(2,2)}(vni(2,3),vni(2,4))=vpar{2}(j);
                model{vni(2,2),vni(2,1)}(vni(2,3),vni(2,4))=vpar{2}(j);
                icode(i,j)=check_admiss(model,iso,s);
            end
        end
    end
    
end
delete(wait)



function icode=check_admiss(model,iso,s);
%% This function checks the admissibility of a peculiar set of parameters
nvar=size(model,1);

icode=1;

% write down spectral densities of the models

nugget='0';  % place holder, not used
expon='at/pi^2./(1+sa.^2).^2';
gaus='at/(8*pi^1.5)*exp(-0.25*sa.^2)';
spher='0.75*at/pi./sa.^3.*besselj(1.5,sa/2).^2';
cubic='210*at./(pi^2*sa.^10).*(6*sa.*cos(sa/2)+(sa.^2-12).*sin(sa/2)).^2';
penta='27720*at./(pi^2*sa.^14).*((sa.^3-60*sa).*cos(sa/2)+(120-12*sa.^2).*sin(sa/2)).^2';
cauchy='at*sa.^(nu-1.5)/(pi^1.5*2^(nu+0.5)*gamma(nu)).*besselk(1.5-nu,sa)';
bess='at*gamma(nu+3/2)/gamma(nu)/pi^1.5./((1+sa.^2).^(nu+1.5))';
Spec=char(nugget,expon,gaus,spher,cubic,penta,cauchy,bess);

if ~iso
   ss2=[-s(end:-1:1) s]; % consider the negative side as well
   [s1,s2,s3]=ndgrid(s,ss2,ss2);
   s=[s1(:),s2(:),s3(:)];
   clear s1 s2 s3
end

fs=zeros(length(s),nvar,nvar);

for i=1:nvar
    for j=i:nvar
        for k=1:size(model{i,j},1);
            if iso
                at=model{i,j}(k,2).^3;
                sa=s*model{i,j}(k,2);
            else
                [~,rot]=trans([0 0 0],model{i,j}(:,1:7),k);
                at=model{i,j}(k,2)*model{i,j}(k,3)*model{i,j}(k,4);
                sa=s*rot;
                sa=sqrt(sa.^2*((model{i,j}(k,2:4)').^2));
            end
            nu=model{i,j}(k,end-1);
            if model{i,j}(k,1)==6   % special treatment for the unstable penta
                sa=max(sa,0.1);     % below 0.1 the model is numerically unstable
            end
            if model{i,j}(k,1)==5   % special treatment for the unstable penta
                sa=max(sa,0.03);     % below 0.04 the model is numerically unstable
            end
            f=eval(Spec(model{i,j}(k,1),:));
            fs(:,i,j)=fs(:,i,j)+f*model{i,j}(k,end);
        end
        fs(:,j,i)=fs(:,i,j);    
    end
end
fs=permute(fs,[2,3,1]);
for i=1:size(s,1)
    ei=min(eig(fs(:,:,i)));
    if ei<-1e-09;  % below -10^-9, it is considered negative, hence not admissible
        icode=0;
        [i,ei,s(i)]; % remove the ; to show at which frequency negative eigenvalue first occurs
        return
    end
end

function [k]=covar(x,x0,model,c,vnu)
%% function to compute the covariances 

% we define the equations for the various covariograms. Any new model
% can be added here after bess (and its spectral density added in
% check_admiss also after bess)

k=[];
nugget='h==0';
expon='exp(-h)';
gaus='exp(-(h).^2)';
spher='1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)';
cubic='1-(7*min(h,1).^2-8.75*min(h,1).^3+3.5*min(h,1).^5-0.75*min(h,1).^7)';
penta='1-22/3*min(h,1).^2+33*min(h,1).^4-77/2*min(h,1).^5+33/2*min(h,1).^7-11/2*min(h,1).^9+5/6*min(h,1).^11';
cauchy='(h.^2+1).^(-nu)';
bess='1/(gamma(nu)*2^(nu-1))*max(h,eps).^nu.*besselk(nu,max(h,eps))';

Gam=char(nugget,expon,gaus,spher,cubic,penta,cauchy,bess);

% definition of some constants
[n1,d]=size(x); 
[n2,d]=size(x0);
[rp,p]=size(c);
r=rp/p;  % number of structures
cx=[x(:,1:d);x0];
nm=size(model,2);

% avoid division by zeros
if nm>2
    model(:,2:1+d)=max(model(:,2:1+d),100*eps);
else
    model(:,2)=max(model(:,2),100*eps);
end

k=zeros(n1*p,n2*p);
for i=1:r,
    
    % calculation of matrix of reduced rotated distances H
    
    [t1]=trans(x(:,1:d),model,i);
    [t2]=trans(x0,model,i);
    h=0;
    for id=1:d
        h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
    end
    h=sqrt(h);
    ji=(i-1)*p+1; js=i*p ;
    nu=vnu(i);
    % evaluation of the current basic structure
    g=eval(Gam(model(i,1),:));
    
    k=k+kron(g,c(ji:js,:));
end

function icode=check_nugget(model);
%% function to check that the matrix of nugget coefficients is positive
% definite (necessary condition)
nvar=size(model,1);
nug=zeros(nvar,nvar);
for i=1:nvar
    for j=i:nvar
        nc=size(model{i,j},1);
        for k=1:nc
            if model{i,j}(k,1)==1;
                nug(i,j)=model{i,j}(k,end);
                nug(j,i)=nug(i,j);
            end
        end
    end
end
ei=min(eig(nug));
icode=ei>= 0;

