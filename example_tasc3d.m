% recreate Fig. 2 and 3 of the paper
% set cas=21 for 2a, cas=22 for 2b, 23 for 2c and 24 for 2d
% set cas=31 for 3a, 32 for 3b, 33 for 3c and 34 for 3d
%
% Reference: Marcotte D., 2015. TASC3D : A program to test the
% admissibility in 3D of non-linear models of coregionalization. Computers
% & Geosciences
fstitle=16;
fslabel=14;
fsaxes=12;
fslegend=12;
s=[];      % this uses the set of default frequencies
global wait

switch cas
    
    case 21 % Fig. 2a
        dy=0.00;
        model=cell(2,2);c=cell(2,2);
        model{1,1}=[7 100 1 1];
        model{1,2}=[7 nan 1 0.9];
        model{2,1}=model{1,2};
        model{2,2}=[7 100 1 1];
        vpar=cell(2,1);
        vpar{1}=[60:0.1:150];
        vpar{2}=[];
        
        [icode9,ics9]=tasc3d(model,vpar,s)
        
        model{1,2}=[7 nan 1 0.7];
        [icode7,ics7]=tasc3d(model,vpar,s)
        
        model{1,2}=[7 nan 1 0.5];
        [icode5,ics5]=tasc3d(model,vpar,s)
        
        vp=vpar{1};
        
        plot(vp,ics9,'-k',vp,ics7+dy,'--k',vp,ics5+2*dy,'-.k',vp,icode9,'-r',vp,icode7+dy,'--r',vp,icode5+2*dy,'-.r','linewidth',2);
        axis([60 150 -0.1 2])
        t=legend('CS, R=0.9','CS, R=0.7','CS, R=0.5','Admiss, R=0.9','Admiss, R=0.7','Admiss, R=0.5')
        
        set(t,'fontsize',fslabel)
        title('a)','fontsize',fstitle)
        xlabel('Range','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        
        
    case 22 % Fig. 2b)
        model=cell(2,2);c=cell(2,2);
        model{1,1}=[5 100 1 1];
        model{1,2}=[5 nan 1 0.9];
        model{2,1}=model{1,2};
        model{2,2}=[5 100 1 1];
        vpar=cell(2,1);
        vpar{1}=[60:0.1:150];
        vpar{2}=[];
       
        s=[1:2:11 11.5:.01:11.6]'; % we need a fine discretisation arond the zero of the cubic at ~11.527
        
        [icode9,ics9]=tasc3d(model,vpar,s)
        
        model{1,2}=[5 nan 1 0.7];
        [icode7,ics7]=tasc3d(model,vpar,s)
        
        model{1,2}=[5 nan 1 0.5];
        [icode5,ics5]=tasc3d(model,vpar,s)
        
        
        vp=vpar{1};
        plot(vp,ics9,'-k',vp,ics7,'--k',vp,ics5,'-.k',vp,icode9,'-r',vp,icode7,'--r',vp,icode5,'-.r','linewidth',2); 
        
        axis([60 150 -0.1 2])
        t=legend('CS, R=0.9','CS, R=0.7','CS, R=0.5','Admiss, R=0.9','Admiss, R=0.7','Admiss, R=0.5')
        set(t,'fontsize',fslegend)
        title('b)','fontsize',fstitle)
        xlabel('Range','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        
    case 23 % Fig 2c)
        model=cell(2,2);c=cell(2,2);
        model{1,1}=[8 100 1.0 1];
        model{1,2}=[8 nan nan 0.7];
        model{2,1}=model{1,2};
        model{2,2}=[8 100 2.0 1];
        vpar=cell(2,1);
        vpar{1}=[30:0.2:150];
        vpar{2}=[0.5:0.01:4];
        [icode,ics]=tasc3d(model,vpar,s)
        
        
        title('c)','fontsize',fstitle)
        xlabel('Range','fontsize',fslabel)
        ylabel('Shape','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        
    case 24 % Fig. 2d
        
        model=cell(2,2);c=cell(2,2);
        model{1,1}=[2 100 0 1];
        model{1,2}=[2 nan 0 nan];
        model{2,1}=model{1,2};
        model{2,2}=[2 100 0 1];
        vpar=cell(2,1);
        vpar{1}=[85:0.05:120];
        vpar{2}=[0.2:0.005:1.1];
        [icode,ics]=tasc3d(model,vpar,s)
        
        title('d)','fontsize',fstitle)
        xlabel('Range','fontsize',fslabel)
        ylabel('Sill','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        
        
    case 31  % Fig. 3a
        
        model=cell(2,2);
        model{1,1}=[2 58 0 0.8;3 40 0 0.2];
        model{2,2}=[3 22 0 0.13;3 53 0 0.65];
        model{1,2}=[3 nan 0 0.06;3 nan 0 0.4];
        model{2,1}=model{1,2};
        
        vpar=cell(2,1);
        
        vpar{1}=[5:0.1:10.5 11.5:0.5:30];
        vpar{2}=[30:0.1:80];
        
        [icodesp,icssp]=tasc3d(model,vpar,s)
        
        title('a)','fontsize',fstitle)
        xlabel('Range comp 1','fontsize',fslabel)
        ylabel('Range comp 2','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        
        caxis([0 2])
        
    case 32     % Fig. 3b)
        vpar=cell(2,1);%
        
        model=cell(2,2);c=cell(2,2);
        model{1,1}=[2 60 25 25 0 0 nan 1 1];
        model{1,2}=[3 100 25 25 0 0 50 1 nan];
        model{2,1}=model{1,2};
        model{2,2}=[3 100 25 25 0 0 50 1 1];
        
        
        vpar{1}=[20:0.5:80];
        vpar{2}=[0.1:0.01:0.98];
        [icode9,ics9]=tasc3d(model,vpar,s)
        
        title('b)','fontsize',fstitle)
        xlabel('Rot_z','fontsize',fslabel)
        ylabel('Correlation','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        caxis([0 2])
        
    case 33   % Fig. 3c)
        
        model=cell(2,2);
        model{1,1}=[1 1 1 0.2;6 60 1 1];
        model{1,2}=[6 60 1 nan];
        model{2,2}=[1 1 1 0.5;6 60 1 1];
        
        vpar=cell(2,1);
        vpar{1}=[0.9:0.001:1 1.000001 1.05];
        vpar{2}=[];
        [icode9,ics9]=tasc3d(model,vpar,s)
        axis([0.9 1.05 -0.1 2])
        grid off
        title('c)','fontsize',fstitle)
        xlabel('Correlation','fontsize',fslabel)
        set(gca,'fontsize',fsaxes)
        t=legend('CS','Admiss')
        set(t,'fontsize',fslegend)
        set(gca,'ytick',[0 1])
        
    case 34  % Fig. 3d)
        
        model=cell(3,3);
        model{1,1}=[4 60 1 1];
        model{1,2}=[4 60 1 nan];
        model{1,3}=[4 60 1 nan];
        model{2,2}=[4 60 1 1];
        model{2,3}=[4 60 1 0.9];
        model{3,3}=[4 60 1 1];
        
        vpar=cell(2,1);
        vpar{1}=[0.6:0.002:0.99 0.991:0.001:1.001 1.1];
        vpar{2}=[0.6:0.002:0.99 0.991:0.001:1.001 1.1];
        [icode9,ics9]=tasc3d(model,vpar,s)
        
        title('d)','fontsize',fstitle)
        xlabel('Correlation 1-2','fontsize',fslegend)
        ylabel('Correlation 1-3','fontsize',fslegend)
        axis equal;set(gca,'fontsize',fsaxes)
        axis([0.6 1.1 0.6 1.1]);
end



