

%simple system
% Xin ---> S ---> Xout
%      e1     e2
% where Xin and Xout are external parameters, e1 e2 are enzyme abundances
% given by

% e1 = 1 + b * c1
% with c1 the distance of the gene from the terminus
% in this way if c=0 (a gene at the terminus) then e=1

Xin=2;
k1=0.25;
k2=0.25;
c1=2000000;
c2=50000;
mx=2000000;
b=2*10^-06;
x=[c1,c2];

%steady state solution
%Sss= (1+b*x(1))/(1+b*x(2)) * k1/k2 * Xin;

%sensitivity to beta
%dS_db = @(x)(k1 / (k2*Sss) * Xin * (x(1)-x(2))*b / (1+b*x(2)).^2) ;

%define exponents for getting beta (beta=10^exponent)
ball=logspace(-7,-4,20);

nrep=10000;

%declare vector of zeros to store values
coefvarS=zeros(nrep,1);
%define gene locations
tmp=10.^linspace(0,6,nrep);

%use the same coordinates for both genes after shuffling
c2s=tmp(randperm(length(tmp)));
c1s=tmp(randperm(length(tmp)));

%define steady state with beta=1 and equal abundance of the enzymes
Sss0=  k1/k2 * Xin;


	for n=1:nrep
		%get coordinate for gene 1 and 2
		x=[c1s(n),c2s(n)];
		tmp=[];
			%now calculate steady state metabolite concentration
			%when beta is changed
			for i=1:length(ball)
    			Sss= (1 + ball(i) * x(1)) / (1 + ball(i) * x(2)) * k1 / k2 * Xin;
    			%accumulate all values
    			tmp(i)=(Sss) ;
			end
			
		coefvarS(n)=std(tmp)/mean(tmp);
	end



CVSS=(coefvarS)
%stuff for surface plotting
xs=size(colormap(parula),1);
range=max((CVSS))-min((CVSS));
step=range/xs;
idx=round((1-min((CVSS)))/step);
CM=colormap(parula);
%CM((idx-10):(idx+9),:)=repmat([1,0,0],20,1);
figure(1)
%subplot(1,2,1)
xv=linspace(min(log10(c1s)),max(log10(c1s)),500);
yv=linspace(min(log10(c2s)),max(log10(c2s)),500);
[X,Y]=meshgrid(xv,yv);
Z=griddata(log10(c1s),log10(c2s),CVSS,X,Y);

surf(X,Y,Z)
xlabel('log10(p_2)');
ylabel('log10(p_1)');
zlabel('std Sss / <Sss>');
colormap(CM)
%set(gca,'xlim',[0,6])
%set(gca,'ylim',[0,6])
shading interp


