
%%%%%%%%%%%%%%
%%% Explanation of the chosen notation
%%%
%%% The variable x containts 
%%% x = [R_2, ..., R_(N-1),\Phi^(B1g)_2, ... , \Phi^(B1g)_(N-1), \Re\phi^(A1g)_2, ..., \Re\phi^(A1g)_(N-1), \Re\phi^(B1g)_2, ..., \Re\phi^(B1g)_(N-1)].
%%%
%%% Note that below we define C \equiv \Phi, and c \equiv \phi.
%%% Also note that for efficiency, we omitted all fields which stay zero throughout the simulation. If need be, they can easily be attached.
%%% We also chose the superconducting order to occur on the real axis, such that \phi^(A1g) and \phi^(B1g) are also real.
%%%%%%%%%%%%%%




close all;clear all;

%Declaring global variables
global edgeR edgeCB1 s0 sx sy sz dz u v w t0h r0 eps Nb N tau0 tauz taux tauy sp sm MA1g MA2g MB1g MB2g mA1g mB1g mB2g FA1g count

%%%%%%%%%%
%%% defining Pauli matrices, and M^n, m^n matrices
%%%%%%%%%%
tau0=[1 0;0 1];
taux=[0 1;1 0];
tauy=[0 -1i;1i 0];
tauz=[1 0;0 -1];
s0=[1 0;0 1];
sx=[0 1;1 0];
sy=[0 -1i;1i 0];
sz=[1 0;0 -1];
sp=1/2*(sx+1i*sy);
sm=1/2*(sx-1i*sy);
MA1g=1/2*kron(s0,tau0);
MB1g=1/2*kron(s0,tauz);
MB2g=1/2*kron(s0,taux);
MA2g=1/2*kron(sz,tauy);
mA1g=kron(sm,tau0);
mB1g=kron(sm,tauz);
mB2g=kron(sm,taux);





N=20;   %Number of lattice points
Nb=4;   %Dimension of Nambu subspace; by default Nb=4


eps=10^(-8);




u=1;
FA1g=1*MA1g;
t0=20/7;
t0h=t0/sqrt(u);

%%%%%%%%%%
%%% Parameter for local charge-4e condensation
%%%%%%%%%%
%v=2.1*u;w=-1*u/20;
%r0=-5.8/t0;edgeR=1.7862/t0;edgeCB1=1.3264/t0;
%r0=-6.3/t0;edgeR=1.7862/t0;edgeCB1=1.3264/t0;
%r0=-6.8;edgeR=1.967;edgeCB1=1.5655;d0=5/7;

%%%%%%%%%%
%%% Parameter for pure nematic phase
%%%%%%%%%%
v=0.1*u;w=-1/2*u;
r0=-4.6/t0;edgeR=2.7862/t0;edgeCB1=2.3264/t0;




%%%%%%%%%%
%%% Variable constraints fed to the minimizer.
%%% (1) lb \leq x \leq ub
%%% (2) A * x \leq b
%%%%%%%%%%
lb=zeros(1,4*(N-2));
for i=1:(N-2)
    lb((N-2)+i)=-1.5*edgeCB1;
    lb(2*(N-2)+i)=-1.5*edgeCB1;
    lb(3*(N-2)+i)=-1.5*edgeCB1;    
end
ub=2*edgeR*ones(1,4*(N-2));
A=zeros(4*(N-2),4*(N-2));
for i=1:(N-2)
    A(i,i)=-1;
    A(i,(N-2)+i)=-1;
    A(i,2*(N-2)+i)=-1;
    A(i,3*(N-2)+i)=-1; 
    A((N-2)+i,i)=-1;
    A((N-2)+i,(N-2)+i)=-1;
    A((N-2)+i,2*(N-2)+i)=1;
    A((N-2)+i,3*(N-2)+i)=1; 
    A(2*(N-2)+i,i)=-1;
    A(2*(N-2)+i,(N-2)+i)=1;
    A(2*(N-2)+i,2*(N-2)+i)=-1;
    A(2*(N-2)+i,3*(N-2)+i)=1; 
    A(3*(N-2)+i,i)=-1;
    A(3*(N-2)+i,(N-2)+i)=1;
    A(3*(N-2)+i,2*(N-2)+i)=1;
    A(3*(N-2)+i,3*(N-2)+i)=-1;    
end
b=zeros(4*(N-2),1);

%%%%%%%%%%
%%% Add additional constraint to enforce the domain wall center 
%%% to be located in the center of the chain.
%%% By default, this is on.
%%%%%%%%%%
forceCB1zeroCenter=1;
if forceCB1zeroCenter==1  %force the domain wall to be centered in the chain
   initial((N-2)+N/2)=0;
   ub((N-2)+N/2)=10^(-6);
   lb((N-2)+N/2)=-10^(-6);
   ub((N-2)+N/2+1:2*(N-2))=0;
   lb((N-2)+1:(N-2)+N/2-1)=0;     
end

%%%%%%%%%%
%%% Initial x - values
%%% The result should not much depend on it.
%%% But to improve efficieny, we choose a profile not too far away from the result.
%%%%%%%%%%
initial = zeros(1,4*(N-2));    % Number of variables
for i=1:(N-2)
    initial(i)=0.96*edgeR+0.0*(rand-0.5);
    initial((N-2)+i)=-0.96*edgeCB1*tanh(i-(N-2)/2)+0.0*(rand-0.5);
    initial(2*(N-2)+i)=0.0*(rand-0.5);
    initial(3*(N-2)+i)=0.0*(rand-0.5);   
end

%%%%%%%%%%
%%% The running minimizer routine.
%%%
%%% Note that we also supply the gradient function for a significantly improved performance.
%%% During the run, the progress is plotted after every step, i.e. x-values plus gradient-x -values.
%%% The gradient values should eventually approach zero.
%%%%%%%%%%
options=optimoptions(@fmincon,'MaxFunctionEvaluations',4000,'SpecifyObjectiveGradient',true,'StepTolerance',10^(-10));
figure; 
nonlcond=@nonlfunc;
[x,fval]=fmincon(@free_energy,initial,A,b,[],[],lb,ub,[],options);

%%%%%%%%%%
%%% After the run, we plot the final configuration,
%%% and save it to the matrix file Xmat.mat
%%%%%%%%%%
count=0;
dispFin=['F=',num2str(fval),' R=',num2str(x(1)),' CB1=',num2str(x((N-2)+1)),...
    ' cA1=',num2str(x(2*(N-2)+1)),' cB1=',num2str(x(3*(N-2)+1))];
disp(dispFin);
%save('Xmat.mat','x');


%%%%%%%%%%
%%% Last variable constraints which is almost similar to the above constraint
%%% (2) A * x \leq b,
%%% just a bit more restrictive.
%%% In particular, we enforce that the system cannot end up in the superconducting regime where the equations are not valid.
%%% I.e. we enforce the inverse superconducting susceptibility \chi >0 to be positive at every lattice site.
%%% Nonlinear constraint reads c(x) \leq 0. 
%%%%%%%%%%
function [c,ceq] = nonlfunc(x)
global N
    c=zeros(2*(N-2),1);
    for i=1:(N-2)
        c(i,1)=-x(i)+x(1*(N-2)+i)+abs(x(2*(N-2)+i)-x(3*(N-2)+i));
        c((N-2)+i,1)=-x(i)-x(1*(N-2)+i)+abs(x(2*(N-2)+i)+x(3*(N-2)+i));
    end
ceq = [];
end

%%%%%%%%%%
%%%%%%%%%%
%%% FUNCTIONS
%%%%%%%%%%
%%%%%%%%%%

%%%%%%%%%%
%%% Local inverse Greens function
%%%%%%%%%%
function Gm1=Gm1_func(R,CB1,CB2,CA2,cA1,cB1,cB2)
    global MA1g MA2g MB1g MB2g mA1g mB1g mB2g
    Gm1=2*R*MA1g+2*CB1*MB1g+2*CB2*MB2g+2*CA2*MA2g+cA1*mA1g'+conj(cA1)*mA1g+cB1*mB1g'+conj(cB1)*mB1g+cB2*mB2g'+conj(cB2)*mB2g;
end
%%%%%%%%%%
%%% Total inverse Greens function
%%%%%%%%%%
function mcalGm1=mcalGm1_func(x)
    global N Nb FA1g edgeR edgeCB1
    mcalGm1=zeros(Nb*N,Nb*N);
    mcalGm1(1:Nb,1:Nb)=Gm1_func(edgeR,edgeCB1,0,0,0,0,0)+FA1g;
    for r=2:(N-1)
        mcalGm1((r-1)*Nb+1:r*Nb,(r-1)*Nb+1:r*Nb)=...
            Gm1_func(x((r-1)),x((N-2)+(r-1)),0,0,x(2*(N-2)+(r-1))+1i*0,x(3*(N-2)+(r-1))+1i*0,0)+FA1g;
        mcalGm1((r-1)*Nb+1:(r)*Nb,(r-2)*Nb+1:(r-1)*Nb)=-0.5*FA1g;
        mcalGm1((r-2)*Nb+1:(r-1)*Nb,(r-1)*Nb+1:(r)*Nb)=-0.5*FA1g;
    end
        mcalGm1((N-1)*Nb+1:N*Nb,(N-1)*Nb+1:N*Nb)=Gm1_func(edgeR,-edgeCB1,0,0,0,0,0)+FA1g;
        mcalGm1((N-1)*Nb+1:(N)*Nb,(N-2)*Nb+1:(N-1)*Nb)=-0.5*FA1g;
        mcalGm1((N-2)*Nb+1:(N-1)*Nb,(N-1)*Nb+1:(N)*Nb)=-0.5*FA1g;   
 end
%%%%%%%%%%
%%% Free energy contribution of a given lattice site, see SM of paper.
%%% Will be summed over inside free_energy(x)
%%%%%%%%%%
function Fadd_r=Fadd_r_func(R,CB1,CB2,CA2,cA1,cB1,cB2,G0A1,G0B1,G0B2,GzA2,GpA1,GpB1,GpB2)
    global r0 u v w t0h
    Fadd_r=2*((r0-R)*G0A1-CB1*G0B1-CB2*G0B2-CA2*GzA2-1/2*(cA1*conj(GpA1)+conj(cA1)*GpA1)-1/2*(cB1*conj(GpB1)+conj(cB1)*GpB1)-1/2*(cB2*conj(GpB2)+conj(cB2)*GpB2))...
        +2/t0h^2*((3*u+v+w)*G0A1^2+(u-v+3*w)*G0B1^2+(u-v-w)*G0B2^2+(u+3*v-w)*GzA2^2+(u-v+w)*abs(GpA1)^2+(u+v+w)*abs(GpB1)^2+(u+v-w)*abs(GpB2)^2);
end


%%%%%%%%%%
%%% Gradient contributions coming from a given lattice site, see SM of paper.
%%% Will be summed over inside free_energy(x)
%%%%%%%%%%
function gradR_add=gradR_add_func(R,CB1,CB2,CA2,cA1,cB1,cB2,Grprp,Grpr,Grrp)
    global r0 u v w t0h MA1g MB1g MB2g MA2g mA1g mB1g mB2g
    gradR_add=2*(r0-R+2/t0h^2*(3*u+v+w)*1/2*trace(MA1g*Grprp))*(-trace(MA1g*Grpr*MA1g*Grrp))...
        +2*(-CB1+2/t0h^2*(u-v+3*w)*1/2*trace(MB1g*Grprp))*(-trace(MB1g*Grpr*MA1g*Grrp))...
        +2*real(2*(-1/2*conj(cA1)+1/t0h^2*(u-v+w)*1/2*trace(mA1g'*Grprp))*(-trace(mA1g*Grpr*MA1g*Grrp)))...
        +2*real(2*(-1/2*conj(cB1)+1/t0h^2*(u+v+w)*1/2*trace(mB1g'*Grprp))*(-trace(mB1g*Grpr*MA1g*Grrp)));
end
function gradCB1_add=gradCB1_add_func(R,CB1,CB2,CA2,cA1,cB1,cB2,Grprp,Grpr,Grrp)
    global r0 u v w t0h MA1g MB1g MB2g MA2g mA1g mB1g mB2g
    gradCB1_add=2*(r0-R+2/t0h^2*(3*u+v+w)*1/2*trace(MA1g*Grprp))*(-trace(MA1g*Grpr*MB1g*Grrp))...
        +2*(-CB1+2/t0h^2*(u-v+3*w)*1/2*trace(MB1g*Grprp))*(-trace(MB1g*Grpr*MB1g*Grrp))...
        +2*real(2*(-1/2*conj(cA1)+1/t0h^2*(u-v+w)*1/2*trace(mA1g'*Grprp))*(-trace(mA1g*Grpr*MB1g*Grrp)))...
        +2*real(2*(-1/2*conj(cB1)+1/t0h^2*(u+v+w)*1/2*trace(mB1g'*Grprp))*(-trace(mB1g*Grpr*MB1g*Grrp)));
end
function gradRecA1_add=gradRecA1_add_func(R,CB1,CB2,CA2,cA1,cB1,cB2,Grprp,Grpr,Grrp)
    global r0 u v w t0h MA1g MB1g MB2g MA2g mA1g mB1g mB2g
    gradRecA1_add=2*(r0-R+2/t0h^2*(3*u+v+w)*1/2*trace(MA1g*Grprp))*(-1/2*trace(MA1g*Grpr*(mA1g'+mA1g)*Grrp))...
        +2*(-CB1+2/t0h^2*(u-v+3*w)*1/2*trace(MB1g*Grprp))*(-1/2*trace(MB1g*Grpr*(mA1g'+mA1g)*Grrp))...
        +2*real(2*(-1/2*conj(cA1)+1/t0h^2*(u-v+w)*1/2*trace(mA1g'*Grprp))*(-1/2*trace(mA1g*Grpr*(mA1g'+mA1g)*Grrp)))...
        +2*real(2*(-1/2*conj(cB1)+1/t0h^2*(u+v+w)*1/2*trace(mB1g'*Grprp))*(-1/2*trace(mB1g*Grpr*(mA1g'+mA1g)*Grrp)));
end
function gradRecB1_add=gradRecB1_add_func(R,CB1,CB2,CA2,cA1,cB1,cB2,Grprp,Grpr,Grrp)
    global r0 u v w t0h MA1g MB1g MB2g MA2g mA1g mB1g mB2g
    gradRecB1_add=2*(r0-R+2/t0h^2*(3*u+v+w)*1/2*trace(MA1g*Grprp))*(-1/2*trace(MA1g*Grpr*(mB1g'+mB1g)*Grrp))...
        +2*(-CB1+2/t0h^2*(u-v+3*w)*1/2*trace(MB1g*Grprp))*(-1/2*trace(MB1g*Grpr*(mB1g'+mB1g)*Grrp))...
        +2*real(2*(-1/2*conj(cA1)+1/t0h^2*(u-v+w)*1/2*trace(mA1g'*Grprp))*(-1/2*trace(mA1g*Grpr*(mB1g'+mB1g)*Grrp)))...
        +2*real(2*(-1/2*conj(cB1)+1/t0h^2*(u+v+w)*1/2*trace(mB1g'*Grprp))*(-1/2*trace(mB1g*Grpr*(mB1g'+mB1g)*Grrp)));
end
function gradImcA1_add=gradImcA1_add_func(R,CB1,CB2,CA2,cA1,cB1,cB2,Grprp,Grpr,Grrp)
    global r0 u v w t0h MA1g MB1g MB2g MA2g mA1g mB1g mB2g
    gradImcA1_add=2*(r0-R+2/t0h^2*(3*u+v+w)*1/2*trace(MA1g*Grprp))*(-1/2*trace(MA1g*Grpr*1i*(mA1g'-mA1g)*Grrp))...
        +2*(-CB1+2/t0h^2*(u-v+3*w)*1/2*trace(MB1g*Grprp))*(-1/2*trace(MB1g*Grpr*1i*(mA1g'-mA1g)*Grrp))...
        +2*real(2*(-1/2*conj(cA1)+1/t0h^2*(u-v+w)*1/2*trace(mA1g'*Grprp))*(-1/2*trace(mA1g*Grpr*1i*(mA1g'-mA1g)*Grrp)))...
        +2*real(2*(-1/2*conj(cB1)+1/t0h^2*(u+v+w)*1/2*trace(mB1g'*Grprp))*(-1/2*trace(mB1g*Grpr*1i*(mA1g'-mA1g)*Grrp)));
end

%%%%%%%%%%
%%% Main function free_energy(x)
%%% evaluates the free energy `y` and the corresponding gradient `grad` for a given configuration `x`.
%%% After every step, the current configuration and the current gradient are being plotted.
%%%%%%%%%%
function [y,grad]=free_energy(x)
    global N MA1g MB1g MB2g MA2g mA1g mB1g mB2g Nb eps count edgeR edgeCB1
    count=count+1;
    grad=zeros(1,4*(N-2));
    mcGm1=mcalGm1_func(x);
    %y=1/2*log(det(mcGm1));
    y=1/2*trace(logm(mcGm1));
    if abs(imag(y))>eps
        disp('shit imaginary det.');
    end
    
    mcG=inv(mcGm1);
    for r=2:(N-1)
        Grr=mcG((r-1)*Nb+1:r*Nb,(r-1)*Nb+1:r*Nb);
        y=y+Fadd_r_func(x((r-1)),x((N-2)+(r-1)),0,0,x(2*(N-2)+(r-1))+1i*0,x(3*(N-2)+(r-1))+1i*0,0,...
            1/2*trace(Grr*MA1g),1/2*trace(Grr*MB1g),1/2*trace(Grr*MB2g),1/2*trace(Grr*MA2g),...
            1/2*trace(Grr*mA1g),1/2*trace(Grr*mB1g),1/2*trace(Grr*mB2g));
        for rp=2:(N-1)
            Grpr=mcG((rp-1)*Nb+1:rp*Nb,(r-1)*Nb+1:r*Nb);
            Grprp=mcG((rp-1)*Nb+1:rp*Nb,(rp-1)*Nb+1:rp*Nb);
            Grrp=mcG((r-1)*Nb+1:r*Nb,(rp-1)*Nb+1:rp*Nb);
            grad(r-1)=grad(r-1)+real(gradR_add_func(x(rp-1),x((N-2)+(rp-1)),0,0,...
                x(2*(N-2)+(rp-1))+1i*0,x(3*(N-2)+(rp-1))+1i*0,0,Grprp,Grpr,Grrp));
            grad((N-2)+(r-1))=grad((N-2)+(r-1))+real(gradCB1_add_func(x(rp-1),x((N-2)+(rp-1)),0,0,...
                x(2*(N-2)+(rp-1))+1i*0,x(3*(N-2)+(rp-1))+1i*0,0,Grprp,Grpr,Grrp));
            grad(2*(N-2)+(r-1))=grad(2*(N-2)+(r-1))+real(gradRecA1_add_func(x(rp-1),x((N-2)+(rp-1)),0,0,...
                x(2*(N-2)+(rp-1))+1i*0,x(3*(N-2)+(rp-1))+1i*0,0,Grprp,Grpr,Grrp));
            grad(3*(N-2)+(r-1))=grad(3*(N-2)+(r-1))+real(gradRecB1_add_func(x(rp-1),x((N-2)+(rp-1)),0,0,...
                x(2*(N-2)+(rp-1))+1i*0,x(3*(N-2)+(rp-1))+1i*0,0,Grprp,Grpr,Grrp));         
        end
            Grpr=mcG((1-1)*Nb+1:1*Nb,(r-1)*Nb+1:r*Nb);
            Grprp=mcG((1-1)*Nb+1:1*Nb,(1-1)*Nb+1:1*Nb);
            Grrp=mcG((r-1)*Nb+1:r*Nb,(1-1)*Nb+1:1*Nb);
            grad(r-1)=grad(r-1)+real(gradR_add_func(edgeR,edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));
            grad((N-2)+(r-1))=grad((N-2)+(r-1))+real(gradCB1_add_func(edgeR,edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));
            grad(2*(N-2)+(r-1))=grad(2*(N-2)+(r-1))+real(gradRecA1_add_func(edgeR,edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));
            grad(3*(N-2)+(r-1))=grad(3*(N-2)+(r-1))+real(gradRecB1_add_func(edgeR,edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));

            Grpr=mcG((N-1)*Nb+1:N*Nb,(r-1)*Nb+1:r*Nb);
            Grprp=mcG((N-1)*Nb+1:N*Nb,(N-1)*Nb+1:N*Nb);
            Grrp=mcG((r-1)*Nb+1:r*Nb,(N-1)*Nb+1:N*Nb);
            grad(r-1)=grad(r-1)+real(gradR_add_func(edgeR,-edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));
            grad((N-2)+(r-1))=grad((N-2)+(r-1))+real(gradCB1_add_func(edgeR,-edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp)); 
            grad(2*(N-2)+(r-1))=grad(2*(N-2)+(r-1))+real(gradRecA1_add_func(edgeR,-edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));
            grad(3*(N-2)+(r-1))=grad(3*(N-2)+(r-1))+real(gradRecB1_add_func(edgeR,-edgeCB1,0,0,0,0,0,Grprp,Grpr,Grrp));            
    end
        Grr=mcG((1-1)*Nb+1:1*Nb,(1-1)*Nb+1:1*Nb);
        y=y+Fadd_r_func(edgeR,edgeCB1,0,0,0,0,0,1/2*trace(Grr*MA1g),1/2*trace(Grr*MB1g),1/2*trace(Grr*MB2g),...
            1/2*trace(Grr*MA2g),1/2*trace(Grr*mA1g),1/2*trace(Grr*mB1g),1/2*trace(Grr*mB2g));
        Grr=mcG((N-1)*Nb+1:N*Nb,(N-1)*Nb+1:N*Nb);
        y=y+Fadd_r_func(edgeR,-edgeCB1,0,0,0,0,0,1/2*trace(Grr*MA1g),1/2*trace(Grr*MB1g),1/2*trace(Grr*MB2g),...
            1/2*trace(Grr*MA2g),1/2*trace(Grr*mA1g),1/2*trace(Grr*mB1g),1/2*trace(Grr*mB2g));
    
    if abs(imag(y))<eps
        y=real(y);
    else
        disp('shit imaginary.');
    end
    if mod(count,1)==0
        Rbar=0;CB1bar=0;RecA1bar=0;RecB1bar=0;gradR=0;gradCB1=0;gradRecA1bar=0;gradRecB1bar=0;
        for i=1:(N-2)
            Rbar=Rbar+x(i)/(N-2);
            CB1bar=CB1bar+x((N-2)+i)/(N-2);
            RecA1bar=RecA1bar+x(2*(N-2)+i)/(N-2);
            RecB1bar=RecB1bar+x(3*(N-2)+i)/(N-2);
            gradR=gradR+abs(grad(i))/(N-2);
            gradCB1=gradCB1+abs(grad((N-2)+i))/(N-2);
            gradRecA1bar=gradRecA1bar+abs(grad(2*(N-2)+i))/(N-2);
            gradRecB1bar=gradRecB1bar+abs(grad(3*(N-2)+i))/(N-2);

        end
        dispX=['F=',num2str(y),' R=',num2str(Rbar),' CB1=',num2str(CB1bar),' cA1=',num2str(RecA1bar),...
            ' cB1=',num2str(RecB1bar),' dR=',num2str(gradR),' dCB1=',num2str(gradCB1),' dcA1=',num2str(gradRecA1bar),' dcB1=',num2str(gradRecB1bar)];
        disp(dispX);
        subplot(2,1,1);
        plot(1:(N-2),x(1:(N-2)),1:(N-2),x((N-2)+1:2*(N-2)),...
            1:(N-2),x(2*(N-2)+1:3*(N-2)),1:(N-2),x(3*(N-2)+1:4*(N-2)));
        legend('R','CB1','cA1','cB1');
        subplot(2,1,2);
        plot(1:(N-2),grad(1:(N-2)),1:(N-2),grad((N-2)+1:2*(N-2)),...
            1:(N-2),grad(2*(N-2)+1:3*(N-2)),1:(N-2),grad(3*(N-2)+1:4*(N-2)));
        drawnow;       
    end
end


