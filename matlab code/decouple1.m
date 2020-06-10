%Initialization
Vm=0;
delta=0;
yload=0;
deltad=0;

%Input taken from Bus
for k=1:nbus
    %For bus
    n=busdata(k,1);
    kb(n)=busdata(k,2);
    Vm(n)=busdata(k,3);
    delta(n)=busdata(k,4);
    
    %For Load
    Pd(n)=busdata(k,5);
    Qd(n)=busdata(k,6);
    
    %For generator
    Pg(n)=busdata(k,7);
    Qg(n) = busdata(k,8);
    Qmin(n)=busdata(k, 9);
    Qmax(n)=busdata(k, 10);
    %Static(For Power factor correction)
    Qsh(n)=busdata(k, 11);
    
    if Vm(n) <= 0 %For error happened
        Vm(n) = 1.0;
        V(n) = 1 + j*0;
    else
        
     %Calculation from given  bus data
        delta(n) = pi/180*delta(n);
        V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
        P(n)=(Pg(n)-Pd(n))/basemva;  %In per-unit
        Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;         
        S(n) = P(n) + j*Q(n);
    end
end

%Data taken from YBUS
Ym = abs(Ybus);
t = angle(Ybus);
ii=0;

for ib=1:nbus
    if kb(ib) == 0 | kb(ib) == 2   %As Slack bus is used as reference bus , michmatch is not 
        ii = ii+1;                  %is not defined in Slack bus
        jj=0;
        for jb=1:nbus
            if kb(jb) == 0 | kb(jb) == 2
                jj = jj+1;
                B1(ii,jj)=imag(Ybus(ib,jb));%B1 matrix for calculating del(angle)
            end                              %B1 is imaginary of YBUS
        end
     end
end

ii=0;
for ib=1:nbus
    if kb(ib) == 0    %Voltage will not change for Generator Bus
        ii = ii+1;
        jj=0;
        for jb=1:nbus
            if kb(jb) == 0
                jj = jj+1;
                B2(ii,jj)=imag(Ybus(ib,jb)); %B2 matrix for calculating del(Voltage)
            end
        end
     end
end

%Initial Value
maxerror = 1;
converge = 1;
iter = 0;
% Start of iterations
while maxerror >= accuracy & iter <= maxiter  %Test for max. power mismatch
    iter = iter+1;
    id=0;   %(Iteraration number for del(angle)
    iv=0;   %(Iteraration number for del(voltage)
    
    for n=1:nbus
        J11=0;    %Jacobian Matrix initialization
        J22=0;
        
        for i=1:nbr   %nbr=total num of lines
            if nl(i) == n | nr(i) == n
                if nl(i) == n, 
                    l = nr(i);
                end
                if nr(i) == n,  
                    l = nl(i);
                end
                
                %Update for Jacobian Matrix(Non-Diagonal)
                J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));  
                J22=J22+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
             end
        end
        
        %Calculated Real & Reactive Power(Diagonal)
        Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J22;
        Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
        
        if kb(n) == 1    %Calculation for Slack Bus
            P(n)=Pk;
            Q(n) = Qk;
        end   
        if kb(n) == 2   %For Generator Bus
            Q(n)=Qk;
            Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
            if Qmax(n) ~= 0
                if iter <= 20                        % Between the 10th & 20th iterations
                    if iter >= 10                    % the Mvar of generator buses are
                        if Qgc  < Qmin(n),           % tested. If not within limits Vm(n)
                            Vm(n) = Vm(n) + 0.005;   % is changed in steps of 0.005 pu to
                        elseif Qgc > Qmax(n),        % bring the generator Mvar within
                            Vm(n) = Vm(n) - 0.005;
                        end                          % the specified limits.
                     end
                end
             end
        end
        if kb(n) ~= 1
            id = id+1;
            DP(id) = P(n)-Pk;           %Calculating del(P)
            DPV(id) = (P(n)-Pk)/Vm(n);  %del(p)/(Voltage)
        end
        
        if kb(n) == 0
            iv=iv+1;
            DQ(iv) = Q(n)-Qk;
            DQV(iv) = (Q(n)-Qk)/Vm(n);  %Calculating del(Q)/(Volt.)
        end
    end
    Dd=-B1\DPV';  %Del(angle)
    DV=-B2\DQV';  %Del(Volt.)
    id=0;
    iv=0;
    for n=1:nbus
        if kb(n) ~= 1      %Not for Slack Bus
            id = id+1;
            delta(n) = delta(n)+Dd(id);   %Updating value of angle
        end
        if kb(n) == 0                %For only Load Bus
            iv = iv+1;
            Vm(n)=Vm(n)+DV(iv);      %Updating the value of voltage
        end
    end
    maxerror=max(max(abs(DP)),max(abs(DQ)));
    if iter == maxiter & maxerror > accuracy
        fprintf('\nWARNING: Iterative solution did not converged after ')
        fprintf('%g', iter), fprintf(' iterations.\n\n')
        fprintf('Press Enter to terminate the iterations and print the results \n')
        converge = 0; pause,
    end
    
end
if converge ~= 1
    tech= ('ITERATIVE SOLUTION DID NOT CONVERGE');
else 
    tech=('Power Flow Solution by Fast Decoupled Method');
end

k=0;
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;

for n = 1:nbus
    if kb(n) == 1    %Calculated Values for Slack Bus
        S(n)=P(n)+j*Q(n);
        Pg(n) = P(n)*basemva + Pd(n);   
        Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
        k=k+1;
        Pgg(k)=Pg(n);
        
    elseif  kb(n) ==2
        S(n)=P(n)+j*Q(n);
        Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
        k=k+1;
        Pgg(k)=Pg(n);
    end
    yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm';
busdata(:,4)=deltad';
Pgt = sum(Pg);
Qgt = sum(Qg);
Pdt = sum(Pd);   %Total
Qdt = sum(Qd);
Qsht = sum(Qsh);

