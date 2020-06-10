j=sqrt(-1);
i = sqrt(-1);

%Input from Linedata
nl = linedata(:,1);%from line
nr = linedata(:,2);%To line
R = linedata(:,3);
X = linedata(:,4);
Bc = j*linedata(:,5); %Shunt Admittance
a = linedata(:, 6); %Tranformer

nbr=length(linedata(:,1));  %Total number of line
nbus = max(max(nl), max(nr)); %Total num of bus
Z = R + j*X; y= ones(nbr,1)./Z;  %Series admittance

for n = 1:nbr
    if a(n) <= 0
        a(n) = 1;
    end
    
    Ybus=zeros(nbus,nbus);     % initialize Ybus to zero
    %Formation of the off diagonal elements
    for k=1:nbr;
        Ybus(nl(k),nr(k))=Ybus(nl(k),nr(k))-y(k)/a(k);
        Ybus(nr(k),nl(k))=Ybus(nl(k),nr(k));
    end
end

%Formation of the diagonal elements
for  n=1:nbus
    for k=1:nbr
        if nl(k)==n
            Ybus(n,n) = Ybus(n,n)+y(k)/(a(k)^2) + Bc(k);
        elseif nr(k)==n
            Ybus(n,n) = Ybus(n,n)+y(k) +Bc(k);
         end
    end
end
