function [info]= compute_info_measures2(pdf)
% Compute either entropy, conditional entropy, mutual
% information, lagged mutual information, or transfer entropy
% use KDE method, Epanichov (Silverman 1986) kernel to compute 1d, 2d, 3d
% pdfs

%1D input pdf: Hx (1 output)
%2D input pdf: H(x1), H(x2), H(x1|x2), H(x2|x1), I(x1;x2) (5 outputs)
%3D input pdf

%Updated 9/8/15 - exclude redundancy when sources are uncorrelated
%Updated 9/30/15 - output as structure instead of individual arguments

dim = length(size(pdf)); %dimension of pdf

v = size(pdf);
if sum(v==1)>0 %if one of the dimensions = 1 (1d pdf)
    dim=dim-1;
end

N = size(pdf,1);

if dim ==1 %compute: entropy H(X)
    Hvect=pdf.*log2(1./pdf);
    Hvect(isnan(Hvect))=0;
    Hx = sum(Hvect);
    info.Hx = Hx;
end

if dim == 2 %compute: H(X),H(Y),H(X|Y),H(Y|X),I(X;Y)
    H_xgy=0; H_ygx=0;
    
    m_i = sum(pdf,2);
    Hivect = m_i.*log2(1./m_i);
    Hivect(isnan(Hivect))=0;
    Hx = sum(Hivect);
    
    m_j = sum(pdf,1);
    Hjvect = m_j.*log2(1./m_j);
    Hjvect(isnan(Hjvect))=0;
    Hy = sum(Hjvect);
    
    for i = 1:N      
        for j = 1:N
            m_ij = pdf(i,j); % joint prob
            
            mj = m_j(j);
            mi = m_i(i);
            
            if m_ij > 0 && mi>0 %H(Y|X)=sum of p(i,j)*log(p(i)/p(i,j))
                H_ygx = H_ygx + m_ij*log2(mi/m_ij);
            end
            
            if m_ij>0 && mj >0 %H(X|Y)=sum of p(i,j)*log(p(j)/p(i,j))
                H_xgy = H_xgy + m_ij*log2(mj/m_ij);
            end
         
        end %end i loop
    end     %end j loop
    
    info.Hx1 = Hx;
    info.Hx2 = Hy;
    info.H_xgy = H_xgy;
    info.H_ygx = H_ygx;
    info.I = min(Hx-H_xgy,Hy-H_ygx);
end

if dim ==3  %compute: H(X), H(Y), H(Z),
    %Data [Xt Yw Yf] --> assume last dimension is the target, first is the source
    I_x1y=0; I_x2y=0; I_x1x2=0; T = 0;
    
    m_jk = reshape(sum(pdf(:,:,:),1),N,N);
    m_ij = reshape(sum(pdf(:,:,:),3),N,N);
    m_ik = reshape(sum(pdf(:,:,:),2),N,N);
    
    m_i = reshape(sum(m_ij(:,:),2),1,N);
    m_j = reshape(sum(m_ij(:,:),1),1,N);
    m_k = reshape(sum(m_jk(:,:),1),1,N);
    
    Hivect = m_i.*log2(1./m_i);
    Hivect(isnan(Hivect))=0;
    Hx1 = sum(Hivect);
    
    Hjvect = m_j.*log2(1./m_j);
    Hjvect(isnan(Hjvect))=0;
    Hx2 = sum(Hjvect);
    
    Hkvect = m_k.*log2(1./m_k);
    Hkvect(isnan(Hkvect))=0;
    Hy = sum(Hkvect);
    
    for i=1:N
        for j=1:N
            for k=1:N
                m_ijk=pdf(i,j,k);
                
                %Transfer Entropy = T = U_i + S
                if m_ijk > eps && m_ij(i,j)> eps && m_jk(j,k) >eps && m_j(j)>eps
                    T_add = m_ijk*log2((m_ijk*m_j(j))/(m_ij(i,j)*m_jk(j,k)));
                    if T_add>0
                        T=T+T_add;
                    end
                end
                
                %Lagged Mutual Info I(xt to yf) = U_i + R
                if j==1
                    I_tau_add = m_ik(i,k)*log2(m_ik(i,k)/(m_i(i)*m_k(k)));
                    if I_tau_add>0
                        I_x1y = I_x1y+I_tau_add;
                    end
                end
                
                %Lagged Mutual Info (yw to yf) = U_j + R
                if i==1
                    I_tau2_add = m_jk(j,k)*log2(m_jk(j,k)/(m_j(j)*m_k(k)));
                    if I_tau2_add>0
                        I_x2y = I_x2y+I_tau2_add;
                    end
                end
                
                %Lagged Mutual Info (xt to yw) (source dependency)
                if k==1
                    I_tau3_add = m_ij(i,j)*log2(m_ij(i,j)/(m_i(i)*m_j(j)));
                    if I_tau3_add>0
                        I_x1x2 = I_x1x2+I_tau3_add;
                    end
                end
        
            end
        end
    end
    
   
    
    dI = T-I_x1y; %interaction information = S-R (if positive: synergy dominates)
    I_tot = dI + I_x1y + I_x2y; %total shared information (Ui + Uj + R + S)
    
    % new formulation of redundancy: normalize I between sources (I_x1x2)
    % multiply R by normalized I (to decrease I when sources independent)
    I_sourcenorm = I_x1x2./min(Hx1,Hx2);
    %R = R .* I_sourcenorm; %if I(S1;S2)=0, R=0
    

% Account for source correlation and keep redundancy within bounds
% I_x1y + I_x2y - I_tot < R < min[I_x1y,I_x2y]
Rmax = min(I_x1y,I_x2y);
Rmin = max(0,I_x1y+I_x2y-I_tot);
dR = Rmax-Rmin;
R = Rmin +dR.*I_sourcenorm;

if R<0
    fprintf('negative R!\n')
    fprintf('Rmin = %4.3f, Rmax = %4.3f\n',Rmin,Rmax)
end

    
    
    %compute synergy and unique contributions from nodes from R
    %I_x1y = U1 + R,  I_x2y = U2 + R,    di = S-R
    U1 = I_x1y-R;
    U2 = I_x2y-R;
    S = I_tot - (U1+U2+R); %(S-R)+R = S
    
    info.Hx1 = Hx1;
    info.Hx2 = Hx2;
    info.Hy = Hy;
    info.I_x1y = I_x1y;
    info.I_x2y = I_x2y;
    info.T =T; %I(Xt;Yf|Yw)
    info.dI = dI;
    info.I_tot = I_tot;
    info.R = R;
    info.S = S;
    info.U1 = U1;
    info.U2 = U2;
    
end

end
% end of function
