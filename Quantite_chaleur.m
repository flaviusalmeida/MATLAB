function Q = Quantite_chaleur(param_teq,Ciment,teq,i,n)

%
C = param_teq(1) * Ciment * 1.E3 ; % paramètre C' [J/m3] : C' = C [J/g] * qté liant [kg/m3] * 1000
Q1 = param_teq(2) * Ciment * 1.E3  ; % paramètre Q1 [J/m3] : C' = C [J/g] * qté liant [kg/m3] * 1000
t1 = param_teq(3) ; % paramètre t1 [heures]
dt1 = param_teq(4) ; % paramètre dt1 [heures]
Q2 = param_teq(5) * Ciment * 1.E3  ; % paramètre Q2 [J/m3] : C' = C [J/g] * qté liant [kg/m3] * 1000
t2 = param_teq(6) ; % paramètre t2 [heures]
dt2 = param_teq(7) ; % paramètre dt2 [heures]

Q = C + ( Q1/(1 + exp ((teq(i,n)/3600. - t1)/dt1) ) ) + ( Q2/(1 + exp ((teq(i,n)/3600. - t2)/dt2) ) );

end