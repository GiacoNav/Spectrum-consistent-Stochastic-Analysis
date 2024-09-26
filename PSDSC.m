function [Gw] = PSDSC(w,G0,ex,wx)
%funzione per il calcolo punto a punto della PSD spettrocompatibile
if w<wx(1), Gw=G0*(wx(1)/wx(2))^ex(2)*(w/wx(1))^ex(1); end
if (w>=wx(1)&&w<wx(2)), Gw=G0*(w/wx(2))^ex(2); end
if (w>=wx(2)&&w<wx(3)), Gw=G0*(w/wx(2))^ex(3); end
if (w>=wx(3)), Gw=G0*(wx(3)/wx(2))^ex(3)*(w/wx(3))^ex(4);
end

