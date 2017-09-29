pumavgvar = [];
rfcavgvar = [];
pepavgvar = [];
pumpatchratio=[];
rfcpatchratio=[];
peppatchratio=[];

pumnumsmall=0;
rfcnumsmall=0;
pepnumsmall=0;
pumnumpatches = 0;
rfcnumpatches = 0;
pepnumpatches = 0;
for n=1:numnames
    
    pumweight = [];
    rfcweight = [];
    pepweight = [];
    pumvariation = [];
    rfcvariation = [];
    pepvariation = [];
    
    
    pumnumsmall=0;
    rfcnumsmall=0;
    pepnumsmall=0;
    
    if toanalyze(n)
        
        if n==11
           n 
        end
        
        for i=1:length(pum_opatches{n})        
            pump = pum_opatches{n}{i};
            if length( pump ) > 4
                pump = pump(pump ~=-100);
                [maxdiff maxind] = max(pump);

                alldists = pdist2(pump, pump(maxind));

                if any( abs( alldists ) > 30 )
%                     var(pump)
                    pump( abs( alldists ) > 30 ) = pump( abs( alldists ) > 30 ) +60;
                    pump( abs( alldists ) <= 30 ) = pump( abs( alldists ) <= 30 ) + alldists(abs( alldists ) <= 30).*2;
%                     var(pump)
                end
                
                pumvariation = [pumvariation var(pump)];
                pumweight = [pumweight length(pump)];
            else
                pumnumsmall = pumnumsmall+1;
            end
        end

        for i=1:length(rfc_opatches{n}) 
            rfcp = rfc_opatches{n}{i};
            if length( rfcp ) > 4
                rfcp = rfcp(rfcp ~=-100);
                % Phase unwrap because we're doing std dev and need direct
                % distances
                [maxdiff maxind] = max(rfcp);

                alldists = pdist2(rfcp, rfcp(maxind));

                if any( abs( alldists ) > 30 )
%                     disp([ 'Before: ' num2str(var(rfcp)) ] );
                    rfcp( abs( alldists ) > 30 ) = rfcp( abs( alldists ) > 30 ) +60;
                    rfcp( abs( alldists ) <= 30 ) = rfcp( abs( alldists ) <= 30 ) + alldists(abs( alldists ) <= 30).*2;
%                     disp([ 'After: ' num2str(var(rfcp)) ] );
                end
                
                rfcvariation = [rfcvariation var(rfcp)];
                rfcweight = [rfcweight length(rfcp)];
            else
                rfcnumsmall = rfcnumsmall+1;
            end
        end

        
        
        for i=1:length(pep_opatches{n}) 
            pepp = pep_opatches{n}{i};
            if length( pepp ) > 4
                pepp = pepp(pepp ~=-100);
                % Phase unwrap because we're doing std dev and need direct
                % distances
                [maxdiff maxind] = max(pepp);

                alldists = pdist2(pepp, pepp(maxind));

                if any( abs( alldists ) > 30 )
                    pepp( abs( alldists ) > 30 ) = pepp( abs( alldists ) > 30 ) +60;
                    pepp( abs( alldists ) <= 30 ) = pepp( abs( alldists ) <= 30 ) + alldists(abs( alldists ) <= 30).*2;
                end
                
                pepvariation = [pepvariation var(pepp)];
                pepweight = [pepweight length(pepp)];
            else
                pepnumsmall = pepnumsmall+1;
            end
        end

        top=0;
        bottom=0;
        for i=1:length(pumvariation)
             top = top + ((pumweight(i)-1) * pumvariation(i));
             bottom = bottom + pumweight(i);
        end
        if length(pumvariation) > 0
            pumavgvar = [pumavgvar; top/bottom];
%             pumpatchratio = [pumpatchratio; pumnumsmall/length(pum_opatches{n})];
%             pumnumpatches = pumnumpatches+length(pum_opatches{n});
        end
        
        top=0;
        bottom=0;
        for i=1:length(rfcvariation)
             top = top + ((rfcweight(i)-1) * rfcvariation(i));
             bottom = bottom + rfcweight(i);
        end

        if length(rfcvariation) > 0
            rfcavgvar = [rfcavgvar; top/bottom];
        end
        
        top=0;
        bottom=0;
        for i=1:length(pepvariation)
             top = top + ((pepweight(i)-1) * pepvariation(i));
             bottom = bottom + pepweight(i);
        end
        if length(pepvariation) > 0
            pepavgvar = [pepavgvar; top/bottom];
%             peppatchratio = [peppatchratio; pepnumsmall/length(pep_opatches{n})];
%             pepnumpatches = pepnumpatches+length(pep_opatches{n});
        end
                 
    end
end


pumpooledvar = mean(pumavgvar)
rfcpooledvar = mean(rfcavgvar)
peppooledvar = mean(pepavgvar)

% pumavgstddevcorrected = mean(pumavgvar.*pumpatchratio)
% rfcavgstddevcorrected = mean(rfcavgvar.*rfcpatchratio)
% pepavgstddevcorrected = mean(pepavgvar.*peppatchratio)

% pumnumpatches
% pumnumsmall
% pumnumsmall/pumnumpatches
% 
% rfcnumpatches
% rfcnumsmall
% rfcnumsmall/rfcnumpatches
% 
% pepnumpatches
% pepnumsmall
% pepnumsmall/pepnumpatches
