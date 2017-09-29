rfcall = [];
pumall = [];
pepall = [];
% Analyze their differences.
for i=1:length(fNames)

    i;
    % For single, rotated mosaic
%     [tok remain]=strtok(fNames{i},'_');
%     [tok remain]=strtok(remain,'_');
% 
%     refangle = str2double(tok);
    
%     rfc = RFC{i}(RFC{i}~=-100);
%     pum = PUM{i}(PUM{i}~=-100);
%     pep = PEPPE{i}(PEPPE{i}~=-100);
 
    % For comparing submosaics and real mosaics
    % Find all cases where they match
    matches = logical((RFC{i}~=-100) & (PUM{i}~=-100) & (PEPPE{i}~=-100));
    rfc = (RFC{i}(matches));
    pum = (PUM{i}(matches));
    pep = (PEPPE{i}(matches));
    
    
    
    diffrfcpum = (rfc-pum);
    diffrfcpep = (rfc-pep);
    diffpeppum = (pep-pum);
    
    diffs = [diffrfcpum diffrfcpep diffpeppum];
    phasewrapped = abs(diffs) > 30;
    
    % Determine our biggest difference- this will be our anchor
    [vals inds]=max( abs( diffs ),[],2);
    
    % Shift the difference to the anchor
    for j=1:length(inds)
        
        if phasewrapped(j, inds(j))
        switch inds(j)
            case 1 % between rfc and pum
                if sign(diffrfcpum(j)) == 1 % rfc is higher, and the anchor
                    pum(j) = pum(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        pep(j) = pep(j)+60;
                    else                                                
                        pep(j) = diffrfcpep(j).*2+pep(j);
                    end
                    
                else % pum higher, and is the anchor- bring rfc to it
                    rfc(j) = rfc(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        pep(j) = pep(j)+60;
                    else
                        pep(j) = -diffpeppum(j).*2+pep(j);
                    end
                    
                end

            case 2 % between rfc and pep
                if sign(diffrfcpep(j)) == 1 % rfc is the anchor
                    pep(j) = pep(j) +60;
                                     
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        pum(j) = pum(j)+60;
                    else                      
                        pum(j) = diffrfcpum(j).*2+pum(j);                        
                    end                   
                    
                else % pep is the anchor
                    rfc(j) = rfc(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        pum(j) = pum(j)+60;
                    else
                        pum(j) = diffpeppum(j).*2+pum(j);
                    end 
                    
                end
            case 3 % between pep and pum
                if sign(diffpeppum(j)) == 1 % pep is the anchor
                    pum(j) = pum(j) +60;
                                       
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        rfc(j) = rfc(j)+60;
                    else
                        rfc(j) = -diffrfcpep(j).*2+rfc(j);
                    end                   
                    
                else % pum is the anchor
                    pep(j) = pep(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        rfc(j) = rfc(j)+60;
                    else
                        rfc(j) = -diffrfcpum(j).*2+rfc(j);
                    end 
                    
                end
        end
        end
    end
    
%     if i == 1
%         i
%     end
% 
%     
%     if any( abs(rfc-pum) > 30 )
%         [val ind]=max(abs(rfc-pum))
%         i;
%     end
%     
%     if any( abs(rfc-pep) > 30 )
%         [val ind]=max(abs(rfc-pep))
%         i;
%     end
%     
%     if any( abs(pep-pum) > 30 )
%         [val ind]=max(abs(pep-pum))
%         i;
%     end
%     

    rfcall = [rfcall; rfc];
    pumall = [pumall; pum];
    pepall = [pepall; pep];
    
    
    
    
    % For single, rotated mosaic
%     tmp = abs(rfc-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     rfcabsdist(i,:) = tmp;
%     rfcabserr(i) = mean( tmp );
%     
%     tmp = abs(pum-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     pumabsdist(i,:) = tmp;
%     pumabserr(i) = mean( abs(pum-refangle) );
%     
%     tmp = abs(pep-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     pepabsdist(i,:) = tmp;
%     pepabserr(i) = mean( abs(pep-refangle) );

 
end

alldata = [pepall, rfcall, pumall];
alldata0 = mean(alldata);
centereddata=bsxfun(@minus,alldata,alldata0);
% [u,s,v]= svd(centereddata);
% 
% direction = -v(:,1);
% 
plot3(centereddata(:,1),centereddata(:,2),centereddata(:,3),'.'); hold on;
% quiver3(0,0,0,direction(1)*60,direction(2)*60,direction(3)*60,'r');
% quiver3(0,0,0,-direction(1)*60,-direction(2)*60,-direction(3)*60,'r');
hold off;
