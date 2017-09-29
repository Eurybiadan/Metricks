for i=1:length(fNames)
    
    fNames{i}
    % For comparing submosaics and real mosaics
    % Find all cases where they match
    matches = logical((RFC{i}~=-100) & (PUM{i}~=-100) & (PEPPE{i}~=-100));
    rfc = (RFC{i}(matches));
    pum = (PUM{i}(matches));
    pep = (PEPPE{i}(matches));        
    
    % For single, rotated mosaic            
    [tok remain]=strtok(fNames{i},'_');
    [tok remain]=strtok(remain,'_');

    refangle = str2double(tok);
    
    rfc = RFC{i}(RFC{i}~=-100);
    pum = PUM{i}(PUM{i}~=-100);
    pep = PEPPE{i}(PEPPE{i}~=-100);
    
    tmp = abs(rfc-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
    
%     if any(tmp >3)
%        tmp; 
%     end
    rfcabsdist{i} = tmp;
    rfcabserr(i) = mean( tmp );
    
    tmp = abs(pum-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     if any(tmp >3)
%        tmp; 
%     end
    pumabsdist{i} = tmp;
    pumabserr(i) = mean( tmp );
    
    tmp = abs(pep-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
    if any(tmp >3)
       tmp; 
    end
    pepabsdist{i} = tmp;
    pepabserr(i) = mean( tmp );

 
end

rfcarray=[];
pumarray=[];
peparray=[];
for i=1:length(rfcabsdist)
    rfcarray = [rfcarray; rfcabsdist{i}];
    pumarray = [pumarray; pumabsdist{i}];
    peparray = [peparray; pepabsdist{i}];
end
    
percentrfcright = (1- (sum(rfcarray(:)>3)/length(rfcarray(:))))*100
percentpumright = (1- (sum(pumarray(:)>3)/length(pumarray(:))))*100
percentpepright = (1- (sum(peparray(:)>3)/length(peparray(:))))*100

avgrfcerr = mean(rfcabserr)
avgpumerr = mean(pumabserr)
avgpeperr = mean(pepabserr)