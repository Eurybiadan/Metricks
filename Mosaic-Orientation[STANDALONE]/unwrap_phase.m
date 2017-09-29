function [ radon_unwrapped,pum_unwrapped,fourier_unwrapped, dists ] = unwrap_phase( RADON, PUM, FOURIER )


% For comparing submosaics and real mosaics
    % Find all cases where they match
    matches = logical((RADON~=-100) & (PUM~=-100) & (FOURIER~=-100));
    radon_unwrapped = (RADON(matches));
    pum_unwrapped = (PUM(matches));
    fourier_unwrapped = (FOURIER(matches));
    
    diffrfcpum = (radon_unwrapped-pum_unwrapped);
    diffrfcpep = (radon_unwrapped-fourier_unwrapped);
    diffpeppum = (fourier_unwrapped-pum_unwrapped);
    
    diffs = [diffrfcpum diffrfcpep diffpeppum];
    phasewrapped = abs(diffs) > 30;
    
%     Determine our biggest difference- this will be our anchor
    [vals inds]=max( abs( diffs ),[],2);
    
    % Shift the difference to the anchor
    for j=1:length(inds)
        
        if phasewrapped(j, inds(j))
        switch inds(j)
            case 1 % between rfc and pum
                if sign(diffrfcpum(j)) == 1 % rfc is higher, and the anchor
                    pum_unwrapped(j) = pum_unwrapped(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        fourier_unwrapped(j) = fourier_unwrapped(j)+60;
                    else                                                
                        fourier_unwrapped(j) = diffrfcpep(j).*2+fourier_unwrapped(j);
                    end
                    
                else % pum higher, and is the anchor- bring rfc to it
                    radon_unwrapped(j) = radon_unwrapped(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        fourier_unwrapped(j) = fourier_unwrapped(j)+60;
                    else
                        fourier_unwrapped(j) = -diffpeppum(j).*2+fourier_unwrapped(j);
                    end
                    
                end

            case 2 % between rfc and pep
                if sign(diffrfcpep(j)) == 1 % rfc is the anchor
                    fourier_unwrapped(j) = fourier_unwrapped(j) +60;
                                     
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        pum_unwrapped(j) = pum_unwrapped(j)+60;
                    else                      
                        pum_unwrapped(j) = diffrfcpum(j).*2+pum_unwrapped(j);                        
                    end                   
                    
                else % pep is the anchor
                    radon_unwrapped(j) = radon_unwrapped(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        pum_unwrapped(j) = pum_unwrapped(j)+60;
                    else
                        pum_unwrapped(j) = diffpeppum(j).*2+pum_unwrapped(j);
                    end 
                    
                end
            case 3 % between pep and pum
                if sign(diffpeppum(j)) == 1 % pep is the anchor
                    pum_unwrapped(j) = pum_unwrapped(j) +60;
                                       
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        radon_unwrapped(j) = radon_unwrapped(j)+60;
                    else
                        radon_unwrapped(j) = -diffrfcpep(j).*2+radon_unwrapped(j);
                    end                   
                    
                else % pum is the anchor
                    fourier_unwrapped(j) = fourier_unwrapped(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        radon_unwrapped(j) = radon_unwrapped(j)+60;
                    else
                        radon_unwrapped(j) = -diffrfcpum(j).*2+radon_unwrapped(j);
                    end 
                    
                end
        end
        end
        
        dists(j,:) = abs([radon_unwrapped(j)-pum_unwrapped(j) radon_unwrapped(j)-fourier_unwrapped(j) pum_unwrapped(j)-fourier_unwrapped(j)]);
    end

end

