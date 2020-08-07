function vectorRightHandSide =  compute_vector_rhs(simulation,E_0,E,E_pts,fieldType)
diff_array = zeros(size(E),'single');

switch fieldType
    case 'x'
        diff_array(:,1) = -(abs(E_0(:,1)).^2-abs(E(:,1)).^2).*conj(E);
    case 'y'
        diff_array(:,2) = -(abs(E_0(:,2)).^2-abs(E(:,2)).^2).*conj(E);
    case 'z'
        diff_array(:,3) = -(abs(E_0(:,3)).^2-abs(E(:,3)).^2).*conj(E);
    case 'xy'
        diff_array(:,1) = abs(E_0(:,1)).^2-abs(E(:,1)).^2;
        diff_array(:,2) = abs(E_0(:,2)).^2-abs(E(:,2)).^2;
        diff_array = -diff_array.*conj(E);        
    case 'yz'
        diff_array(:,2) = abs(E_0(:,2)).^2-abs(E(:,2)).^2;
        diff_array(:,3) = abs(E_0(:,3)).^2-abs(E(:,3)).^2;
        diff_array = -diff_array.*conj(E);        
    case 'xz'
        diff_array(:,1) = abs(E_0(:,1)).^2-abs(E(:,1)).^2;
        diff_array(:,3) = abs(E_0(:,3)).^2-abs(E(:,3)).^2;
        diff_array = -diff_array.*conj(E);
    case 'xyz'
        diff_array = -(abs(E_0).^2-abs(E).^2).*conj(E);
    otherwise
        print('specify case');
end

vectorRightHandSide = compute_adjoint_rhs(simulation,diff_array,E_pts);

