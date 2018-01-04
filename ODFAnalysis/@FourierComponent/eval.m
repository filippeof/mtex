function f = eval(component,ori,varargin)
% evaluate odf using NSOFT
%
% Input
%  odf - @ODF
%  ori - @orientation
% 
% Output
%  f - double
%

if isempty(ori), f = []; return; end

% TODO: this can be done better!!!
if component.antipodal || check_option(varargin,'antipodal')
  varargin = delete_option(varargin,'antipodal');
  component.antipodal = false;
  f = 0.5 * (eval(component,ori,varargin{:}) + eval(component,inv(ori),varargin{:}));
  return
end

% set parameter
L = dim2deg(length(component.f_hat));
L = int32(min(L,get_option(varargin,'bandwidth',L)));
Ldim = deg2dim(double(L+1));

% export to Euler angle
g = Euler(ori,'nfft'); %[alpha';beta';gamma'];
	
f_hat = [real(component.f_hat(1:Ldim)),imag(component.f_hat(1:Ldim))].';

% run NFSOFT
if getMTEXpref('extern',false)

  f = call_extern('fc2odf','intern',L,'EXTERN',g,f_hat);
  
else

  L = double(L);
  % create plan
  plan = nfsoftmex('init',L,length(ori),0,0,6,1000,3*L);
    
  % set nodes
  nfsoftmex('set_x',plan,flipud(g));

  % node-dependent precomputation
  nfsoftmex('precompute',plan);

  fhat = component.f_hat(1:Ldim);
  for l = 1:L
    
    [k1,k2] = meshgrid(-l:l,-l:l);
    k1(k1>0) = 0;
    k2(k2>0) = 0;
    s = (-1).^k1 .* (-1).^k2;
    
    ind = (deg2dim(l)+1):deg2dim(l+1);
    fhat(ind) = s.*reshape(fhat(ind),2*l+1,2*l+1);
  end
  
  % set Fourier coefficients
  nfsoftmex('set_f_hat',plan,fhat);
  
  % transform
  nfsoftmex('trafo',plan);

  % get function values
  f = real(nfsoftmex('get_f',plan));

  % kill plan
  nfsoftmex('finalize',plan);
  
end
