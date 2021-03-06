function mtexColorMap(name,varargin)
% define an MTEX colormap

if ischar(name)
  if isempty(which([name '.m']))
    name = [name,'ColorMap'];
    if isempty(which(name))
      error('unknown colormap name');      
    end
    map = feval(name);
  else
    map = colormap(name);
  end
else
  map = name;
end

if isappdata(gcf,'mtexFig')
  mtexFig = getappdata(gcf,'mtexFig');
  for i = 1:numel(mtexFig.children)
    colormap(mtexFig.children(i),map);
  end  
else
  colormap(map);
end

  
