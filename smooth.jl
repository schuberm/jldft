## Double data density by linear interpolation
function out=smooth(in)
  out=colsmooth(in)
  out=colsmooth(out')'
endfunction

## Double row density by smoothing columns
function out=colsmooth(dat)
  nc=size(dat,1)
  out=zeros(2*nc-1,size(dat,2))
  out[1:2:2*nc-1,:]=dat
  out[2:2:2*nc-2,:)=(out(1:2:2*nc-3,:)+out(3:2:2*nc-1,:)]/2
endfunction
