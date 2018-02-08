function [z]=beulerfun(w)
  global glob_h glob_t glob_y glob_odefun;
  z=w-glob_y-glob_h*feval(glob_odefun,glob_t,w);
end