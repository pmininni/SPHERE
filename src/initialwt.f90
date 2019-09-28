! Calculates the angular velocity and its derivative.  wparam0-9 are
! used as input parameters. The injt parameter tells at which stat the
! precession starts working. wx0, wy0, wz0 are the intial conditions

      ! Constant magnitude angular velocity, rotating around x axis
      ! Initial condition should be wx=0,wy=0,wz!=0

      if (o.eq.ord) then
         tmp = (t-INT(injt*bstep)-1)*dt
      else
         tmp = (t-INT(injt*bstep)-1)*dt + dt/(o+1)
      endif
      wx  = 0
      wy  = wz0*SIN(wparam0*tmp)
      wz  = wz0*COS(wparam0*tmp)
      wtx = 0
      wty =  wparam0*wz0*COS(wparam0*tmp)
      wtz = -wparam0*wz0*SIN(wparam0*tmp)
